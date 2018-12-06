#!/usr/bin/env python

""" PIV main script

JC Gonzalez Sanchez, 2018
"""


import os, sys, re
import gzip, json
import pandas as pd
import datetime
import int2graph
from collections import defaultdict
from pymongo import MongoClient
from flask import render_template, url_for
from flask_debugtoolbar_lineprofilerpanel.profile import line_profile


def open_file(input_file, mode="r"):
    """ Open file Zipped or not
    """

    if re.search(".gz$", input_file):
        infile = gzip.open(input_file, mode)
    else:
        infile = open(input_file, mode)
    return infile

def hasNumbers(inputString):
    return any(char.isdigit() for char in inputString)


def get_protein_ids(prot_data_file):
    D = defaultdict(dict)
    D["doms"] = defaultdict(dict)
    with open_file(prot_data_file) as f:
        for line in f:
            t = line.rstrip().split("\t")
            uni_id, uni_ac, gn, length = t[:4]
            length = int(length)

            D["AC"][uni_ac], D["AC"][uni_id], D["AC"][gn] = uni_ac, uni_ac, uni_ac
            D["AC"][uni_ac.upper()], D["AC"][uni_id.upper()], D["AC"][gn.upper()] = uni_ac, uni_ac, uni_ac
            D["ID"][uni_ac], D["ID"][uni_id], D["ID"][gn] = uni_id.upper(), uni_id.upper(), uni_id.upper()
            D["ID"][uni_ac.upper()], D["ID"][uni_id.upper()], D["ID"][gn.upper()] = uni_id.upper(), uni_id.upper(), uni_id.upper()
            D["GN"][uni_ac], D["GN"][uni_id], D["GN"][gn] = gn.upper(), gn.upper(), gn.upper()
            D["GN"][uni_ac.upper()], D["GN"][uni_id.upper()], D["GN"][gn.upper()] = gn.upper(), gn.upper(), gn.upper()

    return D

def get_all_BioGrid_ints(data, gene):
    ints = defaultdict(set)
    for cursor in data.find( {"$or":
                                [{"Official Symbol Interactor A": gene},
                                {"Official Symbol Interactor B": gene}]},
                                {"_id": 0, "#BioGRID Interaction ID": 1,
                                "Pubmed ID": 1,
                                "Official Symbol Interactor A": 1,
                                "Official Symbol Interactor B": 1}):

        if gene == cursor["Official Symbol Interactor A"]:
            interactor = cursor["Official Symbol Interactor B"]
        else:
            interactor = cursor["Official Symbol Interactor A"]

        ints[interactor].add(cursor["Pubmed ID"])

    return ints

def parse_input(input_text, prot_dict, max_prots, biogrid_data):
    protein_set = set()
    custom_pairs = []
    for line in input_text.split("\n"):

        if line.strip() and line[0] != "#":
            vals = line.rstrip().upper().split()
            prots = []
            for v in vals:
                if v in prot_dict["AC"]:
                    prots.append(prot_dict["AC"][v])

            # Types of input:
            # 1. Single protein
            if len(prots) == 1:
                prot = prots[0]
                gene = prot_dict["GN"][prot]
                protein_set.add(prot)
                interactors = []
                ints = get_all_BioGrid_ints(biogrid_data, gene)
                for int in sorted(ints, key=lambda int: len(ints[int]),
                                                                reverse=True):
                    if len(interactors) < max_prots and int in prot_dict["AC"]:
                        interactors.append(prot_dict["AC"][int])
                protein_set = protein_set | set(interactors)

            # 2. Pair of proteins
            elif len(prots) == 2:
                custom_pairs.append( prots )
                for prot in prots:
                    protein_set.add(prot)

    return protein_set, custom_pairs

def parse_mutation_input(input_text, prot_dict, protein_set):
    mutations = defaultdict(lambda: defaultdict(set))
    for line in input_text.split("\n"):
        if line.strip() and line[0] != "#":
            prot, mut = line.rstrip().split()[0].split("/")
            pos = int(re.search("([\d]+)", mut).group(1))
            if prot in prot_dict["AC"]:
                prot = prot_dict["AC"][prot]
            if prot in protein_set:
                mutations[prot][pos].add(mut)
    return mutations

@line_profile
def main(client, query_prots, query_muts, max_prots="", query_lmd2="",
		 sps="Hsa", max_pval=999, main_dir=""):

    if hasNumbers(max_prots):
        max_prots = int(re.search("(\d+)", max_prots).group(1))
    else:
        max_prots = 5

    ## Set MongoDB databases & collections ( client[database][collection] )
    protein_data = client['protein_data'][sps]
    biogrid_data = client['interactions_'+sps]['biogrid_'+sps]
    iprets_data = client['interactions_'+sps]['iprets_'+sps]
    dom_prop_data = client['interactions_'+sps]['domain_propensities_'+sps]
    db3did_data = client['interactions_common']['db3did']
    elm_int_data = client['interactions_common']['elm_dom']
    elm_classes = client['other_data']['elm_classes']

    ## Data Directories & Files
    data_dir = main_dir+"static/data/"
    sps_dir = data_dir+"species/"+sps+"/"
    proteome_file = sps_dir+"uniprot_sprot_human_complete_20303prts_March_2018.fasta.gz"
    protein_data_file = sps_dir+"swissprot_human_parsed_data.tsv.gz"
    protein_data_json_file = sps_dir+"protein_data_new_Hsa_normal.json.gz"
    dd_prop_file = sps_dir + "dom_dom_lo.txt"
    output_dir = "static/output/"

    log_file = main_dir+"static/log.txt"
    sys.stdout = open(log_file, 'a')
    st = datetime.datetime.now()
    print "\n[{}] Running \"{}\"".format(st, "PIV pipeline.py")

    ## Get protein dictionary & sequences for species
    prot_ids = get_protein_ids(protein_data_file)

    ## Read Input/Query Values & Check Format
    # 1. Name of the query
    # if query_name == "":
    #     query_name = "Untitled Graph"

    # 2. Query proteins
    input_proteins, custom_pairs = parse_input(query_prots, prot_ids, max_prots,
                                                                biogrid_data)
    if len(input_proteins) == 0:
        print "ERROR: no proteins were found in input!"
        return render_template("input_error.html")

    ### One protein only-input option.
    # 1. Get all possible interactors for this protein in BioGRID.
    # 2. Sort them by evidence (number of publications)
    # 3. Use only "max_prots" of those.
    ########

    # 3. Query Mutations
    input_mutations = parse_mutation_input(query_muts, prot_ids, input_proteins)

    print "[{}] Received input from HTML form".format(st)

    # lmd2_file = ""
    # if query_lmd2 != "":
    #     tf_lmd2 = tempfile.NamedTemporaryFile(delete=False)
    #     tf_lmd2.write(query_lmd2)
    #     tf_lmd2.close()
    #     lmd2_file = tf_lmd2.name

    ## Get output file names
    for i in range(0, 500):
      number = "0"*(4-len(str(i)))+str(i)
      graph_json = "graph_elements"+number+".json"
      graph_path = main_dir+output_dir+graph_json
      ints_json = "interaction_table"+number+".json"
      ints_path = main_dir+output_dir+ints_json
      if not os.path.isfile(graph_path):
        break

    ## Run int2graph
	print "[{}] Running int2graph...".format(st)

    int2graph.main(input_proteins, custom_pairs, protein_data, input_mutations,
            biogrid_data, iprets_data, db3did_data, dom_prop_data, elm_int_data, elm_classes,
            max_prots, graph_path, ints_path)

    print "[{}] ...done!. Created files \"{}\" and \"{}\"".format(st, graph_json,
                                                                   ints_json)

    ## Print HTML
    print "[{}] Printing HTML output".format(st)
    sys.stdout = sys.__stdout__

    return render_template("results_page.html",
                           graph_json = "output/"+graph_json,
                           ints_json = "output/"+ints_json)
