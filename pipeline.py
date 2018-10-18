#!/usr/bin/env python

""" PIV main script

JC Gonzalez Sanchez, 2018
"""


import os, sys, re
import gzip, json
import pandas as pd
import datetime
import int2mech
import make_graph_elements
from collections import defaultdict
from flask import render_template, url_for
from flask_debugtoolbar_lineprofilerpanel.profile import line_profile
from pymongo import MongoClient

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

def parse_fasta(proteome_file):
    """Get protein sequences from FASTA file
    and the convertion between the protein IDs: accession, uniprot ID and gene
    """

    D = defaultdict(dict)
    seqs, masks = {}, {}

    with open_file(proteome_file) as f:
        for line in f:
            if line[0] == ">":
                uni_ac, uni_id = line.split()[0].split("|")[1:]
                gn = uni_id
                if "GN=" in line:
                    gn = re.search(r"GN=([^\s]*)\s",line).group(1)

                D["AC"][uni_ac], D["AC"][uni_id], D["AC"][gn] = uni_ac.upper(), uni_ac.upper(), uni_ac.upper()
                D["AC"][uni_ac.upper()], D["AC"][uni_id.upper()], D["AC"][gn.upper()] = uni_ac.upper(), uni_ac.upper(), uni_ac.upper()
                D["ID"][uni_ac], D["ID"][uni_id], D["ID"][gn] = uni_id.upper(), uni_id.upper(), uni_id.upper()
                D["ID"][uni_ac.upper()], D["ID"][uni_id.upper()], D["ID"][gn.upper()] = uni_id.upper(), uni_id.upper(), uni_id.upper()
                D["GN"][uni_ac], D["GN"][uni_id], D["GN"][gn] = gn.upper(), gn.upper(), gn.upper()
                D["GN"][uni_ac.upper()], D["GN"][uni_id.upper()], D["GN"][gn.upper()] = gn.upper(), gn.upper(), gn.upper()

                seqs[uni_ac.upper()] = ""
                # masks[uni_ac] = ""

            else:
                seqs[uni_ac.upper()] += line.rstrip()
                # masks[uni_ac] += "0" * len(line.rstrip())

    return D, seqs #, masks

def get_protein_data(prot_data_file):
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
            # D["len"][uni_ac] = length

            # if "domains" in line:
            #     for dom in t[-1].rstrip().split("//")[:-1]:
            #         coor,name = dom.split(":")
            #         start, end = coor.split("-")
            #         start, end = int(start), int(end)
            #         D["doms"][uni_ac][(start,end)] = name
    return D

def parse_input(input_text, prot_dict, max_prots):
    # interaction_pairs = set()
    protein_set = set()
    for line in input_text.split("\n"):

        if line.strip() and line[0] != "#":
            vals = line.rstrip().upper().split()

            for v in vals:
                if v in prot_dict["AC"]:
                    prot = prot_dict["AC"][v]

                    if len(protein_set) < max_prots:
                        protein_set.add(prot)

    # for prot_a in protein_set:
    #     for prot_b in protein_set:
    #         if prot_a != prot_b:
    #             if prot_a < prot_b:
    #                 interaction_pairs.add((prot_a,prot_b))
    #             else:
    #                 interaction_pairs.add((prot_b,prot_a))
    # return interaction_pairs, protein_set
    return protein_set

def parse_mutation_input(input_text, prot_dict, protein_set):
    mutations = defaultdict(lambda: defaultdict(set))
    for line in input_text.split("\n"):
        if line.strip() and line[0] != "#":
            prot, mut = line.rstrip().split("/")
            pos = int(re.search("([\d]+)", mut).group(1))
            if prot in prot_dict["AC"]:
                prot = prot_dict["AC"][prot]
            if prot in protein_set:
                mutations[prot][pos].add(mut)
    return mutations

@line_profile
def main(query_prots, query_muts, max_prots="", query_lmd2="",
		 sps="Hsa", max_pval=999, main_dir=""):

    client = MongoClient('localhost', 27017)

    if hasNumbers(max_prots):
        max_prots = int(re.search("(\d+)", max_prots).group(1))
    else:
        max_prots = 20

    ## Data Directories & Files
    data_dir = main_dir+"static/data/"
    sps_dir = data_dir+"species/"+sps+"/"
    proteome_file = sps_dir+"uniprot_sprot_human_complete_20303prts_March_2018.fasta.gz"
    protein_data_file = sps_dir+"swissprot_human_parsed_data.tsv.gz"
    protein_data_json_file = sps_dir+"protein_data_new_Hsa_normal.json.gz"
    output_dir = "static/output/"
    log_file = main_dir+"static/log.txt"

    sys.stdout = open(log_file, 'a')
    st = datetime.datetime.now()
    print "\n[{}] Running \"{}\"".format(st, "PIV pipeline.py")

    ## Get protein dictionary & sequences for species
    prot_ids = get_protein_data(protein_data_file)
    protein_data = json.load(open_file(protein_data_json_file))

    ## Read Input/Query Values & Check Format

    # 1. Name of the query
    # if query_name == "":
    #     query_name = "Untitled Graph"

    # 2. Query proteins
    input_proteins = parse_input(query_prots, prot_ids, max_prots)
    if len(input_proteins) == 0:
        sys.exit("ERROR: no proteins were found in your input!")
        #FIX: I could render a template for an error html where this message is printed.
        # better yet would be to check the input before submission (no idea how)

    example = 0
    if "#Example" in query_prots:
        example = 1

    # 3. Query Mutations
    input_mutations = parse_mutation_input(query_muts, prot_ids, input_proteins)

    print "[{}] Received input from HTML form".format(st)

    lmd2_file = ""
    # if query_lmd2 != "":
    #     tf_lmd2 = tempfile.NamedTemporaryFile(delete=False)
    #     tf_lmd2.write(query_lmd2)
    #     tf_lmd2.close()
    #     lmd2_file = tf_lmd2.name

    ## Get output file names
    for i in range(0, 500):
      number = "0"*(4-len(str(i)))+str(i)
      outfile_int = "interactions"+number+".tsv"
      outfile_json = "graph_elements"+number+".json"
      outfile_table_json = "interaction_table"+number+".json"
      if not os.path.isfile(main_dir+output_dir+outfile_int):
        break

    ## Run int2mech: create interaction data file
    if example == 0:
    	print "[{}] Running int2mech. Using {} as protein data...".format(st,protein_data_file )
        int_file =  main_dir+output_dir+outfile_int
    	int2mech.main(input_proteins, prot_ids, protein_data,
    	              int_file, data_dir, sps, max_prots, client)
    	print "[{}] ...interaction file created in \"{}\"".format(st, int_file)
        tsv_file = "output/" + outfile_int
    else:
        if "TCF3" in query_prots:
            int_file = main_dir + "static/examples/ints_example_TCF3.tsv"
            tsv_file = "examples/ints_example_TCF3.tsv"
        elif "SMARCA4" in query_prots:
            int_file = main_dir + "static/examples/ints_example_SMARCA4.tsv"
            tsv_file = "examples/ints_example_SMARCA4.tsv"
        elif "SORBS3" in query_prots:
            int_file = main_dir + "static/examples/ints_example_SORBS3.tsv"
            tsv_file = "examples/ints_example_SORBS3.tsv"

    ## Create cytoscape's graph elements
    print "[{}] Running graph creation tool...".format(st)
    make_graph_elements.main(input_proteins, protein_data, input_mutations,
             		int_file, lmd2_file="",
                    output_file= main_dir+output_dir+outfile_json,
             		max_pval=max_pval)
    print "[{}] ...graph created in \"{}\"".format(st,
     										main_dir+output_dir+outfile_json)

    # Create JSON file from TSV interaction file
    data = pd.read_csv(int_file, sep='\t', index_col=False)
    data.to_json(main_dir+output_dir+outfile_table_json, orient='split')

    ## Print HTML
    print "[{}] Printing HTML output".format(st)
    sys.stdout = sys.__stdout__

    return render_template("results_page.html",
                           elements_json = "output/" + outfile_json,
                           ints_json = "output/" + outfile_table_json)
