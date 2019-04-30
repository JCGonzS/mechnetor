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


# def get_protein_ids(prot_data_file):
#     D = defaultdict(dict)
#     D["doms"] = defaultdict(dict)
#     with open_file(prot_data_file) as f:
#         for line in f:
#             t = line.rstrip().split("\t")
#             uni_id, uni_ac, gn, length = t[:4]
#             length = int(length)
#
#             D["AC"][uni_ac], D["AC"][uni_id], D["AC"][gn] = uni_ac, uni_ac, uni_ac
#             D["AC"][uni_ac.upper()], D["AC"][uni_id.upper()], D["AC"][gn.upper()] = uni_ac, uni_ac, uni_ac
#             D["ID"][uni_ac], D["ID"][uni_id], D["ID"][gn] = uni_id.upper(), uni_id.upper(), uni_id.upper()
#             D["ID"][uni_ac.upper()], D["ID"][uni_id.upper()], D["ID"][gn.upper()] = uni_id.upper(), uni_id.upper(), uni_id.upper()
#             D["GN"][uni_ac], D["GN"][uni_id], D["GN"][gn] = gn.upper(), gn.upper(), gn.upper()
#             D["GN"][uni_ac.upper()], D["GN"][uni_id.upper()], D["GN"][gn.upper()] = gn.upper(), gn.upper(), gn.upper()
#
#     return D

def get_protein_ids(prot_data):
    D = defaultdict(dict)
    for c in prot_data.find({}, {"_id":0, "uniprot_acc":1, "uniprot_id":1,
                            "gene":1, "data_class":1}):
        uni_ac = c["uniprot_acc"]
        uni_id = c["uniprot_id"]
        gn = c["gene"]
        dc = c["data_class"]
        if dc == "Reviewed":
            # print uni_ac, uni_id, gn
            # if uni_ac in D["AC"]:
            #     print ">",uni_ac
            # if uni_id in D["AC"]:
            #     print ">",uni_id
            # if gn in D["AC"]:
            #     print ">",gn
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

def main(client, query_prots, query_muts, mode, sps="Hsa", max_prots="",
        output_dir=""):

    if hasNumbers(max_prots):
        max_prots = int(re.search("(\d+)", max_prots).group(1))
    else:
        max_prots = 5

    ## DATA: Set MongoDB databases & collections ( client[database][collection] )
    protein_data = client['protein_data'][sps]
    cosmic_data = client['cosmicv87']['genome_screens']
    biogrid_data = client['interactions_'+sps]['biogrid_'+sps]
    iprets_data = client['interactions_'+sps]['iprets_'+sps]
    # dom_prop_data = client['interactions_'+sps]['domain_propensities_'+sps]
    dd_ass_data = client['interactions_'+sps]['dom_dom_ass_'+sps]
    db3did_data = client['interactions_common']['db3did']
    elm_int_data = client['interactions_common']['elm_dom']
    elm_classes = client['other_data']['elm_classes']

    st = datetime.datetime.now()
    print "\n[{}] Running PIV \"{}\"".format(st, sys.argv[0])

    ## Get protein dictionary & sequences for species
    prot_ids = get_protein_ids(protein_data)

    # 2. Query proteins
    input_proteins, custom_pairs = parse_input(query_prots, prot_ids, max_prots,
                                                                biogrid_data)
    print "\n[{}] Input read: {} proteins".format(st, len(input_proteins))
    sys.exit()
    if len(input_proteins) == 0:
        sys.exit("ERROR: no proteins were found in your input!")

    # 3. Query Mutations
    input_mutations = parse_mutation_input(query_muts, prot_ids, input_proteins)


    ## Run int2graph
    graph_ele, lines = int2graph.main(input_proteins, custom_pairs,
        protein_data, cosmic_data, input_mutations, biogrid_data, iprets_data,
        db3did_data, dd_ass_data, elm_int_data, elm_classes, max_prots)

    int_table = {}
    int_table["columns"] = ["#Gene(A)","Accession(A)","Gene(B)","Accession(B)",
        "Type", "F(A)","Start-End(A)", "Mutations(A)",
        "F(B)", "Start-End(B)", "Mutations(B)",
        "Info", "Source"]
    int_table["index"] = range(len(lines))
    int_table["data"] = lines

    ## Get output file names
    ide = "test"
    graph_json = output_dir+"graph_elements_"+ide+".json"
    table_json = output_dir+"interaction_table_"+ide+".json"
    table_tsv = output_dir+"interaction_table_"+ide+".tsv"

    if mode=="int":
        with open(table_tsv, "w") as output:
            output.write("\t".join(int_table["columns"])+"\n")
            for l in lines:
                output.write("\t".join(l)+"\n")
    else:
        ## Print graph as JSON file
        with open(graph_json, "w") as output:
            json.dump(graph_ele, output, indent=4)

        ## Print interactions as JSON file
        with open(table_json,"w") as output:
            json.dump(int_table, output)

    print "[{}] Job Completed".format(st)

    return

if __name__ == "__main__":
    # query_prots = "CREBBP ID2\nCREBBP MYOG\nMYOG ID2"
    # query_prots = "DLG5 SMARCA4\nSMARCA4 SMARCA2"
    # query_prots = "SMARCA4\nTP53"
    # query_muts = "BRCA2/T23C\nFANCA/S15P"
    query_prots = open_file(sys.argv[1]).read()
    if len(sys.argv) > 2:
        query_muts = open_file(sys.argv[2]).read()
    else:
        query_muts = ""
    client = MongoClient('localhost', 27017)
    main(client, query_prots, query_muts, mode="ints", sps="Hsa", max_prots="5",
            output_dir="static/examples/")
