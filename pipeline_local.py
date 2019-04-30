#!/usr/bin/env python

""" PIV main script

JC Gonzalez Sanchez, 2018
"""


import os, sys, re
import gzip, json, random, string
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

def get_protein_ids(prot_data):
    D = defaultdict(lambda: defaultdict(set))

    for c in prot_data.find({}, {"_id":0, "uniprot_acc":1, "uniprot_id":1,
                            "gene":1, "data_class":1}):
        uni_ac = c["uniprot_acc"]
        uni_id = c["uniprot_id"]
        gn = c["gene"]
        dc = c["data_class"]

        # if dc == "Reviewed":
        for dic, val in zip(["AC", "ID", "GN"], [uni_ac, uni_id, gn]):
            for key in [uni_ac, uni_id, gn]:
                D[dic][key].add(val.upper())
                D[dic][key.upper()].add(val.upper())
    return D

def get_all_BioGrid_ints(data, gene):
    ints = defaultdict(set)
    for cursor in data.find( {"$or":
                                [{"Official Symbol Interactor A": gene},
                                {"Official Symbol Interactor B": gene}]
                             },
                             {"_id": 0,
                             "Pubmed ID": 1,
                             "Official Symbol Interactor A": 1,
                             "Official Symbol Interactor B": 1}):

        if gene == cursor["Official Symbol Interactor A"]:
            interactor = cursor["Official Symbol Interactor B"]
        else:
            interactor = cursor["Official Symbol Interactor A"]
        ints[interactor].add(cursor["Pubmed ID"])

    return ints

def get_unique_uni_id(protein, uni_ids):

    for uni_id in sorted(list(uni_ids), reverse=True):
        if protein in uni_id:
            break

    return uni_id

def parse_input(input_text, prot_dict, max_prots, protein_data, biogrid_data):
    input_prots = set()
    all_proteins = set()
    not_found = set()
    multihits = {}
    custom_pairs = []
    tot_ints = 0
    for line in input_text.split("\n"):

        ## Make proteins UniProtID centric
        if line.strip() and line[0] != "#":
            vals = line.rstrip().upper().split()
            uni_acs = set()
            for v in vals:
                v = v.upper()
                if v in prot_dict["ID"]:
                    uni_id = get_unique_uni_id(v, prot_dict["ID"][v])
                    uni_ac = protein_data.find_one({"uniprot_id": uni_id},
                                    {"_id": 0, "uniprot_acc": 1})["uniprot_acc"]
                    uni_acs.add(uni_ac)

                    if len(prot_dict["ID"][v])>1:
                        multihits[v] = (uni_id, uni_ac)

                else:
                    not_found.add("'"+v+"'")

            input_prots = input_prots|uni_acs

            ## Types of input:
            ## 1. Pair of proteins
            if len(uni_acs) == 2:
                custom_pairs.append( list(uni_acs) )
                for uni_ac in uni_acs:
                    all_proteins.add(uni_ac)

            ## 2. Single protein
            elif len(uni_acs)>0:
                uni_ac = list(uni_acs)[0]
                all_proteins.add(uni_ac)

                ## Get X interactors
                interactors = []
                gene = protein_data.find_one({"uniprot_acc": uni_ac},
                                               {"_id": 0, "gene": 1})["gene"]
                ints = get_all_BioGrid_ints(biogrid_data, gene)
                tot_ints += len(ints)
                for int_gene in sorted(ints, key=lambda int: len(ints[int]),
                                                                reverse=True):
                    int_gene = int_gene.upper()
                    if (len(interactors) < max_prots and int_gene in prot_dict["ID"]):
                        uni_id = get_unique_uni_id(int_gene, prot_dict["ID"][int_gene])
                        uni_ac = protein_data.find_one({"uniprot_id": uni_id},
                                        {"_id": 0, "uniprot_acc": 1})["uniprot_acc"]
                        interactors.append(uni_ac)
                all_proteins = all_proteins | set(interactors)

    return input_prots, all_proteins, custom_pairs, not_found, tot_ints

def parse_mutation_input(input_text, prot_dict, protein_set, protein_data):
    mutations = defaultdict(lambda: defaultdict(set))
    for line in input_text.split("\n"):
        if line.strip() and line[0] != "#":
            prot, mut = line.rstrip().split()[0].split("/")
            pos = int(re.search("([\d]+)", mut).group(1))
            if prot in prot_dict["ID"]:
                uni_id = get_unique_uni_id(prot, prot_dict["ID"][prot])
                uni_ac = protein_data.find_one({"uniprot_id": uni_id},
                                {"_id": 0, "uniprot_acc": 1})["uniprot_acc"]
            if uni_ac in protein_set:
                mutations[uni_ac][pos].add(mut)

    return mutations

def main(client, query_prots, query_muts, mode, show,
            sps="Hsa", max_prots="", output_dir=""):

    if mode=="ints":
        max_prots = 9999
    elif hasNumbers(max_prots):
        max_prots = int(re.search("(\d+)", max_prots).group(1))
    else:
        max_prots = 5

    ## DATA: Set MongoDB databases & collections ( client[database][collection] )
    protein_data = client['protein_data'][sps+"_2019_allIntELMs"]
    cosmic_data = client['cosmicv87']['genome_screens']
    biogrid_data = client['interactions_'+sps]['biogrid_'+sps]
    iprets_data = client['interactions_'+sps]['iprets_'+sps]
    dd_ass_data = client['interactions_'+sps]['dom_dom_ass_'+sps]
    db3did_data = client['interactions_common']['db3did']
    elm_int_data = client['interactions_common']['elm_dom']
    elm_classes = client['other_data']['elm_classes']

    st = datetime.datetime.now()
    print "\n[{}] Running PIV".format(st)

    ## Get protein dictionary & sequences for species
    prot_ids = get_protein_ids(protein_data)

    # 2. Query proteins
    input_prots, all_prots, custom_pairs, not_found, tot_ints = parse_input(query_prots,
                                                        prot_ids, max_prots,
                                                        protein_data,
                                                        biogrid_data)
    if len(input_prots) == 0:
        sys.exit("ERROR: no proteins were found in your input!")

    print "[{}] Your input contains {} proteins. Plus a maximum of {} interactors for each, the total is {} proteins and {} interactions".format(st, len(input_prots), max_prots, len(all_prots), tot_ints)
    if len(not_found)>0:
        print "[{}] The following input protein(s) could not be identified:".format(st),"; ".join(not_found)


    # 3. Query Mutations
    input_muts = parse_mutation_input(query_muts, prot_ids, input_prots,
                                            protein_data)

    ## Run int2graph
    graph_ele, lines = int2graph.main(mode,
                    all_prots, input_prots, custom_pairs, input_muts,
                    protein_data, cosmic_data,  biogrid_data, iprets_data,
                    db3did_data, dd_ass_data, elm_int_data, elm_classes)

    int_table = {}
    int_table["columns"] = ["#Gene(A)","Accession(A)","Gene(B)","Accession(B)",
        "Type", "F(A)","Start-End(A)", "Mutations(A)",
        "F(B)", "Start-End(B)", "Mutations(B)",
        "Info", "Source"]
    int_table["index"] = range(len(lines))
    int_table["data"] = lines

    ## Get output file names
    ide = ''.join(random.choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) for _ in range(8))
    graph_json = output_dir+"graph_elements_"+ide+".json"
    table_json = output_dir+"interaction_table_"+ide+".json"
    table_tsv = output_dir+"interaction_table_"+ide+".tsv"

    if mode=="ints":
        if show:
            print "\t".join(int_table["columns"])
            for l in int_table["data"]:
                print "\t".join(l)
        elif:
            with open(table_tsv, "w") as output:
                output.write("\t".join(int_table["columns"])+"\n")
                for l in int_table["data"]:
                    output.write("\t".join(l)+"\n")
        print "[{}] Job Completed. Find your results in: {}".format(st,
                                                                    table_tsv)

    else:
        ## Print graph as JSON file
        with open(graph_json, "w") as output:
            json.dump(graph_ele, output, indent=4)

        ## Print interactions as JSON file
        with open(table_json,"w") as output:
            json.dump(int_table, output)

    return

if __name__ == "__main__":
    args = sys.argv[1:]
    show = False
    if "-show" in args:
        show = True
        args.remove("-show")

    query_prots, query_muts = "", ""
    if len(args) > 1:
        query_prots = open_file(args[1]).read()
    if len(args) > 2:
        query_muts = open_file(args[2]).read()

    # query_prots = "CREBBP ID2\nCREBBP MYOG\nMYOG ID2"
    # query_prots = "DLG5 SMARCA4\nSMARCA4 SMARCA2"
    # query_prots = "SMARCA4\nTP53"
    # query_prots = "ID2  CSK"
    # query_muts = "Q02363/T32L\nP41240/R16V"
    client = MongoClient('localhost', 27017)
    main(client, query_prots, query_muts, mode="ints", show=show,
        sps="Hsa", max_prots="5", output_dir="static/output/")
