#!/usr/bin/env python

""" PIV main script

JC Gonzalez Sanchez, 2018
"""

import os, sys, re
import gzip, json, random, string
import datetime, argparse
import pandas as pd
import int2graph, run_interprets
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
    for c in prot_data.find({"data_class": "Reviewed"},
                            {"_id":0, "uni_ac":1, "uni_id":1, "gene":1
                            }):
        uni_ac = c["uni_ac"]
        uni_id = c["uni_id"]
        gn = c["gene"]

        for dic, val in zip(["AC", "ID", "GN"], [uni_ac, uni_id, gn]):
            for key in [uni_ac, uni_id, gn]:
                D[dic][key].add(val.upper())
                D[dic][key.upper()].add(val.upper())

    for c in prot_data.find({"data_class": "Unreviewed"},
                            {"_id":0, "uni_ac":1, "uni_id":1, "gene":1
                             }):
        uni_ac = c["uni_ac"]
        uni_id = c["uni_id"]
        gn = c["gene"]

        for dic, val in zip(["AC", "ID", "GN"], [uni_ac, uni_id, gn]):
            for key in [uni_ac, uni_id, gn]:
                if key not in D[dic]:
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
    input_seqs = defaultdict(str)
    all_proteins = set()
    not_found = set()
    multihits = {}
    custom_pairs = []
    tot_ints = 0
    flag = 0
    # for record in SeqIO.parse(input_text.split("\n"), "fasta"):
    #     print record.id

    for line in input_text.split("\n"):

        ## Make proteins UniProtID centric
        if line.strip() and line[0] != "#":

            if line[0]==">":
                header = str(line.rstrip().split()[0].split(">")[1])
                flag = 1
            elif (flag==1 and len(line.split())==1 and not re.search("\d",line)):
                input_seqs[header]+=line.rstrip()

            else:
                flag = 0
                vals = line.rstrip().upper().split()
                uni_acs = set()
                for v in vals:
                    v = v.upper()
                    if v in prot_dict["ID"]:
                        uni_id = get_unique_uni_id(v, prot_dict["ID"][v])
                        uni_ac = protein_data.find_one({"uni_id": uni_id},
                                        {"_id": 0, "uni_ac": 1})["uni_ac"]
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
                    gene = protein_data.find_one({"uni_ac": uni_ac},
                                                   {"_id": 0, "gene": 1})["gene"]
                    ints = get_all_BioGrid_ints(biogrid_data, gene)
                    tot_ints += len(ints)
                    for int_gene in sorted(ints, key=lambda int: len(ints[int]),
                                                                    reverse=True):
                        int_gene = int_gene.upper()
                        if (len(interactors) < max_prots and int_gene in prot_dict["ID"]):
                            uni_id = get_unique_uni_id(int_gene, prot_dict["ID"][int_gene])
                            uni_ac = protein_data.find_one({"uni_id": uni_id},
                                            {"_id": 0, "uni_ac": 1})["uni_ac"]
                            interactors.append(uni_ac)
                    all_proteins = all_proteins | set(interactors)
        else:
            flag = 0

    remove_keys = []
    if len(input_seqs)>0:
        for key in input_seqs:
            if len(input_seqs[key])==0:
                remove_keys.append(key)
    for key in remove_keys:
        input_seqs.pop(key, "None")

    return input_prots, input_seqs, all_proteins, custom_pairs, not_found, tot_ints

def parse_mutation_input(input_text, prot_dict, protein_set, protein_data):
    mutations = defaultdict(lambda: defaultdict(set))
    for line in input_text.split("\n"):
        if line.strip() and line[0] != "#":
            prot, mut = line.rstrip().split()[0].split("/")
            pos = int(re.search("([\d]+)", mut).group(1))
            if prot in prot_dict["ID"]:
                uni_id = get_unique_uni_id(prot, prot_dict["ID"][prot])
                uni_ac = protein_data.find_one({"uni_id": uni_id},
                                {"_id": 0, "uni_ac": 1})["uni_ac"]
                if uni_ac in protein_set:
                    mutations[uni_ac][pos].add(mut)

    return mutations

def get_unique_random_identifier(main_output_dir):
    flag = 0
    while flag == 0:
        ide = ''.join(random.choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) for _ in range(8))
        directory = main_output_dir+"job_"+ide+"/"
        if not os.path.exists(directory):
            flag = 1
    return ide

def get_pfam_info(data):
    d = {}
    for c in data.find():
        d[c["Acc"]] = (c["Ide"], c["Des"])
    return d

def get_elm_info(data):
    ## Collection fields:
    #    Accession       ELMIdentifier   FunctionalSiteName      Description
    #    Regex   Probability     #Instances      #Instances_in_PDB
    d = {}
    for c in data.find():
        d[c["Accession"]] = {
            "ide":   c["ELMIdentifier"],
            "name":  c["FunctionalSiteName"],
            "des":   c["Description"],
            "regex": c["Regex"],
            "prob":  c["Probability"]}
    return d

def parse_blast(blastout_file):
    hits = defaultdict(list)
    with open_file(blastout_file) as f:
        for line in f:
            if line[0]!="#":
                t = line.rstrip().split("\t")
                query = t[0]
                dc, subj = t[1].split("|")[:2]
                hits[query].append({ "acc": subj,
                              "ide": float(t[2]),
                              "e-val": float(t[10]),
                              "q-start": int(t[6]),
                              "q-end": int(t[7])
                })
    return hits

def parse_pfamscan(pfam_file, data_dict):
    with open_file(pfam_file) as f:
        for line in f:
            if line[0] != "#" and line.strip():
                t = line.rstrip().split()
                query = t[0]
                start, end = int(t[3]), int(t[4])
                pfam_ac, pfam_name = t[5].split(".")[0], t[6]
                e_val = float(t[12])
                data_dict[query]["pfams"].append({"acc": pfam_ac,
                                                  "start": int(start),
                                                  "end": int(end),
                                                  "e-val": float(e_val)
                                                 })
    return data_dict

def find_all_elms(fasta_file, data_dict, ide, elm_info, fpm2_script = "fpm2.pl"):

    tmp_file = "/tmp/elm_"+ide+".txt"
    for elm_acc in elm_info:
        elm_ide = elm_info[elm_acc]["ide"]
        elm_regex = elm_info[elm_acc]["regex"]
        elm_prob = elm_info[elm_acc]["prob"]

        mode = "cat"
        if ".gz" in fasta_file:
            mode = "zcat"
        com = "{} {} | ./{} \"{}\" -m 1 > {}".format(mode, fasta_file, fpm2_script, elm_regex, tmp_file)
        os.system( com )
        with open_file(tmp_file) as f:
            for line2 in f:
                if line2[0]==">":
                    label = line2.split()[0].split("/")[0].replace(">","")
                    start, end = line2.split()[0].split("/")[1].split("-")
                else:
                    seq = line2.rstrip()
                    if "elms" in data_dict[label]:
                        data_dict[label]["elms"][elm_acc].append(
                                                    {"start" : int(start),
                                                     "end" :   int(end),
                                                     "seq" :   seq
                                                     })
                    else:
                        data_dict[label]["elms"]=defaultdict(list)
    os.unlink(tmp_file)
    return data_dict

def main(client, query_prots, query_muts,
            sps="Hsa", max_prots="",
            blastdb_dir="", main_output_dir="",
            make_graph=True, show=False):

    if not make_graph:
        max_prots = 9999
    elif hasNumbers(max_prots):
        max_prots = int(re.search("(\d+)", max_prots).group(1))
    else:
        max_prots = 5

    ## DATA: Set MongoDB databases & collections ( client[database][collection] )
    db3did_data = client["common"]["3did"]
    elm_int_data = client["common"]["elm_int_dom"]
    elm_classes = client["common"]["elm_classes"]
    pfam_data = client["common"]["pfamA_data"]

    protein_data = client[sps]["protein_data"]
    biogrid_data = client[sps]["biogrid"]
    dd_ass_data = client[sps]["dom_dom_ass"]
    cosmic_data, iprets_data = "no", "no"
    blastdb = blastdb_dir+"uniprot_"+sps
    if sps in ["Hsa", "Dme"]:
        iprets_data = client[sps]["interprets_biogrid"]
    if sps == "Hsa":
        cosmic_data = client["cosmic_v87"]["genome_screens"]

    ## Get output file names
    ide = get_unique_random_identifier(main_output_dir)
    output_dir = main_output_dir+"job_"+ide+"/"
    if not os.path.exists(output_dir):
        try:
            os.mkdir(output_dir)
        except OSError:
            st = datetime.datetime.now()
            print "[{}] Creation of the directory {} failed".format(st, output_dir)
    fasta_file = output_dir+"seqs_"+ide+".fasta"
    blastout_file = output_dir+"blastout_"+ide+".tsv.gz"
    pfamout_file = output_dir+"pfamscan_"+ide
    graph_json = "graph_elements_"+ide+".json"
    table_json = "interaction_table_"+ide+".json"
    table_tsv  = "interaction_table_"+ide+".tsv.gz"

    st = datetime.datetime.now()
    print "\n[{}] Running PIV. Job identifier: {}".format(st, ide)

    pfam_info = get_pfam_info(pfam_data)
    elm_info = get_elm_info(elm_classes)

    # 1. Get protein dictionary & sequences for species
    prot_ids = get_protein_ids(protein_data)

    # 2. Parse protein input
    input_prots, input_seqs, all_prots, custom_pairs, not_found, tot_ints = parse_input(query_prots,
                                                                            prot_ids, max_prots,
                                                                            protein_data,
                                                                            biogrid_data)

    if len(not_found)>0:
        st = datetime.datetime.now()
        print "[{}] The following input protein(s) could not be identified:".format(st),"; ".join(not_found)

    total_n_prots = len(input_prots)+len(input_seqs.keys())
    if total_n_prots == 0:
        st = datetime.datetime.now()
        sys.exit("[{}] ERROR: no proteins were found in your input!".format(st))


    # 3. If sequence(s) in input -> compute data:
    fasta_data = defaultdict(lambda: defaultdict(list))
    fasta_iprets = {}
    if len(input_seqs)>0:

        # 3-1. Print sequences in file
        with open_file(fasta_file, "w") as out:
            for head in input_seqs:
                out.write(">"+head+"\n")
                out.write(input_seqs[head]+"\n")

        # 3-2. Run BLAST to identify similar proteins to input sequence(s).
        # # Blast parameters:
        # psiblast_path = "psiblast"
        # e_val = "0.01" #Default = 10
        # ite = "1"  #Default = 1
        # outfmt = "7"
        # # Run blast through bash.
        # command = "{} -query {} -db {} -evalue {} -num_iterations {} -outfmt {} | gzip > {}".format(psiblast_path, fasta_file, blastdb, e_val, ite, outfmt, blastout_file)
        # st = datetime.datetime.now()
        # print "[{}] Running PSIBLAST: {}".format(st, command)
        # os.system(command)
        # # Parse results
        # hits = parse_blast(blastout_file)


        # 3-3. Run PfamScan to identify Pfam domains on input sequence(s).
        # pfamscan.py parameters:
        pfamscan_path = "/var/www/flask_apps/jc_test/jc_app/pfamscan.py"
        evalue = "1.0"
        email = "juan-carlos.gonzalez@bioquant.uni-heidelberg.de"
        command = "python {} --sequence {} --database pfam-a --evalue {} --format txt --outfile {} --email {} --quiet".format(pfamscan_path, fasta_file, evalue, pfamout_file, email)
        st = datetime.datetime.now()
        print "[{}] Running PfamScan: {}".format(st, command)
        os.system(command)
        os.system("gzip "+pfamout_file+".out.txt")
        os.unlink(pfamout_file+".sequence.txt")
        # os.system("gzip "+pfamout_file+".sequence.txt")
        # Parse results
        fasta_data = parse_pfamscan(pfamout_file+".out.txt.gz", fasta_data)

        # 3-4. Run ELM search on input sequence(s).
        st = datetime.datetime.now()
        print "[{}] Searching for ELMs".format(st)
        fasta_data = find_all_elms(fasta_file, fasta_data, ide, elm_info)

        # 3-5. Run InterPreTS
        all_seqs=input_seqs.copy()
        db_seqs = {}
        for uni_ac in input_prots:
            db_seqs[uni_ac] = protein_data.find_one({"uni_ac": uni_ac},
                                                   {"_id": 0, "sequence": 1})["sequence"]
        all_seqs.update(db_seqs)
        st = datetime.datetime.now()
        print "[{}] Running InterPrets".format(st)
        fasta_iprets = run_interprets.main(all_seqs, output_dir, "", print_output=False)

    print "[{}] Your input contains {} proteins. Plus a maximum of {} interactors for each, the total is {} proteins and {} interactions".format(st, total_n_prots, max_prots, len(all_prots), tot_ints)

    # 4. Query Mutations
    input_muts = parse_mutation_input(query_muts, prot_ids, input_prots,
                                            protein_data)

    # 5. Run int2graph
    graph_ele, lines = int2graph.main(make_graph, sps, all_prots,
                    input_prots, custom_pairs, input_seqs, input_muts, fasta_data,
                    protein_data, cosmic_data, biogrid_data, iprets_data, fasta_iprets,
                    db3did_data, dd_ass_data, elm_int_data, pfam_info, elm_info )

    int_table = {}
    int_table["columns"] = ["#Gene(A)","Accession(A)","Gene(B)","Accession(B)",
        "Type", "F(A)","Start-End(A)", "Mutations(A)",
        "F(B)", "Start-End(B)", "Mutations(B)",
        "Info", "Source"]
    int_table["index"] = range(len(lines))
    int_table["data"] = lines


    if make_graph:
        ## Print graph as JSON file
        with open(output_dir+graph_json, "w") as output:
            json.dump(graph_ele, output, indent=4)

        ## Print interactions as JSON file
        with open(output_dir+table_json,"w") as output:
            json.dump(int_table, output)

    else:
        if show:
            print "\t".join(int_table["columns"])
            for l in int_table["data"]:
                print "\t".join(l)
        else:
            with open_file(output_dir+table_tsv, "w") as output:
                output.write("\t".join(int_table["columns"])+"\n")
                for l in int_table["data"]:
                    output.write("\t".join(l)+"\n")

            st = datetime.datetime.now()
            print "[{}] Job Completed. Find your results in: {}".format(st,
                                                        output_dir+table_tsv)


def true_or_false(x, arg):
    if arg[x]:
        return True
    else:
        return False

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-p", "--prots", required=True,
    	help="File with protein input (protein identifier or FASTA sequences)")
    ap.add_argument("-m", "--muts", required=False,
    	help="File with protein input (protein identifier or FASTA sequences)")
    ap.add_argument("-sp", "--species", required=False, default="Hsa",
        help="Organism (default: Hsa)")
    ap.add_argument("-show", required=False, action="store_true",
        help="Prints results on screen")
    args = vars(ap.parse_args())


    query_prots, query_muts = "", ""
    if args["prots"]:
        query_prots = open_file(args["prots"]).read()
    if args["muts"]:
        query_muts = open_file(args["muts"]).read()

    client = MongoClient('localhost', 27017)
    main(   client, query_prots, query_muts,
            sps=args["species"], max_prots="5",
            blastdb_dir = "static/data/blastdb/",
            main_output_dir="static/output/",
            make_graph=False, show=args["show"],)
    sys.exit()
