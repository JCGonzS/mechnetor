#!/usr/bin/env python

""" PIV main script

JC Gonzalez Sanchez, 2018
"""


import os, sys, re
import gzip, json, random, string, datetime, pymongo
import pandas as pd
import run_interprets, int2graph
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

def make_protein_id_dictionary(prot_data):
    D = defaultdict(lambda: defaultdict(list))
    for doc in prot_data.find({},
                {"_id":0, "uni_ac":1, "uni_id":1, "genes":1, "alt_ids":1,
                "data_class":1}).sort("data_class", pymongo.ASCENDING):
                                        # First "Reviewed", then "Unreviewed"
        uni_ac = doc["uni_ac"]
        uni_id = doc["uni_id"]
        genes = doc["genes"]
        alt_ids = doc["alt_ids"]

        for dic, val in zip(["AC", "ID"], [uni_ac, uni_id]):
            for key in [uni_ac, uni_id]+genes+alt_ids:
                for k in [key, key.upper()]:
                    if val not in D[dic][k]:
                        D[dic][k].append(val)

        D["GN"][uni_ac] = genes

    return D

# def get_unique_uni_id(protein, uni_ids):
#
#     for uni_id in sorted(list(uni_ids), reverse=True):
#         if protein in uni_id:
#             break
#
#     return uni_id

def parse_input(input_text, prot_dict, max_prots, protein_data):
    input_prots = set()
    input_seqs = defaultdict(str)
    all_proteins = set()
    not_found = set()
    custom_pairs = []
    flag = 0

    for line in input_text.split("\n"):

        ## Make proteins UniProtID centric
        if line.strip() and line[0] != "#":
            if line[0] == ">":
                header = str(line.rstrip().split()[0].split(">")[1])
                header = header.replace("/","-")
                flag = 1
            elif (flag == 1 and len(line.split()) == 1
                  and not re.search("\d",line)):
                input_seqs[header]+=line.rstrip()
            else:
                flag = 0
                vals = line.rstrip().upper().split()
                uni_acs = set()
                for val in vals:
                    if val in prot_dict["AC"]:
                        uni_ac = prot_dict["AC"][val][0]
                        uni_acs.add(uni_ac)
                    else:
                        not_found.add("'"+val+"'")

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
                    gene = prot_dict["GN"][uni_ac][0]

                    if max_prots > 0:
                        ## Get X interactors
                        doc = protein_data.find_one({"uni_ac": uni_ac},
                                                {"_id":0, "biogrid_interactors":1})
                        interactors = []
                        for interactor_ac in doc["biogrid_interactors"]:
                            if (len(interactors) < max_prots):
                                interactors.append(interactor_ac)

                        all_proteins = all_proteins | set(interactors)
        else:
            flag = 0

    remove_keys = []
    if len(input_seqs) > 0:
        for key in input_seqs:
            if len(input_seqs[key]) == 0:
                remove_keys.append(key)
    for key in remove_keys:
        input_seqs.pop(key, "None")

    return input_prots, input_seqs, all_proteins, custom_pairs, not_found

def parse_mutation_input(input_text, prot_dict, protein_set, protein_data):
    mutations = defaultdict(lambda: defaultdict(set))
    for line in input_text.split("\n"):
        if line.strip() and line[0] != "#":
            prot, mut = line.rstrip().split()[0].split("/")
            pos = int(re.search("([\d]+)", mut).group(1))
            if prot in prot_dict["AC"]:
                uni_ac = prot_dict["AC"][prot][0]
            	if uni_ac in protein_set:
                    mutations[uni_ac][pos].add(mut)

    return mutations

def get_unique_random_identifier(output_dir):
    flag = 0
    while flag == 0:
        ide = ''.join(random.choice(
                string.ascii_uppercase + string.ascii_lowercase + string.digits
                ) for _ in range(8))
        outfile = output_dir+"interaction_table_"+ide+".tsv"
        if not os.path.isfile(outfile):
          flag = 1
    return ide

def get_pfam_info(data):
    d = {}
    for c in data.find():
        d[c["Accession"]] = (c["Identifier"], c["Description"])
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

def parse_blast(blastout_file, mask, data_dict, prot_id,
                max_ol=0.25, max_cov=0.9):
    prots = set()
    sp, tr = defaultdict(list), defaultdict(list)
    with open_file(blastout_file) as f:
        for line in f:
            if line[0]!="#":
                t = line.rstrip().split("\t")
                query = t[0]
                dc, subj = t[1].split("|")[:2]
                prots.add(query)
                if dc=="sp":
                    sp[query].append({ "acc": subj,
                                  "ide": float(t[2]),
                                  "e-val": float(t[10]),
                                  "q-start": int(t[6]),
                                  "q-end": int(t[7])
                                  })
                elif dc=="tr":
                    tr[query].append({ "acc": subj,
                                  "ide": float(t[2]),
                                  "e-val": float(t[10]),
                                  "q-start": int(t[6]),
                                  "q-end": int(t[7])
                                  })

    for prot in prots:
        mask_cov = 0
        final = []
        for hit in sp[prot]+tr[prot]:
            ol = calculate_overlap(hit["q-start"], hit["q-end"], mask[prot])
            if ol < max_ol:
                mask[prot] = fill_mask(hit["q-start"], hit["q-end"], mask[prot])
                final.append(hit)
            mask_cov = mask[prot].count("1")/float(len(mask[prot]))
            if mask_cov >= max_cov:
                break
        for hit in sorted(final, key=lambda i: i["q-start"]):
            gene = prot_id["GN"][hit["acc"]][0]
            l = str(hit["q-start"])+"-"+str(hit["q-end"])+"|"
            l += hit["acc"]+"|"+gene+"|("+"{:1.0e}".format(hit["e-val"])+", "+str(hit["ide"])+"%)"
            data_dict[prot]["blast"].append(l)
    return data_dict

def parse_pfamscan(pfam_file, data_dict):
    # pfams = defaultdict(lambda: defaultdict(set) )
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
        com = "{} {} | perl {} \"{}\" -m 1 > {}".format(mode, fasta_file, fpm2_script, elm_regex, tmp_file)
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

def calculate_overlap(start, end, mask):
    sub_mask = mask[start-1:end]
    overlap = sub_mask.count("1") / float(len(sub_mask))
    return overlap

def fill_mask(start, end, mask):
    length = end - start + 1
    for i in range(start-1, end):
        mask[i] = "1"
    return mask

@line_profile
def main(CLIENT, query_prots, query_muts, SP="Hsa",
         make_graph=True, hide_no_int=True, MAX_PROTS="",
         MAIN_OUTPUT_DIR="", temp_dir="temp/", blastdb_dir="",
         error_template="templates/input_error.html.jinja2",
         results_template="templates/results.html.jinja2"):

    LOG_FILE = MAIN_OUTPUT_DIR+"log.txt"
    sys.stdout = open(LOG_FILE, 'a')

    if make_graph == False:
        MAX_PROTS = 9999
    elif hasNumbers(MAX_PROTS):
        MAX_PROTS = int(re.search("(\d+)", MAX_PROTS).group(1))
    	# if MAX_PROTS > 50:
        #     MAX_PROTS = 50
    elif MAX_PROTS == "no":
        MAX_PROTS = 999
    else:
        MAX_PROTS = 5

    ## DATA: Set MongoDB databases & collections (CLIENT[database][collection])
    PFAM_DATA       = CLIENT["common"]["pfamA_data"]
    DB3DID_DATA     = CLIENT["common"]["db_3did"]
    ELM_INT_DATA    = CLIENT["common"]["elm_int_dom"]
    ELM_CLASSES     = CLIENT["common"]["elm_classes"]
    PROTEIN_DATA    = CLIENT[SP]["protein_data"]
    # BIOGRID_DATA    = CLIENT[SP]["biogrid"]
    PPI_DATA        = CLIENT[SP]["ppi_db"]
    ASS_PROB_DATA   = CLIENT[SP]["association_probabilities"]
    COSMIC_DATA, IPRETS_DATA = "no", "no"
    # if SP in ["Hsa", "Dme"]:
    IPRETS_DATA = CLIENT[SP]["interprets"]
    if SP == "Hsa":
        COSMIC_DATA = CLIENT["cosmic_v87"]["genome_screens"]
    blastdb_uni = blastdb_dir+"uniprot_"+SP
    blastdb_pdb = blastdb_dir+"pdbaa_2019"

    ## Get output file names
    IDE = get_unique_random_identifier(MAIN_OUTPUT_DIR)
    output_dir = MAIN_OUTPUT_DIR+"job_"+IDE+"/"
    if not os.path.exists(output_dir):
        try:
            os.mkdir(output_dir)
        except OSError:
            st = "[{}]".format(datetime.datetime.now())
            print st, "Creation of the directory {} failed".format(output_dir)
    # if not os.path.exists(temp_dir):
    #     st = "[{}]".format(datetime.datetime.now())
    #     print st, "Error. The specified 'temp' directory"+
    #           "'{}' does not exist'".format(temp_dir)

    fasta_file = output_dir+"seqs_"+IDE+".fasta"
    blastout_file = output_dir+"blastout_"+IDE+".tsv.gz"
    pfamout_file = output_dir+"pfamscan_"+IDE
    iprets_file = output_dir+"i2_"+IDE+".tsv.gz"
    graph_json = "graph_elements_"+IDE+".json"
    table_json = "interaction_table_"+IDE+".json"
    table_tsv  = "interaction_table_"+IDE+".tsv.gz"

    st = "[{}]".format(datetime.datetime.now())
    print st, "Running PIV. Job identifier: {}".format(IDE)

    pfam_info = get_pfam_info(PFAM_DATA)
    elm_info = get_elm_info(ELM_CLASSES)

    ### 1. Get protein dictionary & sequences for species
    prot_ids = make_protein_id_dictionary(PROTEIN_DATA)

    ### 2. Parse protein input
    (input_prots, input_seqs, all_prots, custom_pairs,
        not_found) = parse_input(query_prots, prot_ids, MAX_PROTS,
                                 PROTEIN_DATA)

    # print "input prots", input_prots
    # print all_prots

    if len(not_found) > 0:
        st = "[{}]".format(datetime.datetime.now())
        print st, "The following input protein(s) could not be identified:", "; ".join(not_found)

    total_n_prots = len(input_prots)+len(input_seqs.keys())
    if total_n_prots == 0:
        st = "[{}]".format(datetime.datetime.now())
        print st, "ERROR: no valid proteins found in input!"
        return render_template(error_template)

    ### 3. If sequence(s) in input -> compute data:
    fasta_data = defaultdict(lambda: defaultdict(list))
    fasta_link = {}
    fasta_iprets = {}

    if len(input_seqs) > 0:
        mask = defaultdict(list)

        ### 3-1. Print sequences in file
        with open_file(fasta_file, "w") as out:
            for head in input_seqs:
                fasta_link[head] = "static/jobs/job_"+IDE+"/seqs/"+head+".fa"
                mask[head] = ["0"] * len(input_seqs[head])
                out.write(">"+head+"\n")
                out.write(input_seqs[head]+"\n")

        ### 3-2. Run BLAST to identify similar proteins to input sequence(s)
        ## Blast parameters:
        psiblast_path = "psiblast"
        e_val = "0.01" #Default = 10
        ite = "1"  #-num_iterations (Default = 1)
        outfmt = "7"
        ## Run blast through bash.
        command = "{} -query {} -db {} -evalue {} -outfmt {} | gzip > {}".format(psiblast_path, fasta_file, blastdb_uni, e_val, outfmt, blastout_file)
        st = "[{}]".format(datetime.datetime.now())
        print st, "Running PSIBLAST: {}".format(command)
        os.system(command)
        ## Parse results
        fasta_data = parse_blast(blastout_file, mask, fasta_data, prot_ids)

        ### 3-3. Run PfamScan to identify Pfam domains on input sequence(s).
        # pfamscan.py parameters:
        pfamscan_path = "/var/www/flask_apps/jc_test/jc_app/pfamscan.py"
        evalue = "0.1"
        email = "juan-carlos.gonzalez@bioquant.uni-heidelberg.de"
        command = "python {} --sequence {} --database pfam-a --evalue {} --format txt --outfile {} --email {} --quiet".format(pfamscan_path, fasta_file, evalue, pfamout_file, email)
        st = "[{}]".format(datetime.datetime.now())
        print st, "Running PfamScan: {}".format(command)
        os.system(command)
        os.system("gzip "+pfamout_file+".out.txt")
        os.unlink(pfamout_file+".sequence.txt")
        # os.system("gzip "+pfamout_file+".sequence.txt")
        # Parse results
        fasta_data = parse_pfamscan(pfamout_file+".out.txt.gz", fasta_data)

        # 3-4. Run ELM search on input sequence(s).
        fpm2_path = "/var/www/flask_apps/jc_test/jc_app/fpm2.pl"
        print "[{}]".format(datetime.datetime.now()), "Searching for ELMs"
        fasta_data = find_all_elms(fasta_file, fasta_data, IDE, elm_info,
                                   fpm2_script=fpm2_path)

        # 3-5. Run InterPreTS
        all_seqs = input_seqs.copy()
        db_seqs = {}
        for uni_ac in input_prots:
            db_seqs[uni_ac] = PROTEIN_DATA.find_one({"uni_ac": uni_ac},
                                         {"_id": 0, "sequence": 1})["sequence"]
        all_seqs.update(db_seqs)
        print "[{}]".format(datetime.datetime.now()), "Running InterPrets"
        fasta_iprets = run_interprets.main(
                   all_seqs, output_dir, iprets_file, IDE,
                   mode="psiblast", psiblast=psiblast_path, blastdb=blastdb_pdb,
                   print_output=True, temp_dir=temp_dir)

    st = "[{}]".format(datetime.datetime.now())
    print st, "Your input contains {} proteins.".format(total_n_prots)
    # " Plus a maximum of {} interactors for each, the total is {} proteins and {} interactions".format(MAX_PROTS, len(all_prots), tot_ints)

    # Query Mutations
    input_muts = parse_mutation_input(query_muts, prot_ids, input_prots,
                                      PROTEIN_DATA)

    ## Run int2graph
    graph_ele, lines, no_int_prots = int2graph.main(SP, all_prots,
            input_prots, custom_pairs, input_seqs, input_muts,
            fasta_data, fasta_link, PROTEIN_DATA, COSMIC_DATA, PPI_DATA,
            IPRETS_DATA, fasta_iprets, DB3DID_DATA, ASS_PROB_DATA, ELM_INT_DATA,
            pfam_info, elm_info,
            make_graph=make_graph, hide_no_int=hide_no_int)

    if len(no_int_prots) > 0:
        st = "[{}]".format(datetime.datetime.now())
        print st,"No interactions found for: {}".format("; ".join(no_int_prots))

    st = "[{}]".format(datetime.datetime.now())
    print st, "int2graph.py run successfully"

    int_table = {}
    int_table["columns"] = ["#Gene(A)","Accession(A)","Gene(B)","Accession(B)",
        "Type", "F(A)","Start-End(A)", "Mutations(A)",
        "F(B)", "Start-End(B)", "Mutations(B)",
        "Info", "Source"]
    int_table["index"] = range(len(lines))
    int_table["data"] = lines

    ## Print graph as JSON file
    with open(output_dir+graph_json, "w") as output:
        json.dump(graph_ele, output, indent=4)

    ## Print interactions as JSON file
    with open(output_dir+table_json,"w") as output:
        json.dump(int_table, output)

    st = "[{}]".format(datetime.datetime.now())
    print st, "...done!. Created files \"{}\" and \"{}\"".format(graph_json,
                                                                 table_json)
    ## Print HTML
    sys.stdout = sys.__stdout__

    return render_template(results_template,
                           graph_json="jobs/"+"job_"+IDE+"/"+graph_json,
                           ints_json="jobs/"+"job_"+IDE+"/"+table_json,
                           not_identified=not_found,
                           no_int_prots=no_int_prots)
