#!/usr/bin/env python

""" PIV main script

JC Gonzalez Sanchez, 2018
"""


import os, sys, re
import gzip, json, random, string, datetime
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
                header = header.replace("/","-")
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

def get_unique_random_identifier(output_dir):
    flag = 0
    while flag == 0:
        ide = ''.join(random.choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) for _ in range(8))
        outfile = output_dir+"interaction_table_"+ide+".tsv"
        if not os.path.isfile(outfile):
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

def parse_blast(blastout_file, mask, data_dict, prot_id, max_ol=0.25, max_cov=0.9):
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
            gene = list(prot_id["GN"][hit["acc"]])[0]
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
def main(client, query_prots, query_muts, make_graph,
         sps="Hsa", max_prots="",
         main_output_dir="", temp_dir="temp/", blastdb_dir=""):

    log_file = main_output_dir+"log.txt"
    sys.stdout = open(log_file, 'a')

    if make_graph==False:
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
    blastdb_uni = blastdb_dir+"uniprot_"+sps
    blastdb_pdb = blastdb_dir+"pdbaa_2019"
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
    # if not os.path.exists(temp_dir):
    #     st = datetime.datetime.now()
    #     print "[{}] Error. The specified 'temp' directory '{}' does not exist'".format(st, temp_dir)

    fasta_file = output_dir+"seqs_"+ide+".fasta"
    blastout_file = output_dir+"blastout_"+ide+".tsv.gz"
    pfamout_file = output_dir+"pfamscan_"+ide
    iprets_file = output_dir+"i2_"+ide+".tsv.gz"
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
        print "[{}] ERROR: no proteins found in input!".format(st)
        return render_template("input_error.html")

    # 3. If sequence(s) in input -> compute data:
    fasta_data = defaultdict(lambda: defaultdict(list))
    fasta_link = {}
    fasta_iprets = {}
    if len(input_seqs)>0:
        mask = defaultdict(list)
        # 3-1. Print sequences in file
        with open_file(fasta_file, "w") as out:
            for head in input_seqs:
                fasta_link[head] = "static/jobs/job_"+ide+"/seqs/"+head+".fa"
                mask[head] = ["0"]*len(input_seqs[head])
                out.write(">"+head+"\n")
                out.write(input_seqs[head]+"\n")

        # 3-2. Run BLAST to identify similar proteins to input sequence(s).
        # Blast parameters:
        psiblast_path = "psiblast"
        e_val = "0.01" #Default = 10
        ite = "1"  #-num_iterations (Default = 1)
        outfmt = "7"
        # Run blast through bash.
        command = "{} -query {} -db {} -evalue {} -outfmt {} | gzip > {}".format(psiblast_path, fasta_file, blastdb_uni, e_val, outfmt, blastout_file)
        st = datetime.datetime.now()
        print "[{}] Running PSIBLAST: {}".format(st, command)
        os.system(command)
        # Parse results
        fasta_data = parse_blast(blastout_file, mask, fasta_data, prot_ids)

        # 3-3. Run PfamScan to identify Pfam domains on input sequence(s).
        # pfamscan.py parameters:
        pfamscan_path = "/var/www/flask_apps/jc_test/jc_app/pfamscan.py"
        evalue = "0.1"
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
        fpm2_path = "/var/www/flask_apps/jc_test/jc_app/fpm2.pl"
        st = datetime.datetime.now()
        print "[{}] Searching for ELMs".format(st)
        fasta_data = find_all_elms(fasta_file, fasta_data, ide, elm_info, fpm2_script=fpm2_path)

        # 3-5. Run InterPreTS
        all_seqs=input_seqs.copy()
        db_seqs = {}
        for uni_ac in input_prots:
            db_seqs[uni_ac] = protein_data.find_one({"uni_ac": uni_ac},
                                                   {"_id": 0, "sequence": 1})["sequence"]
        all_seqs.update(db_seqs)
        st = datetime.datetime.now()
        print "[{}] Running InterPrets".format(st)
        fasta_iprets = run_interprets.main(all_seqs, output_dir, iprets_file, ide,
                       mode="psiblast", psiblast=psiblast_path, blastdb=blastdb_pdb,
                       print_output=True, temp_dir=temp_dir)

    print "[{}] Your input contains {} proteins.".format(st, total_n_prots)
    # " Plus a maximum of {} interactors for each, the total is {} proteins and {} interactions".format(max_prots, len(all_prots), tot_ints)

    # Query Mutations
    input_muts = parse_mutation_input(query_muts, prot_ids, input_prots,
                                      protein_data)

    ## Run int2graph
    graph_ele, lines = int2graph.main(make_graph, sps, all_prots,
                    input_prots, custom_pairs, input_seqs, input_muts,
                    fasta_data, fasta_link, protein_data, cosmic_data, biogrid_data,
                    iprets_data, fasta_iprets, db3did_data, dd_ass_data, elm_int_data,
                    pfam_info, elm_info)

    st = datetime.datetime.now()
    print "[{}] int2graph.py run successfully".format(st)

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

    st = datetime.datetime.now()
    print "[{}] ...done!. Created files \"{}\" and \"{}\"".format(st,
                                                                  graph_json,
                                                                  table_json)
    ## Print HTML
    sys.stdout = sys.__stdout__

    return render_template("results_page.html",
                           graph_json = "jobs/"+"job_"+ide+"/"+graph_json,
                           ints_json = "jobs/"+"job_"+ide+"/"+table_json)
