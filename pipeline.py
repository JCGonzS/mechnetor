#!/usr/bin/env python

""" PIV main script

JC Gonzalez Sanchez, 2018
"""

import os, sys, re, itertools
import gzip, json, csv, random, string, datetime, pymongo, argparse
import pandas as pd
import run_interprets, int2graph, find_all_slims
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

def make_protein_id_dictionary(prot_data, org):
    D = defaultdict(lambda: defaultdict(list))
    for doc in prot_data.find({"organism": org},
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

# def function(input_text, prot_dict):
#     for line in input_text.split("\n"):

#         ## Make proteins UniProtID centric
#         if line.strip() and line[0] != "#":

def parse_protein_input(input_text, sp, protein_data, prot_dict):
    input_prots, done = set(), set()
    input_seqs = defaultdict(str)
    input_to_uniac = {}
    not_found = set()
    custom_pairs = []
    sp_map = {}
    flag = 0

    for line in input_text.split("\n"):

        ## Make proteins UniProtID centric
        if line.strip() and line[0] != "#":
            ## Types of input:
            ## 1. FASTA sequence
            if line[0] == ">":
                flag = 1
                header = str(line.rstrip().split()[0].split(">")[1])
                header = header.replace("/","-")
                if sp!="any":
                    sp_map[header] = sp
            elif (flag == 1 and len(line.split()) == 1
            and not re.search("\d",line)):
                input_seqs[header] += line.rstrip()
            else:
                ## 2. Protein identifier(s)
                flag = 0
                line = line.replace(","," ")
                vals = line.rstrip().upper().split()
                found = []
                for val in vals:
                    ac, org = None, None
                    if val in input_to_uniac:
                        ac = input_to_uniac[val]
                    elif sp != "any":
                        if val in prot_dict[sp]["AC"]:
                            ac = prot_dict[sp]["AC"][val][0]
                            input_to_uniac[val] = ac
                            sp_map[ac] = sp
                            input_prots.add(ac)
                    else:
                        doc = protein_data.find_one(
                                {"$or": [{"uni_id": val},{"uni_ac": val}]},
                                {"_id": 0, "uni_ac": 1, "organism": 1})
                        if doc:
                            ac = doc["uni_ac"]
                            input_to_uniac[val] = ac
                            sp_map[ac] = doc["organism"]
                            input_prots.add(ac)
                    if ac:
                        found.append(ac)
                    else:
                        not_found.add("'"+val+"'")
                ## If pair of proteins in line, save as custom pair.
                if len(found) == 2:
                    custom_pairs.append(sorted(found))
        else:
            flag = 0

    ## Remove FASTA headers with empty sequences
    remove_keys = []
    if len(input_seqs) > 0:
        for key in input_seqs:
            if len(input_seqs[key]) == 0:
                remove_keys.append(key)
    for key in remove_keys:
        input_seqs.pop(key, "None")

    return input_prots, input_seqs, custom_pairs, list(not_found), sp_map, input_to_uniac

def get_additional_interactors(proteins, sp_map, add_ints, protein_data):
    all_proteins = proteins

    add_all_ints = False
    if add_ints in ["all", "All", "ALL"]:
        add_all_ints = True
    else:
        try:
            add_ints = int(add_ints)
        except:
            add_ints = 0

    if (add_ints > 0 or add_all_ints == True):
        for uni_ac in proteins:
            doc = protein_data.find_one({"uni_ac": uni_ac},
                                        {"_id":0, "biogrid_interactors":1})

            if add_all_ints==True:
                interactors = doc["biogrid_interactors"]
            else:
                interactors = []
                for uni_ac2 in doc["biogrid_interactors"]:
                    interactors.append(uni_ac2)
                    if (len(interactors) >= int(add_ints)):
                        break
            for prot in interactors:
                sp_map[prot] = sp_map[uni_ac]
            all_proteins = all_proteins | set(interactors)

    return all_proteins, sp_map, add_ints

def parse_mutation_input(input_text, input_to_uniac, input_seqs):
    mutations = defaultdict(lambda: defaultdict(set))

    for line in input_text.split("\n"):
        if line.strip() and line[0] != "#":
            prot, mut = line.rstrip().split()[0].split("/")
            pos = int(re.search("([\d]+)", mut).group(1))
            if prot in input_to_uniac:
                uni_ac = input_to_uniac[prot]
                mutations[uni_ac][pos].add(mut)
            elif prot in input_seqs:
                mutations[prot][pos].add(mut)
            ## add protein if not present in input
            # else:
            #   (...)

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
        d[c["Accession"]] = {
            "ide": c["Identifier"],
            "des": c["Description"],
            "type": c["Type"]}
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
            "prob":  c["Probability"]
            }
    return d

def get_dmi_info(data):
    info = {}
    dmi = defaultdict(dict)
    for c in data.find():
        info["3DID:"+c["MOTIF"]] = {
            "ide": c["MOTIF"],
            "regex": c["REGEX"],
            "prob": "-"
        }
        dmi[c["MOTIF"]][c["DOMAIN"]] = c["PDBS"]
    return info, dmi

def parse_blast(blastout_file, mask, data_dict,
                max_ol=0.25, max_cov=0.9):
    prots = set()
    sp, tr = defaultdict(list), defaultdict(list)
    with open_file(blastout_file) as f:
        for line in f:
            if line[0]!="#":
                t = line.rstrip().split("\t")
                query = t[0]
                dc, acc, uni_id = t[1].split("|")
                prots.add(query)
                if dc=="sp":
                    sp[query].append({ "acc": uni_id,
                                  "ide": float(t[2]),
                                  "e-val": float(t[10]),
                                  "q-start": int(t[6]),
                                  "q-end": int(t[7])
                                  })
                elif dc=="tr":
                    tr[query].append({ "acc": uni_id,
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
            l = str(hit["q-start"])+"-"+str(hit["q-end"])+"|"
            l += hit["acc"]+"|"+"("+"{:1.0e}".format(hit["e-val"])+", "+str(hit["ide"])+"%)"
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

# def find_all_slims(fasta_file, data_dict, ide, elm_info, 
#                     fpm2_script = "fpm2.pl"):

#     tmp_file = "lm_search_"+ide+".txt"

#     for elm_acc in elm_info:
#         elm_ide = elm_info[elm_acc]["ide"]
#         elm_regex = elm_info[elm_acc]["regex"]
#         elm_prob = elm_info[elm_acc]["prob"]

#         mode = "cat"
#         if ".gz" in fasta_file:
#             mode = "zcat"
#         com = "{} {} | perl {} \"{}\" -m 1 > {}".format(mode, fasta_file, fpm2_script, elm_regex, tmp_file)
#         os.system( com )

#         with open_file(tmp_file) as f:
#             for line2 in f:
#                 if line2[0]==">":
#                     label = line2.split()[0].split("/")[0].replace(">","")
#                     start, end = line2.split()[0].split("/")[1].split("-")
#                     gene = "-"
#                     if "GN=" in line.split("/")[1]:
#                         gene = re.search("GN=([^\s]+)", line.split("/")[1]).group(1)
#                 else:
#                     seq = line2.rstrip()
#                     if "elms" in data_dict[label]:
#                         data_dict[label]["elms"][elm_acc].append(
#                                                     {"start" : int(start),
#                                                      "end" :   int(end),
#                                                      "seq" :   seq
#                                                      })
#                     else:
#                         data_dict[label]["elms"]=defaultdict(list)
#     # os.unlink(tmp_file)
#     return data_dict

def extract_linear_motifs(lm_hits_file, data_dict, max_overlap=2, max_eval=1):

    with open_file(lm_hits_file) as f:
        for line in f:
            if line[0] != "#" and line.strip():
                t = line.rstrip().split("\t")
                elm_name, elm_acc = t[0], t[1]
                label = t[2]
                start, end = t[4], t[5]
                seq = t[9]

                group = "elms"
                if elm_acc.startswith("3D"):
                    group = "3dlms"
                    elm_acc = elm_acc.split(":")[1]

                if group in data_dict[label]:
                    data_dict[label][group][elm_acc].append(
                                                { "start": int(start),
                                                  "end": int(end),
                                                  "seq": seq
                                                })
                else:
                    data_dict[label][group] = defaultdict(list)
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

def print_log(ide, msg):
    st = "[{}]".format(datetime.datetime.now())+" [JOB ID: "+ide+"]"
    print st, msg

def get_sorted_lists(d):
    keys, vals = [], []
    for k in sorted(d, key=lambda k: len(d[k])):
        v = len(d[k])
        keys.append(k)
        vals.append(v)
    return keys, vals

def dict_from_set_len(d):
    new_d = {}
    for k in d:
        new_d[k] = len(d[k])
    return new_d

def get_stats(lines, p):
    colors = {
        "PROT::PROT": "#5F6A6A",
        "DOM::DOM": "#16A085",
        "iDOM::iDOM": "#D4AC0D",
        "ELM::DOM": "#AF7AC5",
        "LM::DOM" : "#27ae60",
        "InterPreTS": "#E74C3C",
        "UniProt region": "#EC7063"
    }
    names = {
        "PROT::PROT": "Binary",
        "DOM::DOM": "Domain-Domain",
        "iDOM::iDOM": "(in)Domain-Domain",
        "ELM::DOM": "Domain-Motif",
        "LM::DOM": "Domain-Motif 3did",
        "InterPreTS": "3D structure",
        "UniProt region": "UniProt region"
    }
    types = {"DOM::DOM": "DDI",
             "iDOM::iDOM": "iDDI",
             "ELM::DOM": "DMI",
             "LM::DOM": "DMI-3D"
             }
    prot_ints = defaultdict(set)
    int_types = defaultdict(set)
    dom_ints_tot = defaultdict(set)
    dom_ints_per_type = defaultdict(lambda: defaultdict(list))
    for line in lines:
        # "#Gene(A)","Accession(A)","Gene(B)","Accession(B)",
        #     "Type", "F(A)","Start-End(A)", "Mutations(A)",
        #     "F(B)", "Start-End(B)", "Mutations(B)",
        #     "Info", "Source"]
        gene_a, gene_b = str(line[2]), str(line[3])
        int_type = str(line[6])
        ele_a, ele_b = str(line[7]), str(line[10])
        prot_ints[gene_a].add(gene_b)
        prot_ints[gene_b].add(gene_a)
        a, b = gene_a, gene_b
        if b < a:
            a, b = gene_b, gene_a
        int_types[int_type].add((a, b))

        if "DOM" in int_type:
            ele_a = ele_a.split()[0]
            ele_b = ele_b.split()[0]
            int_type = types[int_type]
            if gene_a not in dom_ints_per_type[gene_b+":"+ele_b][int_type]:
                dom_ints_per_type[gene_b+":"+ele_b][int_type].append(gene_a)
            dom_ints_tot[gene_b+":"+ele_b].add(gene_a)
            if int_type in ["DDI", "iDDI"]:
                if gene_b not in dom_ints_per_type[gene_a+":"+ele_a][int_type]:
                    dom_ints_per_type[gene_a+":"+ele_a][int_type].append(gene_b)
                dom_ints_tot[gene_a+":"+ele_a].add(gene_b)

    prot_ints_k, prot_ints_v = get_sorted_lists(prot_ints)
    max_ppi = len([x for x in itertools.combinations(prot_ints_k, 2)])
    int_types_k, int_types_series = [], []
    for k in sorted(int_types, key=lambda k: len(int_types[k])):
        v = len(int_types[k])
        if k == "iELM::DOM":
            continue
        int_types_k.append(str(names[k]))
        int_types_series.append({
                                "value": v,
                                "itemStyle": {"color": colors[k]}
                                })

    int_types_number = dict_from_set_len(int_types)
    dom_ints_k = []
    dom_ints = defaultdict(list)
    for k in sorted(dom_ints_tot, key=lambda k: len(dom_ints_tot[k])):
        dom_ints_k.append(k)
        dom_ints["total"].append(len(dom_ints_tot[k]))
        for int_type in ["DDI", "iDDI", "DMI"]:
            v = 0
            if int_type in dom_ints_per_type[k]:
                v = len(dom_ints_per_type[k][int_type])
            dom_ints[int_type].append(v)

    ## Add parameters to dictionary
    p["prot_ints_k"], p["prot_ints_v"] = prot_ints_k, prot_ints_v
    p["dom_ints_k"], p["dom_ints"] = dom_ints_k, dom_ints
    p["int_types_k"], p["int_types_series"] = int_types_k, int_types_series
    p["max_ppi"] = max_ppi
    p["dom_ints_per_type"] = dom_ints_per_type

    return p

def print_sorted_dict(d):
    for k in sorted(d, key=d.get, reverse=True):
        print k+"\t"+str(d[k])

def percentage(n, tot):
    return "{:3.1f}".format(float(n)/tot*100)

def print_stats_summary(stats_file):
    with open(stats_file, "r") as f:
        d = json.load(f)

    all_prots = d["prot_ints_k"]
    max_ints_per_protein = len(all_prots) -1
    max_ppi = int(d["max_ppi"])

    print "## STATS"
    print "## Proteins ranked by number of interactors. Max:", max_ints_per_protein
    ints_per_protein = {}
    for k, v in reversed(zip(d["prot_ints_k"], d["prot_ints_v"])):
        print "\t".join([k, str(v), percentage(v, max_ints_per_protein)])
        ints_per_protein[k] = v

    print "\n## Interaction types ranked by number of linked protein pairs. Max:", max_ppi
    for k, v in reversed(zip(d["int_types_k"], d["int_types_series"])):
        v = v["value"]
        print "\t".join([k, str(v), percentage(v, max_ppi)])

    # print "\n## Domains ranked"
    # for k, v in enumerate(reversed(zip(d["dom_ints_k"], d["dom_ints_v"])):
    #     prot = k.split(":")[0]
    #     line = [k, str(v), percentage(v, ints_per_protein[prot])]
    #     for t in ["DDI", "iDDI", "DMI"]:
    #
    # #         if t in d["dom_ints"][k]:
    # #             val = len(d["dom_ints"][k][t])
    # #         per = percentage(val, max_ints_prot)
    # #         line.append(str(val))
    # #         line.append(per)
    #     print "\t".join(line)

@line_profile
def main(INPUT_1=None, INPUT_2=None, SP="any", ADDITIONAL_INTERACTORS=0,
         MAIN_OUTPUT_DIR="", CUSTOM_ID=False,
         BLASTDB_DIR="/net/home.isilon/ds-russell/blastdb/",
         CLIENT=MongoClient('localhost', 27017),
         MAKE_NETWORK=True, HIDE_NO_INT=True, TABLE_FORMAT="json",
         CMD_LINE=False, RUN_IPRETS=True):

    error = False
    param = {} # Parameters to print in JSON file

    ## DATA: Set MongoDB databases & collections (CLIENT[database][collection])
    PFAM_DATA       = CLIENT["common"]["pfamA_data"]
    DDI_DATA        = CLIENT["common"]["ddi_db"]
    DMI_3DID        = CLIENT["common"]["dmi_3did"]
    # DB3DID_DATA     = CLIENT["common"]["db_3did"]
    ELM_INT_DATA    = CLIENT["common"]["elm_int_dom"]
    ELM_CLASSES     = CLIENT["common"]["elm_classes"]
    PROTEIN_DATA    = CLIENT["common"]["protein_data"]
    BLASTDB_UNI     = BLASTDB_DIR+"uniprot_sprot_PIV"
    BLASTDB_PDB     = BLASTDB_DIR+"pdbaa_2019"
    COSMIC_DATA     = CLIENT["cosmic_v87"]["genome_screens"]
    if SP != "any":
        # # PROTEIN_DATA    = CLIENT[SP]["protein_data"]
        # ASS_PROB_DATA   = CLIENT[SP]["association_probabilities"]
        # IPRETS_DATA = CLIENT[SP]["interprets"]

        if os.path.isfile(BLASTDB_DIR+"uniprot_"+SP+".pin"):
            BLASTDB_UNI = BLASTDB_DIR+"uniprot_"+SP

    ## Get output file names
    if CUSTOM_ID:
        IDE = CUSTOM_ID
    else:
        IDE = get_unique_random_identifier(MAIN_OUTPUT_DIR)
    OUTPUT_DIR = MAIN_OUTPUT_DIR+"job_"+IDE+"/"
    if not os.path.exists(OUTPUT_DIR):
        try:
            os.mkdir(OUTPUT_DIR)
        except OSError:
            print_log(IDE,
                      "Creation of the directory {} failed".format(OUTPUT_DIR))

    # if not os.path.exists(TEMP_DIR):
    #     st = "[{}]".format(datetime.datetime.now())
    #     print st, "Error. The specified 'temp' directory"+
    #           "'{}' does not exist'".format(TEMP_DIR)

    fasta_file      = OUTPUT_DIR+"seqs_"+IDE+".fasta"
    blastout_file   = OUTPUT_DIR+"blastout_"+IDE+".tsv.gz"
    pfamout_file    = OUTPUT_DIR+"pfamscan_"+IDE
    lmsout_file    = OUTPUT_DIR+"elm_hits_"+IDE+".tsv.gz"
    iprets_file     = OUTPUT_DIR+"i2_"+IDE+".tsv.gz"
    graph_json      = "graph_elements_"+IDE+".json"
    table_file      = "interaction_table_"+IDE+".json"
    if TABLE_FORMAT == "tsv":
        table_file = "interaction_table_"+IDE+".tsv.gz"
    stats_file      = "req_parameters_"+IDE+".json"

    print_log(IDE, "Running PIV")

    ## Load data
    pfam_info = get_pfam_info(PFAM_DATA)
    elm_info = get_elm_info(ELM_CLASSES)
    lm3did_info, dmi_3did = get_dmi_info(DMI_3DID)


    prot_ids = {}
    if SP != "any":
        ### 1. Get protein dictionary & sequences for species
        prot_ids[SP] = make_protein_id_dictionary(PROTEIN_DATA, SP)

    ### 2. Parse protein input
    (input_prots, input_seqs,
    custom_pairs, not_found,
    sp_map, input_to_uniac) = parse_protein_input(INPUT_1, SP,
                                                    PROTEIN_DATA,
                                                    prot_ids)

    all_sp = set(sp_map.values())
    for sp1 in all_sp:
        if sp1 not in prot_ids:
            prot_ids[sp1] = make_protein_id_dictionary(PROTEIN_DATA, sp1)


    total_n_prots = len(input_prots) + len(input_seqs)
    if total_n_prots == 0:
        print_log(IDE, "ERROR: no valid proteins found in input!")
        return True

    if len(not_found) > 0:
        print_log(IDE,
                 ("The following {} input protein(s) could not be identified: ".format(len(not_found))+
                 "; ".join(not_found)))

    print_log(IDE, "Your input contains {} proteins and {} pairs".format(
                        total_n_prots, len(custom_pairs)))

    ### 3. Add additional interaction partners
    all_proteins, sp_map, add = get_additional_interactors(input_prots, sp_map,
                                              ADDITIONAL_INTERACTORS,
                                              PROTEIN_DATA)

    print_log(IDE, ("Added {} additional interactors per protein ".format(add)+
               "for a total of {} proteins".format(
               len(all_proteins)+len(input_seqs))))

    ### 4. If sequence(s) in input -> compute data:
    fasta_data = defaultdict(lambda: defaultdict(list))
    fasta_link = {}
    fasta_iprets = {}

    if len(input_seqs) > 0:
        mask = defaultdict(list)

        ### 4-1. Print sequences in file
        with open_file(fasta_file, "w") as out:
            for head in input_seqs:
                fasta_link[head] = "static/jobs/job_"+IDE+"/seqs/"+head+".fa"
                mask[head] = ["0"] * len(input_seqs[head])
                out.write(">"+head+"\n")
                out.write(input_seqs[head]+"\n")

        ### 4-2. Run BLAST to identify similar proteins to input sequence(s)
        ## Blast parameters:
        psiblast_path = "psiblast"
        e_val = "0.01" #Default = 10
        ite = "1"  #-num_iterations (Default = 1)
        outfmt = "7"
        ## Run blast through bash.
        cmd = "{} -query {} -db {} -evalue {} -outfmt {} | gzip > {}".format(
              psiblast_path, fasta_file, BLASTDB_UNI, e_val, outfmt,
              blastout_file)
        print_log(IDE, "Running PSIBLAST: {}".format(cmd))
        os.system(cmd)
        ## Parse results
        fasta_data = parse_blast(blastout_file, mask, fasta_data)

        ### 4-3. Run PfamScan to identify Pfam domains on input sequence(s).
        # pfamscan.py parameters:
        pfamscan_path = "/var/www/flask_apps/jc_test/jc_app/pfamscan.py"
        evalue = "0.001"
        email = "juan-carlos.gonzalez@bioquant.uni-heidelberg.de"
        cmd = ("python {} --sequence {} --database pfam-a --evalue {}".format(
              pfamscan_path, fasta_file, evalue)+
              " --format txt --outfile {} --email {} --quiet".format(
              pfamout_file, email))
        print_log(IDE,"Running PfamScan: {}".format(cmd))
        os.system(cmd)
        os.system("gzip "+pfamout_file+".out.txt")
        os.unlink(pfamout_file+".sequence.txt")
        # os.system("gzip "+pfamout_file+".sequence.txt")
        ## Parse results
        fasta_data = parse_pfamscan(pfamout_file+".out.txt.gz", fasta_data)

        ### 4-4. Run ELM search on input sequence(s).
        fpm2_path = "/var/www/flask_apps/jc_test/jc_app/fpm2.pl"
        print_log(IDE, "Searching for ELMs")
        all_lms = elm_info.copy()
        all_lms.update(lm3did_info)
        find_all_slims.main(fasta_file, all_lms,
                            print_out=True, outfile=lmsout_file, tmp_file="/tmp/lm_search_"+IDE+".txt")
        fasta_data = extract_linear_motifs(lmsout_file, fasta_data)

        # 4-5. Run InterPreTS
        if RUN_IPRETS:
            all_seqs = input_seqs.copy()
            if len(custom_pairs)==0:
                # between input sequences with each other
                custom_pairs = list(itertools.combinations(input_seqs.keys(), 2))
                # between input sequence with other proteins
                custom_pairs += list(itertools.product(input_seqs.keys(), list(all_proteins)))

                for uni_ac in all_proteins:
                    all_seqs[uni_ac] = PROTEIN_DATA.find_one({"uni_ac": uni_ac},
                                                {"_id": 0, "sequence": 1})["sequence"]
                print_log(IDE, "Running InterPrets")
                fasta_iprets = run_interprets.main(
                                    all_seqs, OUTPUT_DIR, iprets_file, IDE,
                                    these_pairs=custom_pairs,
                                    mode="psiblast", psiblast=psiblast_path, blastdb=BLASTDB_PDB,
                                    print_output=True)
            else:
                for pair in custom_pairs:
                    for prot in pair:
                        if prot in all_proteins:
                            all_seqs[uni_ac] = PROTEIN_DATA.find_one({"uni_ac": uni_ac},
                                                {"_id": 0, "sequence": 1})["sequence"]
                print_log(IDE, "Running InterPrets")
                fasta_iprets = run_interprets.main(
                        all_seqs, OUTPUT_DIR, iprets_file, IDE,
                        these_pairs=custom_pairs,
                        mode="psiblast", psiblast=psiblast_path, blastdb=BLASTDB_PDB,
                        print_output=True)


    ### 5. Parse Query Mutations
    input_muts = parse_mutation_input(INPUT_2, input_to_uniac, input_seqs.keys())

    ### 6. Run int2graph
    graph_ele, lines, no_int_prots = int2graph.main(all_proteins,
            custom_pairs, input_seqs, input_muts, sp_map,
            fasta_data, fasta_link, prot_ids,
            CLIENT, PROTEIN_DATA, COSMIC_DATA,
            fasta_iprets, DDI_DATA, ELM_INT_DATA, DMI_3DID,
            pfam_info, elm_info,
            make_network=MAKE_NETWORK, hide_no_int=HIDE_NO_INT)

    if len(no_int_prots) > 0:
        print_log(IDE,
                "No interaction evidence found for {} proteins:".format(len(no_int_prots)))

    print_log(IDE, "int2graph.py run successfully")
    param["not_found"], param["no_int_prots"] = not_found, no_int_prots

    ### 7. Print Output files
    ## Print graph as JSON file
    if MAKE_NETWORK:
        with open(OUTPUT_DIR+graph_json, "w") as output:
            json.dump(graph_ele, output, indent=4)

        print_log(IDE, "Created \"{}\"".format(graph_json))

    ## Print interactions
    int_table = {}
    int_table["columns"] = ["#Gene(A)","Accession(A)","Gene(B)","Accession(B)",
        "Type", "F(A)","Start-End(A)", "Mutations(A)",
        "F(B)", "Start-End(B)", "Mutations(B)",
        "Info", "Source"]
    int_table["index"] = range(len(lines))
    int_table["data"] = lines
    ## as TSV
    if TABLE_FORMAT=="tsv":
        with open_file(OUTPUT_DIR+table_file, "w") as output:
            tsvwriter = csv.writer(output, delimiter="\t")
            tsvwriter.writerow(int_table["columns"])
            for data in int_table["data"]:
                tsvwriter.writerow(data)
    ## as JSON
    else:
        with open(OUTPUT_DIR+table_file, "w") as output:
            json.dump(int_table, output)
    print_log(IDE, "Created \"{}\"".format(table_file))


    print_log(IDE, "Created \"{}\"".format(stats_file))

    ### 8. Compute Stats from Interaction Table
    param = get_stats(lines, param)
    with open(OUTPUT_DIR+stats_file, "w") as output:
        json.dump(param, output)

    # if CMD_LINE:
    #     print_stats_summary(OUTPUT_DIR+stats_file)

    print_log(IDE, "Job Done!")

    return error


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-p", "--prots", required=True,
    	help="File with protein input as protein identifier or FASTA sequences (Required)")
    ap.add_argument("-m", "--muts", required=False,
    	help="File with mutation input in Mechismo format")
    ap.add_argument("-sp", "--species", required=False, default="any",
        help="Organism (default: Hsa)")
    ap.add_argument("-ai", "--add_ints", required=False, default=0,
        help="Number of additional interactors per protein. Any number or 'all' (default: 0)")
    ap.add_argument("-id", "--job_id", required=False, default=False,
        help="Custom job name (random by default)")
    ap.add_argument("-log", required=False, default=False,
        help="Print log file")
    args = vars(ap.parse_args())

    prots = open_file(args["prots"]).read()
    muts = ""
    if args["muts"]:
        muts = open_file(args["muts"]).read()

    if args["log"]:
        sys.stdout = open(args["log"], 'w')

    main(INPUT_1=prots, INPUT_2=muts, SP=args["species"],
         ADDITIONAL_INTERACTORS=args["add_ints"],
         CUSTOM_ID=args["job_id"],
         MAKE_NETWORK=False, TABLE_FORMAT="tsv", CMD_LINE=True,
         RUN_IPRETS=False)

    sys.exit()
