#!/usr/bin/env python3

""" PIV main script

JC Gonzalez Sanchez, 2020
"""

import os, sys, re, itertools, psycopg2
import gzip, json, csv, random, string, datetime, argparse
from collections import defaultdict
import int2graph, find_all_slims, run_interprets
# from flask import render_template, url_for
# from flask_debugtoolbar_lineprofilerpanel.profile import line_profile

def open_file(input_file, mode="rt"):
    """ Open file Zipped or not
    """

    if re.search(".gz$", input_file):
        infile = gzip.open(input_file, mode)
    else:
        infile = open(input_file, mode)
    return infile

def print_log(ide, msg):
    st = "[{}]".format(datetime.datetime.now())+" [JOB ID: "+ide+"]"
    print(st, msg)

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

def get_pfam_info(cursor, table_name):
    data = {}
    cursor.execute("SELECT accession, type, identifier, description "+
                   "FROM "+table_name+";")
    for row in cursor.fetchall():
        data[row[0]] = {
            "type": row[1],
            "ide":  row[2],
            "des":  row[3]
        }
    return data

def get_elm_info(cursor, table_name):
    data = {}
    cursor.execute("SELECT Accession, Identifier, Name, Description, Regex, Probability "+
                   "FROM "+table_name+";")
    for row in cursor.fetchall():
        data[row[0]] = {
            "ide":   row[1],
            "name":  row[2],
            "des":   row[3],
            "regex": row[4],
            "prob":  row[5]
            }
    return data

def get_3did_dmi(cursor, table_name):
    info = {}
    dmi = defaultdict(dict)
    cursor.execute("SELECT motif, domain_acc, regex, pdb_number "+
                    "FROM "+table_name+";")
    for row in cursor.fetchall():
        info["3DID:"+row[0].upper()] = {
            "ide": row[0],
            "regex": row[2],
            "prob": "-"
        }
        dmi[(row[0].upper(), row[1])] = row[3]
    return info, dmi

def parse_protein_input(input_text, main_org, cursor, id_map):
    input_prots = set()
    input_seqs = defaultdict(str)
    input2uni = {}
    not_found = set()
    custom_pairs = []
    org_map = {}
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
                if org:
                    org_map[header] = org
            
            elif (flag == 1 and len(line.split()) == 1
            and not re.search("\d",line)):
                input_seqs[header] += line.rstrip()
            
            else:
                
                ## 2. Protein identifier(s)
                flag = 0
                line = line.replace(","," ")
                values = line.rstrip().upper().split()
                found = []
                for val in values:
                    uni_id = None
                    if val in input2uni:
                        uni_id = input2uni[val]
                    else:
                        uni_id, org = identify_protein(cursor, id_map, val, main_org)
                        if uni_id:
                            input_prots.add(uni_id)
                            input2uni[val] = uni_id
                            org_map[uni_id] = org

                    if uni_id:
                        found.append(uni_id)
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

    return input_prots, input_seqs, custom_pairs, list(not_found), org_map, input2uni

def parse_mutation_input(input_text, input_prots, input2uni, org_map,
                         main_org, cursor, id_map):
    mutations = defaultdict(lambda: defaultdict(set))

    for line in input_text.split("\n"):
        if line.strip() and line[0] != "#":
            try:
                prot, wt, pos, mt = re.search("(.+)\/(\w)(\d+)(\w+)", line.rstrip()).group(1,2,3,4)
                pos = int(pos)
            except:
                continue

            if prot in input2uni:
                uni_id = input2uni[prot]
                mutations[uni_id][pos].add(wt+str(pos)+mt)
            else: 
                uni_id, org = identify_protein(cursor, id_map, prot, main_org)
                if uni_id:
                    input_prots.add(uni_id)
                    input2uni[prot] = uni_id
                    org_map[uni_id] = org
                else:
                    uni_id = prot
                mutations[uni_id][pos].add(wt+str(pos)+mt) # saved in case they belong to FASTA headers

    return mutations, input_prots, input2uni, org_map

def identify_protein(cursor, id_map, protein_id, org):
    if org:
        cursor.execute("SELECT uniprot_id FROM "+id_map+
                       " WHERE (id = \'"+protein_id+"\' AND organism = \'"+org+"\');")
        results = cursor.fetchone()
        if results:
            return results[0], org
        else:
            return None, None
    else:
        cursor.execute("SELECT uniprot_id, organism FROM "+id_map+
                       " WHERE (id = \'"+protein_id+"\');")
        results = cursor.fetchall()
        if results and len(results) == 1:
            return results[0], results[1]
        else:
            return None, None

def get_protein_data_sql(cursor, sql_table, uni_id):
    cursor.execute("SELECT uniprot_acc, gene, description, length, biogrid_id,"+
                    " sorted_ints, sequence"+
                    " FROM "+sql_table+
                    " WHERE (uniprot_id = \'"+uni_id+"\');")
    row = cursor.fetchone()
    bioids = ""
    if row[4]:
        bioids = row[4]
    ints = []
    if row[5]:
        ints = row[5].split(", ")
    data = {
        "uni_ac":       row[0],
        "genes":        row[1].split("; "),
        "description":  row[2],
        "length":       int(row[3]),
        "biogrid_ids":  bioids,
        "sorted_ints":  ints,
        "seq":          row[6]
    }
    return data

def get_pfam_matches_sql(cursor, sql_table, uni_ac):
    matches = []
    cursor.execute("SELECT env_start, env_end, hmm_ac, e_value"+
                    " FROM "+sql_table+
                    " WHERE (uniprot_acc = \'"+uni_ac+"\');")
    for row in cursor.fetchall():
        matches.append({
            "start": int(row[0]),
            "end": int(row[1]),
            "acc": row[2],
            "e-value": float(row[3])
        })
    return matches

def get_lm_hits_sql(cursor, sql_table, uni_ac):
    data = defaultdict(list)
    cursor.execute("SELECT motif, source, motif_start, motif_end, sequence, status"+
                " FROM "+sql_table+
                " WHERE (uniprot_acc = \'"+uni_ac+"\');")
    for row in cursor.fetchall():
        data[row[0]].append({
            "source": row[1],
            "start":  row[2],
            "end":    row[3],
            "seq":    row[4],
            "status": row[5]
        })
    return data

def get_ptms_sql(cursor, sql_table, uni_ac):
    data = defaultdict(list)
    cursor.execute("SELECT position, residue, ptm"+
                   " FROM "+sql_table+
                   " WHERE (uniprot_acc = \'"+uni_ac+"\');")
    for row in cursor.fetchall():
        data[row[2]].append({
            "pos": int(row[0]),
            "res": row[1]
        })

    return data

def get_uni_feats_sql(cursor, sql_table, uni_ac):
    data = defaultdict(lambda: defaultdict(list))
    for role in ["VARIANT"]:
        cursor.execute("SELECT type, start_pos, end_pos, note"+
                    " FROM "+sql_table+
                    " WHERE (uniprot_acc=\'"+uni_ac+"\' AND type=\'"+role+"\');");
        for row in cursor.fetchall():
            try:
                var, info = re.search("(.+)\((.+)\)", row[3]).group(1, 2)
            except:
                var, info = re.search("(.+)", row[3]).group(1), ""

            if len([x for x in ["cancer","melanoma","carcinoma"] if x in info])==0:
                dis = re.findall("in ([A-Z]+)", info)
                if len(dis)>0:
                    var = var.replace(" ","")
                    start = int(row[1])
                    end = int(row[2])
                    data[row[0]][(start,end)].append({
                        "var": var,
                        "info": info
                    })
    return data

def get_cosmic_mutations_sql(cursor, sql_table, uni_ac, min_sample_size=1):
    data = defaultdict(list)
    cursor.execute("SELECT cosmic_id, cds_mut, aa_mut, sample_num, cancer_types"+
                   " FROM "+sql_table+
                   " WHERE (uniprot_acc = \'"+uni_ac+"\');");
    for row in cursor.fetchall():
        pos = re.search("(\d+)", row[2]).group(1)
        if row[3] >= min_sample_size:
            data[pos].append({
                "cosmic_id": row[0],
                "cds_mut": row[1],
                "aa_mut": row[2],
                "sample_num": row[3],
                "cancer_types": row[4]
            })
    return data  

def get_DDI_sql(cursor, sql_table):
    data = {}
    cursor.execute("SELECT Pfam_Acc_A, Pfam_Acc_B, Source, PDBs "+
                    "FROM "+sql_table+";") 
    for row in cursor.fetchall():
        pfam_a, pfam_b, source = row[0], row[1], row[2]
        pdbs = []
        if row[3]:
            pdbs = str(row[3]).split(";")
        data[(pfam_a, pfam_b)] = {"dbs":source, "pdbs": pdbs}
        data[(pfam_b, pfam_a)] = {"dbs":source, "pdbs": pdbs}
    return data

def get_ELM_dom_sql(cursor, sql_table):
    data = defaultdict(list)
    cursor.execute(
        "SELECT elm_ac, elm_id, domain_ids, present_in_taxon, not_present_in_taxon, elm_containing_genes, "+
        "dom_containing_genes_hsa, dom_containing_genes, phosphosites, observations, other_elm_required "+
        "FROM "+sql_table+";")

    for row in cursor.fetchall():
        elm = row[0].upper()
        doms = string2list_fix(row[2].upper())
        for dom in doms:
            data[(elm, dom)].append( {
                "in_taxon"        : string2list_fix(row[3]),
                "not_in_taxon"    : string2list_fix(row[4]),
                "elm_genes"       : string2list_fix(row[5]),
                "human_dom_genes" : string2list_fix(row[6]),
                "dom_genes"       : string2list_fix(row[7]),
                "phos"            : string2list_fix(row[8]),
                "obs"             : format_none(row[9]),
                "req_elms"        : string2list_fix(row[10])
            })
    return data

def get_interprets_sql(cursor, sql_table, id_a, id_b):
    data = []

    cursor.execute("SELECT pdb_a, blast_eval_a, blast_pcid_a, "+
                    "uni_start_a, uni_end_a, pdb_start_a, pdb_end_a, "+
                    "pdb_b, blast_eval_b, blast_pcid_b, "+
                    "uni_start_b, uni_end_b, pdb_start_b, pdb_end_b, "+
                    "z_score, p_value"+
                    " FROM "+sql_table+
                    " WHERE (uniprot_id_a=\'"+id_a+"\' AND uniprot_id_b=\'"+id_b+"\');");
    for row in cursor.fetchall():
        if row[0]==None:
            data.append(None)
            break
        else:
            a, b = row[:7], row[7:-2] 
            z_value, p_value = float(row[-2]), float(row[-1])
            data.append({
                "a": a, "b": b, "z": z_value, "p": p_value
            })

    return data

def format_none(string):
    if string==None:
        return ""
    else:
        return string

def string2list_fix(string):
    return [str(x) for x in str(string).upper().replace(" ","").split(",") if x!="" and x!="NONE"]

def review_interprets_hits(new_hits, protein_data):
    mask, mask2 = {}, {}
    final_hits = defaultdict(list)
    for pair in new_hits:
        (id_a, id_b) = pair
        len_a, len_b = protein_data[id_a]["length"], protein_data[id_b]["length"]
        for z in sorted(new_hits[pair], reverse=True):
            for hit in new_hits[pair][z]:
                if pair not in final_hits:
                    final_hits[pair].append(hit)
                else:
                    start_a, end_a = int(hit["info_a"][3]), int(hit["info_a"][4])
                    start_b, end_b = int(hit["info_b"][3]), int(hit["info_b"][4])
        
                    mask[id_a] = ["0"]*len_a
                    mask[id_b] = ["0"]*len_b
                    mask2[id_a] = fill_mask(start_a, end_a, ["0"]*len_a)
                    mask2[id_b] = fill_mask(start_b, end_b, ["0"]*len_b)
                    
                    for fhit in final_hits[pair]:
                        fstart_a, fend_a = int(fhit["info_a"][3]), int(fhit["info_a"][4])
                        fstart_b, fend_b = int(fhit["info_b"][3]), int(fhit["info_b"][4])
                        mask[id_a] = fill_mask(fstart_a, fend_a, mask[id_a])
                        mask[id_b] = fill_mask(fstart_b, fend_b,mask[id_b])
                        ov2_a = calculate_overlap(fstart_a, fend_a, mask2[id_a])
                        ov2_b = calculate_overlap(fstart_b, fend_b, mask2[id_b])
                        if ov2_a>=0.75 and ov2_b>=0.75:
                            break
                    else:      
                        ov_a = calculate_overlap(start_a, end_a, mask[id_a])
                        ov_b = calculate_overlap(start_b, end_b, mask[id_b])
                        if ov_a<0.75 and ov_b<0.75:
                            final_hits[pair].append(hit)
    return final_hits

def fill_mask(start, end, mask):
    for i in range(start-1, end):
        mask[i] = "1"
    return mask

def calculate_overlap(start, end, mask):
    sub_mask = mask[start-1:end]
    overlap = sub_mask.count("1") / float(len(sub_mask))
    return overlap
###

def hasNumbers(inputString):
    return any(char.isdigit() for char in inputString)

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

def get_additional_interactors(protein_data, n_ints):
    interactors = defaultdict(list)

    for uni_id, data in protein_data.items():
        if n_ints == "all":
            interactors[uni_id] = data["sorted_ints"]
        else:
            n = 0
            for uni_id_2 in data["sorted_ints"]:
                if uni_id_2 not in protein_data:
                    interactors[uni_id].append(uni_id_2)
                    n += 1
                    if n == int(n_ints):
                        break

    return interactors
    
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

def get_stats(columns, lines, p):
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
        (a, b) = sorted([gene_a, gene_b])
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
    p["table_columns"] = columns
    p["prot_ints_k"], p["prot_ints_v"] = prot_ints_k, prot_ints_v
    p["dom_ints_k"], p["dom_ints"] = dom_ints_k, dom_ints
    p["int_types_k"], p["int_types_series"] = int_types_k, int_types_series
    p["max_ppi"] = max_ppi
    p["dom_ints_per_type"] = dom_ints_per_type

    return p

def print_sorted_dict(d):
    for k in sorted(d, key=d.get, reverse=True):
        print(k+"\t"+str(d[k]))

def percentage(n, tot):
    return "{:3.1f}".format(float(n)/tot*100)

def print_stats_summary(stats_file):
    with open(stats_file, "r") as f:
        d = json.load(f)

    all_prots = d["prot_ints_k"]
    max_ints_per_protein = len(all_prots) -1
    max_ppi = int(d["max_ppi"])

    print("## STATS")
    print("## Proteins ranked by number of interactors. Max:", max_ints_per_protein)
    ints_per_protein = {}
    for k, v in reversed(zip(d["prot_ints_k"], d["prot_ints_v"])):
        print("\t".join([k, str(v), percentage(v, max_ints_per_protein)]))
        ints_per_protein[k] = v

    print("\n## Interaction types ranked by number of linked protein pairs. Max:", max_ppi)
    for k, v in reversed(zip(d["int_types_k"], d["int_types_series"])):
        v = v["value"]
        print("\t".join([k, str(v), percentage(v, max_ppi)]))

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

# @line_profile
def main(INPUT_1=None, INPUT_2=None, ORG="HUMAN", 
         MAIN_OUTPUT_DIR="", CUSTOM_ID=False,
         DATA_DIR="static/data/",
         BLASTDB_DIR="/net/home.isilon/ds-russell/blastdb/",
         PSQL_USER="bq_jgonzalez", PSQL_DB="piv", 
         MAKE_NETWORK=True, HIDE_NO_INT=True, TABLE_FORMAT="json",
         ADDITIONAL_INTERACTORS=0, ONLY_INT_PAIRS=True,
         CMD_LINE=False, RUN_IPRETS=True):

    try:
        add_n_ints = int(ADDITIONAL_INTERACTORS)
    except:
        add_n_ints = 0
        if str(ADDITIONAL_INTERACTORS).lower()=="all":
            add_n_ints = "all"

    # organisms = {
    # "ARATH": "ARATH", "Ath": "ARATH", "3702": "ARATH",  
    # "Cel": ["CAEEL", "6239"],
    # "Dme": ["DROME", "7227"],
    # "Dre": ["DANRE", "7955"],
    # "Hsa": ["HUMAN", "9606"],
    # "Mmu": ["MOUSE", "10090"],
    # "Sce": ["YEAST", "559292"],
    # "Xla": ["XENLA", "8355"],
    # "Xtr": ["XENTR", "8364"]
    # }
    # if SP!="any":

    error = False
    param = {} # Parameters to print in JSON file

    ## DATA: Set PSQL tables
    ID_MAP_TABLE       = "id_mapping"
    PROTEIN_DATA_TABLE = "protein_data"
    UNI_FEAT_TABLE     = "uniprot_features"
    PFAM_DATA_TABLE    = "pfam_a_data"
    PFAM_MATCHES_TABLE = "pfam_a_matches"  
    LM_HITS_TABLE      = "lm_hits"             
    PTM_TABLE          = "ptms"
    ELM_CLASSES_TABLE  = "elm_classes"
    COSMIC_TABLE       = "cosmic_genome_screens"
    PPI_DATA           = "ppi_db"
    DDI_TABLE          = "ddi_db"
    ELM_DOM_TABLE      = "elm_int_dom"
    DMI_3DID_TABLE     = "dmi_3did"
    IPRETS_TABLE       = "interprets_pdb2019"
    conn = psycopg2.connect(database=PSQL_DB, 
                            user=PSQL_USER)
    cursor = conn.cursor()

    BLASTPGP_PATH = "blastpgp"
    PSIBLAST_PATH = "psiblast"
    BLASTDB_PDB   = BLASTDB_DIR+"pdbaa_2019"#"pdbseq"
    # BLASTDB_UNI     = BLASTDB_DIR+"uniprot_sprot_PIV"
    # if SP != "any":
        # if os.path.isfile(BLASTDB_DIR+"uniprot_"+SP+".pin"):
        #     BLASTDB_UNI = BLASTDB_DIR+"uniprot_"+SP

    ### Create Job directory
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
    
    ### Set output files names
    # fasta_file      = OUTPUT_DIR+"seqs_"+IDE+".fasta"
    # blastout_file   = OUTPUT_DIR+"blastout_"+IDE+".tsv.gz"
    # pfamout_file    = OUTPUT_DIR+"pfamscan_"+IDE
    # lmsout_file     = OUTPUT_DIR+"elm_hits_"+IDE+".tsv.gz"
    iprets_file     = OUTPUT_DIR+"i2_"+IDE+".tsv.gz"
    graph_json      = "graph_elements_"+IDE+".json"
    table_file      = "interaction_table_"+IDE+".json"
    if TABLE_FORMAT == "tsv":
        table_file = "interaction_table_"+IDE+".tsv.gz"
    stats_file      = "req_parameters_"+IDE+".json"

    print_log(IDE, "Running PIV")

    ### Parse protein input
    (input_prots, input_seqs,
        custom_pairs, not_found,
        org_map, input2uni) = parse_protein_input(INPUT_1, ORG, 
                                                  cursor, ID_MAP_TABLE)
    
    ### Parse Query Mutations
    (input_muts, input_prots, input2uni, 
        org_map) = parse_mutation_input(INPUT_2, input_prots, input2uni, org_map,
                                        ORG, cursor, ID_MAP_TABLE)
    
    all_proteins = set(input_prots)
    if custom_pairs:
        protein_pairs = custom_pairs
    else:
        protein_pairs = list(itertools.combinations(input_prots, 2))
    
    if input_prots:
        print_log(IDE, "The input contains {} proteins and {} pairs".format(
                        len(input_prots), len(protein_pairs)))
        if not_found:
            print_log(IDE,
                 ("The following {} input protein(s) could not be identified: ".format(len(not_found))+
                 "; ".join(not_found)))
    else:
        print_log(IDE, "ERROR: no valid proteins found in input!")
        return True
    
    protein_data, pfam_matches, lms, ptms = {}, {}, {}, {}
    uni_feats, cosmic_muts = {}, {}
    for uni_id in input_prots:
        protein_data[uni_id] = get_protein_data_sql(cursor, PROTEIN_DATA_TABLE, uni_id)
    
    ### Get additional interaction partners
    if (add_n_ints == "all" or add_n_ints > 0):
        ad_interactors = get_additional_interactors(protein_data, add_n_ints)
        extra_prots = set()
        for uni_id in ad_interactors:
            for uni_id_2 in ad_interactors[uni_id]:
                extra_prots.add(uni_id_2)
                org_map[uni_id_2] = org_map[uni_id]
                if uni_id_2 not in protein_data:
                    protein_data[uni_id_2] = get_protein_data_sql(cursor, PROTEIN_DATA_TABLE, uni_id_2)
                if ONLY_INT_PAIRS:
                    protein_pairs.append( sorted([uni_id, uni_id_2]) )
        
        all_proteins = set(input_prots) | extra_prots
        if not ONLY_INT_PAIRS:
            protein_pairs = list(itertools.combinations(all_proteins, 2))
        
        print_log(IDE, ("Added {} additional interactors per protein ".format(add_n_ints)+
                "for a total of {} proteins".format(
                len(all_proteins)+len(input_seqs))))


    ## Load data
    pfam_info = get_pfam_info(cursor, PFAM_DATA_TABLE)
    elm_info = get_elm_info(cursor, ELM_CLASSES_TABLE)
    ddi = get_DDI_sql(cursor, DDI_TABLE)
    dmi_elm = get_ELM_dom_sql(cursor, ELM_DOM_TABLE)
    lm3did_info, dmi_3did = get_3did_dmi(cursor, DMI_3DID_TABLE)

    for uni_id in all_proteins:
        uni_ac = protein_data[uni_id]["uni_ac"]
        pfam_matches[uni_id] = get_pfam_matches_sql(cursor, PFAM_MATCHES_TABLE, uni_ac)        
        lms[uni_id] = get_lm_hits_sql(cursor, LM_HITS_TABLE, uni_ac)
        ptms[uni_id] = get_ptms_sql(cursor, PTM_TABLE, uni_ac)
        uni_feats[uni_id] = get_uni_feats_sql(cursor, UNI_FEAT_TABLE, uni_ac)
        if org_map[uni_id] == "HUMAN":
                    cosmic_muts[uni_id] = get_cosmic_mutations_sql(cursor, COSMIC_TABLE, uni_ac)


    #### PUT THIS IN A DIFFERENT SCRIPT
    ### 4. If sequence(s) in input -> compute data:
    fasta_data = defaultdict(lambda: defaultdict(list))
    fasta_link = {}
    fasta_iprets = {}
    # input_seqs = {} ##################### DEACTIVATED
    # if len(input_seqs) > 0:
    #     mask = defaultdict(list)

    #     ### 4-1. Print sequences in file
    #     with open_file(fasta_file, "w") as out:
    #         for head in input_seqs:
    #             fasta_link[head] = "static/jobs/job_"+IDE+"/seqs/"+head+".fa"
    #             mask[head] = ["0"] * len(input_seqs[head])
    #             out.write(">"+head+"\n")
    #             out.write(input_seqs[head]+"\n")

    #     ### 4-2. Run BLAST to identify similar proteins to input sequence(s)
    #     ## Blast parameters:
    #     psiblast_path = "psiblast"
    #     e_val = "0.01" #Default = 10
    #     ite = "1"  #-num_iterations (Default = 1)
    #     outfmt = "7"
    #     ## Run blast through bash.
    #     cmd = "{} -query {} -db {} -evalue {} -outfmt {} | gzip > {}".format(
    #           psiblast_path, fasta_file, BLASTDB_UNI, e_val, outfmt,
    #           blastout_file)
    #     print_log(IDE, "Running PSIBLAST: {}".format(cmd))
    #     os.system(cmd)
    #     ## Parse results
    #     fasta_data = parse_blast(blastout_file, mask, fasta_data)

    #     ### 4-3. Run PfamScan to identify Pfam domains on input sequence(s).
    #     # pfamscan.py parameters:
    #     pfamscan_path = "/var/www/flask_apps/jc_test/jc_app/pfamscan.py"
    #     evalue = "0.001"
    #     email = "juan-carlos.gonzalez@bioquant.uni-heidelberg.de"
    #     cmd = ("python {} --sequence {} --database pfam-a --evalue {}".format(
    #           pfamscan_path, fasta_file, evalue)+
    #           " --format txt --outfile {} --email {} --quiet".format(
    #           pfamout_file, email))
    #     print_log(IDE,"Running PfamScan: {}".format(cmd))
    #     os.system(cmd)
    #     os.system("gzip "+pfamout_file+".out.txt")
    #     os.unlink(pfamout_file+".sequence.txt")
    #     # os.system("gzip "+pfamout_file+".sequence.txt")
    #     ## Parse results
    #     fasta_data = parse_pfamscan(pfamout_file+".out.txt.gz", fasta_data)

    #     ### 4-4. Run ELM search on input sequence(s).
    #     fpm2_path = "/var/www/flask_apps/jc_test/jc_app/fpm2.pl"
    #     print_log(IDE, "Searching for ELMs")
    #     all_lms = elm_info.copy()
    #     all_lms.update(lm3did_info)
    #     find_all_slims.main(fasta_file, all_lms,
    #                         print_out=True, outfile=lmsout_file, tmp_file="/tmp/lm_search_"+IDE+".txt")
    #     fasta_data = extract_linear_motifs(lmsout_file, fasta_data)

    #     # 4-5. Run InterPreTS
    #     if RUN_IPRETS:
    #         all_seqs = input_seqs.copy()
            # if len(custom_pairs)==0:
            #     # between input sequences with each other
            #     custom_pairs = list(itertools.combinations(input_seqs.keys(), 2))
            #     # between input sequence with other proteins
            #     custom_pairs += list(itertools.product(input_seqs.keys(), list(all_proteins)))

            #     for uni_ac in all_proteins:
            #         all_seqs[uni_ac] = PROTEIN_DATA.find_one({"uni_ac": uni_ac},
            #                                     {"_id": 0, "sequence": 1})["sequence"]
            #     print_log(IDE, "Running InterPrets")
            #     fasta_iprets = run_interprets.main(
            #                         all_seqs, OUTPUT_DIR, iprets_file, IDE,
            #                         these_pairs=custom_pairs,
            #                         mode="psiblast", psiblast=psiblast_path, blastdb=BLASTDB_PDB,
            #                         print_output=True)
            # else:
            #     for pair in custom_pairs:
            #         for prot in pair:
            #             if prot in all_proteins:
            #                 all_seqs[uni_ac] = PROTEIN_DATA.find_one({"uni_ac": uni_ac},
            #                                     {"_id": 0, "sequence": 1})["sequence"]
            #     print_log(IDE, "Running InterPrets")
            #     fasta_iprets = run_interprets.main(
            #             all_seqs, OUTPUT_DIR, iprets_file, IDE,
            #             these_pairs=custom_pairs,
            #             mode="psiblast", psiblast=psiblast_path, blastdb=BLASTDB_PDB,
            #             print_output=True)

    iprets_hits = defaultdict(list)
    run_seqs, run_pairs = {}, []
    for pair in protein_pairs:
        (uni_id_a, uni_id_b) = sorted(pair)
        hits = get_interprets_sql(cursor, IPRETS_TABLE, uni_id_a, uni_id_b)
        if hits:
            iprets_hits[(uni_id_a, uni_id_b)] = hits
        else:
            run_pairs.append((uni_id_a, uni_id_b))
            for uni_id in (uni_id_a, uni_id_b):
                run_seqs[uni_id] = protein_data[uni_id]["seq"]
    
    if run_pairs:
        print_log(IDE, "Running InterPrets for {} pairs".format(len(run_pairs)))
        new_hits = run_interprets.main( run_seqs, OUTPUT_DIR, iprets_file, IDE, org_map,
                        these_pairs=run_pairs,
                        mode="blastpgp", blastpgp=BLASTPGP_PATH, psiblast=PSIBLAST_PATH, 
                        blastdb=BLASTDB_PDB, data_dir=DATA_DIR, print_output=False)
        
        final_hits = review_interprets_hits(new_hits, protein_data)

        for pair in run_pairs:
            (id_a, id_b) = pair
            ac_a, ac_b = protein_data[id_a]["uni_ac"], protein_data[id_b]["uni_ac"]
            if pair in final_hits:
                for hit in final_hits[pair]:
                    vals = ["\'"+x+"\'" for x in [id_a, ac_a]+hit["info_a"]+[id_b, ac_b]+hit["info_b"]+hit["scores"] ]
                    cursor.execute("INSERT INTO "+IPRETS_TABLE+" VALUES ("+", ".join(vals)+");")
                    conn.commit()

                    iprets_hits[pair].append({
                        "a": hit["info_a"],
                        "b": hit["info_b"],
                        "z": float(hit["scores"][-4]),
                        "p": float(hit["scores"][-3])
                    })
            else:
                vals = ["\'"+x+"\'" for x in [id_a, ac_a, id_b, ac_b] ]
                cursor.execute("INSERT INTO "+IPRETS_TABLE+
                               " (uniprot_id_a, uniprot_ac_a, uniprot_id_b, uniprot_ac_b)"+
                               " VALUES ("+", ".join(vals)+");")
                conn.commit()
  
    ### 6. Run int2graph
    graph_ele, columns, lines, no_int_prots = int2graph.main(
            all_proteins, protein_pairs, input_seqs, input_muts, org_map,
            conn, protein_data, pfam_matches, lms, ptms, uni_feats, cosmic_muts,
            PPI_DATA, ddi, dmi_elm, dmi_3did, pfam_info, elm_info, iprets_hits,
            fasta_data, fasta_link, fasta_iprets,
            make_network=MAKE_NETWORK, hide_no_int=HIDE_NO_INT)
    
    conn.close()

    if len(no_int_prots) > 0:
        print_log(IDE, "No interaction evidence found for {} proteins.".format(
                  len(no_int_prots)))

    print_log(IDE, "int2graph.py run successfully")
    param["not_found"], param["no_int_prots"] = not_found, no_int_prots

    ### 7. Print Output files
    ## Print graph as JSON file
    if MAKE_NETWORK:
        with open(OUTPUT_DIR+graph_json, "wt") as output:
            json.dump(graph_ele, output, indent=4)
        print_log(IDE, "Created \"{}\"".format(graph_json))

    ## Print interactions
    ## as TSV
    if TABLE_FORMAT=="tsv":
        with open_file(OUTPUT_DIR+table_file, "wt") as output:
            tsvwriter = csv.writer(output, delimiter="\t")
            tsvwriter.writerow(columns)
            for data in lines:
                tsvwriter.writerow(data)
    ## as JSON
    else:
        int_table = {}
        int_table["columns"] = columns
        int_table["index"] = list(range(len(lines)))
        int_table["data"] = lines
        with open(OUTPUT_DIR+table_file, "wt") as output:
            json.dump(int_table, output)
    print_log(IDE, "Created \"{}\"".format(table_file))

    ### 8. Compute Stats from Interaction Table
    param = get_stats(columns, lines, param)
    with open(OUTPUT_DIR+stats_file, "wt") as output:
        json.dump(param, output)
    print_log(IDE, "Created \"{}\"".format(stats_file))

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
    ap.add_argument("-sp", "--species", required=False, default="HUMAN",
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

    main(INPUT_1=prots, INPUT_2=muts, ORG=args["species"],
         ADDITIONAL_INTERACTORS=args["add_ints"],
         CUSTOM_ID=args["job_id"],
         MAKE_NETWORK=True, TABLE_FORMAT="json", CMD_LINE=True,
         RUN_IPRETS=False)

    sys.exit()
