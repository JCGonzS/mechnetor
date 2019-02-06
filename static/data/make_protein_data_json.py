#!/usr/bin/env python

"""Preparation of data for PIV tool

This scripts reads and parses different databases (flat files)
and generates a single JSON file that is required.

You only need to run this script after updating the database files

Author: Juan Carlos Gonzalez-Sanchez, 2018
"""

### Add COSMIC mutation data

import sys, re
import json, csv, gzip, pprint
from collections import defaultdict
from Bio import SwissProt


## Data Files
data_dir = ""
common_data_dir = data_dir+"common/"
pfam_data_file = common_data_dir + "Pfam-A.hmm_r32.0.dat.gz"
species = "Hsa"
sp_data_dir = data_dir+"species/"+species+"/"
proteome_file =     sp_data_dir + "uniprot_sprot_human_complete_20303prts_March_2018.fasta.gz"
uniprot_text_file = sp_data_dir + "uniprot_homo_sapiens_proteome_73112prts_Aug2018_data.txt.gz"
pfam_matches_file = sp_data_dir + "pfamA_matches_9606.tsv.gz"
elm_hits_file =     sp_data_dir + "elm_parsed_info_noOverlap.tsv.gz"
psp_file =          sp_data_dir + "PSP_ptms_human.tsv.gz"


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

def assign_Number(inputString):
	if ";" in inputString:
		return inputString.split(";")[0]
	else:
		return inputString

def get_protein_data_from_uniprot_text(uniprot_file):
    """From SwissProt"""

    ## this funcion should also print the
    ## "parsed_swissprot" version that pipeline_v2 uses
    ##
    ## !!!!

    D = defaultdict(dict)
    D["doms"] = defaultdict(dict)
    for record in SwissProt.parse(open_file(uniprot_file)):
        dc = record.data_class
        # if dc != "Reviewed":
        #     continue
        uni_id = record.entry_name
        accs = record.accessions
        uni_ac = accs[0].upper()
        gn = record.gene_name
        if record.gene_name != "":
            gn = re.search("Name[s]?=([^;]+)", record.gene_name ).group(1).split()[0]
        des = re.search("Name: Full=([^;]+)", record.description).group(1)
        seq = record.sequence
        length = record.sequence_length

        D["AC"][uni_ac], D["AC"][uni_id], D["AC"][gn] = uni_ac.upper(), uni_ac.upper(), uni_ac.upper()
        D["AC"][uni_ac.upper()], D["AC"][uni_id.upper()], D["AC"][gn.upper()] = uni_ac.upper(), uni_ac.upper(), uni_ac.upper()
        D["ID"][uni_ac], D["ID"][uni_id], D["ID"][gn] = uni_id.upper(), uni_id.upper(), uni_id.upper()
        D["ID"][uni_ac.upper()], D["ID"][uni_id.upper()], D["ID"][gn.upper()] = uni_id.upper(), uni_id.upper(), uni_id.upper()
        D["GN"][uni_ac], D["GN"][uni_id], D["GN"][gn] = gn.upper(), gn.upper(), gn.upper()
        D["GN"][uni_ac.upper()], D["GN"][uni_id.upper()], D["GN"][gn.upper()] = gn.upper(), gn.upper(), gn.upper()
        D["seq"][uni_ac.upper()] = seq
        D["dc"][uni_ac.upper()] = dc
        D["des"][uni_ac.upper()] = des

        for ac in accs:
            D["AC"][ac] = uni_ac.upper()
            D["AC"][ac.upper()] = uni_ac.upper()

        for feat in record.features:
            if feat[0] == "DOMAIN":
                name = feat[-2].split(".")[0]
                if re.search("([\d]+)", str(feat[1])) and re.search("([\d]+)", str(feat[2])):
                    start = int(re.search("([\d]+)", str(feat[1])).group(1))
                    end = int(re.search("([\d]+)", str(feat[2])).group(1))
                    D["doms"][uni_ac][(start,end)] = name
    return D

def get_pfam_doms(pfam_file, prot_id, max_eval=1):
    """Pfam-A matches in species proteome. File downloaded from PFAM.
    """
    pfams = defaultdict(lambda: defaultdict(set) )
    pfam_sets = defaultdict(set)
    pfam_names = {}

    with open_file(pfam_file) as f:
        for line in f:
            if line[0] != "#":
                t = line.rstrip().split("\t")
                uni_ac = t[0].upper()
                start, end = int(t[3]), int(t[4])
                pfam_ac, pfam_name, domain_e_val = t[5], t[6], float(t[12])

                if uni_ac in prot_id and domain_e_val <= max_eval:
                    uni_ac = prot_id[uni_ac]
                    pfams[uni_ac][pfam_ac].add((start, end, domain_e_val))
                    pfam_sets[pfam_ac].add(uni_ac)
                    pfam_names[pfam_ac] = pfam_name

    return pfams, pfam_sets, pfam_names

# def get_pfam_doms(pfam_hits_file, prot_id, max_eval=1):
#     """ Reads pre-generated file with domains already filtered
#         (best E-value and avoiding overlapping)
#     """
#     pfams = defaultdict(lambda: defaultdict(set) )
#     pfam_sets = defaultdict(set)
#     pfam_names = {}
#     with open_file(pfam_hits_file) as f:
#         for line in f:
#             if line[0] != "#" and line.strip():
#                 t = line.rstrip().split("\t")
#                 uni_ac = t[0].split("|")[0].upper()
#                 pfam_ac, pfam_name = t[1].split("|")
#                 start, end, e_val = int(t[2]), int(t[3]), float(t[4])
#
#                 if e_val <= max_eval and uni_ac in prot_id:
#                     uni_ac = prot_id[uni_ac]
#                     pfams[uni_ac][pfam_ac].add((start, end, e_val))
#                     pfam_sets[pfam_ac].add(uni_ac)
#                     pfam_names[pfam_ac] = pfam_name
#     return pfams, pfam_sets, pfam_names

def pfam_descriptions(pfam_data_file):
    d = {}
    with open_file(pfam_data_file) as f:
        for line in f:
            if "#=GF AC" in line:
                ac = line.rstrip().split()[2].split(".")[0]
            elif "#=GF DE" in line:
                des = " ".join(line.rstrip().split()[2:])
            elif line == "//\n":
                d[ac] = des
    return d

def get_linear_motifs(elm_hits_file, prot_id, max_eval=1):
    """ Reads pre-generated file with ELMs annotated for each protein
        If there's any filtering for these, it has already been done. All the
        info of this file is considered equally valid (take all)
    """
    elms = defaultdict(lambda: defaultdict(set) )
    elm_sets = defaultdict(set)
    elm_names = {}
    with open_file(elm_hits_file) as f:
        for line in f:
            if line[0] != "#" and line.strip():
                t = line.rstrip().split("\t")
                uni_ac = t[0].split("|")[0].upper()
                elm_ac, elm_name = t[1].split("|")
                elm_start, elm_end = int(t[2]), int(t[3])
                some_score = float(t[4])

                if some_score <= max_eval and uni_ac in prot_id:
                    uni_ac = prot_id[uni_ac]
                    elms[uni_ac][elm_ac].add((elm_start, elm_end, some_score))
                    elm_sets[elm_ac].add(uni_ac)
                    elm_names[elm_ac] = elm_name
    return elms, elm_sets, elm_names

def get_ptms(psp_file, prot_id):
    psp = defaultdict(set)
    done_pos = defaultdict(set)
    with open_file(psp_file) as f:
        for line in f:
            t = line.rstrip().split("\t")
            uni_ac = t[2]
            if uni_ac in prot_id:
                uni_ac = prot_id[uni_ac]
                pos, mod = t[4].split("-")
                res, pos = re.search("(\w)([\d]+)", pos).group(1,2)
                lt_lit="."
                ms_lit="."
                ms_cst="."
                cst_cat="."
                if len(t)>10 and hasNumbers(t[10]):
                    lt_lit = assign_Number(t[10])
                if len(t)>11 and hasNumbers(t[11]):
                    ms_lit = assign_Number(t[11])
                if len(t)>12 and hasNumbers(t[12]):
                    ms_cst = assign_Number(t[12])
                if len(t)>13 and hasNumbers(t[13]):
                    cst_cat = assign_Number(t[13])

                if pos not in done_pos[(uni_ac,mod)]:
                    psp[uni_ac].add( (res, pos, mod) )
                    done_pos[(uni_ac,mod)].add(pos)
    return psp

def get_ppi_pairs(biogrid_file):
	with open_file(biogrid_file) as f:
		for line in f:
			t = line.rstrip().split("\t")
			if line[0] != "#":
				a = [t[7]] + t[9].split("|")
				b = [t[8]] + t[10].split("|")

	return

biogrid_file = sp_data_dir + "BIOGRID-ORGANISM-3.5.165.tab2.txt.gz"
get_ppi_pairs(biogrid_file)
sys.exit
####
error_msg = """
ERROR
Usage:
~$ python {} 'mode'

where mode is either 'normal' or 'mongo'
""".format(sys.argv[0])

if len(sys.argv) < 2:
    sys.exit(error_msg)
mode = sys.argv[1]
if mode not in ["mongo", "normal"]:
    sys.exit(error_msg)


## 1. Get Uniprot Accession, IDs and Gene names conversions ###
# prot_id_dict, sequences = parse_fasta(proteome_file)
prot_dict = get_protein_data_from_uniprot_text(uniprot_text_file)

## 2. Get Pfam domains and ELMs
pfams, pfam_sets, pfam_names = get_pfam_doms(pfam_matches_file,
                                                prot_dict["AC"])

pfam_des = pfam_descriptions(pfam_data_file)

elms, elm_sets, elm_names = get_linear_motifs(elm_hits_file,
                                            prot_dict["AC"])
## 3. Get PTMs
ptms = get_ptms(psp_file, prot_dict["AC"])

## 4. Output
pp = pprint.PrettyPrinter(indent=4)
outfile_name = sp_data_dir+"protein_data_"+species+"_"+mode+".json.gz"
if mode == "normal":
    protein_data = {}
    for uni_ac in prot_dict["seq"]:
        protein_data[uni_ac] = {
                "uniprot_id" : prot_dict["ID"][uni_ac],
                "gene" : prot_dict["GN"][uni_ac],
                "description" : prot_dict["des"][uni_ac],
                "data_class" : prot_dict["dc"][uni_ac],
                "length" : len(prot_dict["seq"][uni_ac]),
                "sequence" : prot_dict["seq"][uni_ac],
                "pfams" : {},
                "elms": {},
                "phosphorylation" : {},
                "acetylation" : {},
                "cosmic_muts" : []
        }

        if uni_ac in pfams:
            for pfam_ac in sorted(pfams[uni_ac]):
                des = ""
                if pfam_ac.split(".")[0] in pfam_des:
                    des = pfam_des[pfam_ac.split(".")[0]]
                pfam_name = pfam_names[pfam_ac]
                if pfam_name not in protein_data[uni_ac]["pfams"]:
                    protein_data[uni_ac]["pfams"][pfam_name] = []
                for (start, end, e_val) in sorted(pfams[uni_ac][pfam_ac],
                                                key=lambda x: int(x[0])):
                    protein_data[uni_ac]["pfams"][pfam_name].append(
                        {
                            "acc" : pfam_ac,
                            "des": des,
                            "start" : int(start),
                            "end" : int(end),
                            "e-val" : float(e_val)
                        })

        if uni_ac in elms:
            for elm_ac in sorted(elms[uni_ac]):
                elm_name = elm_names[elm_ac]
                if elm_name not in protein_data[uni_ac]["elms"]:
                    protein_data[uni_ac]["elms"][elm_name] = []
                for (start, end, score) in sorted(elms[uni_ac][elm_ac],
                                                key=lambda x: int(x[0])):
                    protein_data[uni_ac]["elms"][elm_name].append(
                        {
                            "acc" : elm_ac,
                            "start" : int(start),
                            "end" : int(end),
                            "score" : float(score)
                        })

        if uni_ac in ptms:
            for (res, pos, mod) in sorted(ptms[uni_ac], key=lambda x: int(x[1])):
                if mod == "p":
                    ptm_type = "phosphorylation"
                elif mod == "ac":
                    ptm_type = "acetylation"
                if pos not in protein_data[uni_ac][ptm_type]:
                    protein_data[uni_ac][ptm_type][int(pos)] = res

    with open_file(outfile_name, "w") as out:
        json.dump(protein_data, out)

elif mode == "mongo":
    protein_data = {}
    with open_file(outfile_name, "w") as out:
        for uni_ac in prot_dict["seq"]:
            protein_data = {
                    "uniprot_acc" : uni_ac,
                    "uniprot_id" : prot_dict["ID"][uni_ac],
                    "gene" : prot_dict["GN"][uni_ac],
                    "description" : prot_dict["des"][uni_ac],
                    "data_class" : prot_dict["dc"][uni_ac],
                    "length" : len(prot_dict["seq"][uni_ac]),
                    "sequence" : prot_dict["seq"][uni_ac],
                    "pfams" : [],
                    "elms" : [],
                    "phosphorylation" : [],
                    "acetylation" : [],
                    "cosmic_muts" : []
            }

            if uni_ac in pfams:
                for pfam_ac in sorted(pfams[uni_ac]):
                    des = ""
                    if pfam_ac.split(".")[0] in pfam_des:
                        des = pfam_des[pfam_ac.split(".")[0]]
                    pfam_name = pfam_names[pfam_ac]
                    for (start, end, e_val) in sorted(pfams[uni_ac][pfam_ac],
                                                      key=lambda x: int(x[0])):
                        protein_data["pfams"].append(
                            {
                                "acc" : pfam_ac,
                                "des": des,
                                "name" : pfam_name,
                                "start" : int(start),
                                "end" : int(end),
                                "e-val" : float(e_val)
                            })

            if uni_ac in elms:
                for elm_ac in sorted(elms[uni_ac]):
                    elm_name = elm_names[elm_ac]
                    for (start, end, score) in sorted(elms[uni_ac][elm_ac],
                                                    key=lambda x: int(x[0])):
                        protein_data["elms"].append(
                            {
                                "acc" : elm_ac,
                                "name" : elm_name,
                                "start" : int(start),
                                "end" : int(end),
                                "seq" : prot_dict["seq"][uni_ac][start-1:end],
                                "score" : float(score)
                            })

            if uni_ac in ptms:
                for (res, pos, mod) in sorted(ptms[uni_ac],
                                             key=lambda x: int(x[1])):
                    if mod == "p":
                        ptm_type = "phosphorylation"
                    elif mod == "ac":
                        ptm_type = "acetylation"
                    protein_data[ptm_type].append(
                        {
                            "pos" : int(pos),
                            "res" : res
                        })

            out.write(str(protein_data)+"\n")

sys.exit()
