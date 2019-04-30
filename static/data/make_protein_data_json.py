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
pfam_data_file = common_data_dir + "Pfam-A.hmm.dat.gz"
species = "Hsa"
sp_data_dir = data_dir+"species/"+species+"/"
uniprot_text_file = sp_data_dir + "uniprot_homo_sapiens_proteome_73928prts_Mar2019_data.txt.gz"
pfam_matches_file = sp_data_dir + "pfamA_matches_9606_Aug18.tsv.gz"
elm_hits_file =     sp_data_dir + "elm_hits.tsv.gz"
# elm_hits_file =     sp_data_dir + "elm_parsed_info_Overlapping.tsv.gz"
elm_dom_file = common_data_dir +  "elm_interaction_domains_edited_Jan18.tsv"
psp_file =          sp_data_dir + "PSP_ptms_human_Mar2019.tsv.gz"


def open_file(input_file, mode="r"):
    """ Open file zipped or not
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

def merge_gene_names(genes):
    """For a list of similar gene names ("GENE1A", "GENE2B", "GENE3C"),
    returns a string following the format: "GENE(1A,2B,3C)"
    """
    pattern = ""
    flag = 0
    for i in range(1, len(genes[0])+1):
        reg = genes[0][:i]
        for g in genes[1:]:
            if reg != g[:i]:
                flag = 1
        if flag == 0:
            pattern = reg
        else:
            break
    var = []
    for g in genes:
        v = g.split(pattern)[1]
        if v.strip():
            var.append(v)

    return pattern+"("+",".join(var)+")"


def get_protein_data_from_uniprot_text(uniprot_file):
    """Extracts information from a UniProt text file.
    Requires Bio::SwissProt module
    """

    D = defaultdict(lambda: defaultdict(set))
    masks = {}
    for record in SwissProt.parse(open_file(uniprot_file)):
        dc  = record.data_class
        uni_id = record.entry_name
        accs = record.accessions
        uni_ac = accs[0].upper()
        des = re.search("Name: Full=([^;|{]+)", record.description).group(1)
        gns, syns = [uni_id], []
        if record.gene_name.strip():
            names = [match.split()[0] for match in re.findall("Name=([^;]+)", record.gene_name)]
            if len(names)>0:
                gns = []
                for name in names:
                    name = name.split("_")[0]
                    if name not in gns:
                        gns.append(name)

            syns = [match.split()[0].strip(",") for match in re.findall("Synonyms=([^;]+)", record.gene_name)]

        if len(gns)>1:
            gn = merge_gene_names(gns)
        else:
            gn = gns[0]

        # if dc != "Reviewed":
        #     continue
        # if gn in D["AC"]:
        #     print uni_ac, gn, D["AC"][gn]
        # else:
        #     D["AC"][gn] = uni_ac.upper()
        for dic, val in zip(["AC", "ID", "GN"], [uni_ac, uni_id, gn]):
            for key in [uni_ac, uni_id, gn]:
                D[dic][key].add(val.upper())
                D[dic][key.upper()].add(val.upper())

        for ac in accs[1:]:
            # if ac!=uni_ac and ac in D["AC"]:
            #     print uni_ac, ac, D["AC"][ac]
            D["AC"][ac].add(uni_ac.upper())
            D["AC"][ac.upper()].add(uni_ac.upper())
        for g in gns[1:]:
            D["AC"][g].add(uni_ac.upper())
            D["AC"][g.upper()].add(uni_ac.upper())

        D["des"][uni_ac.upper()] = des
        D["dc"][uni_ac.upper()] = dc
        D["seq"][uni_ac.upper()] = record.sequence
        masks[uni_ac] = "0" * len(record.sequence)

    return D, masks

def get_pfam_doms(pfam_file, prot_dict, max_eval=999):
    """Pfam-A matches in species proteome. File downloaded from PFAM.
    """
    pfams = defaultdict(lambda: defaultdict(set) )
    pfams_temp = defaultdict(lambda: defaultdict(set) )
    pfam_sets = defaultdict(set)
    pfam_names = {}
    with open_file(pfam_file) as f:
        for line in f:
            if line[0] != "#":
                t = line.rstrip().split("\t")
                uni_ac = t[0].upper()
                start, end = int(t[4]), int(t[5])
                pfam_ac, pfam_name, domain_e_val = t[6], t[7], float(t[13])

                # if domain_e_val <= max_eval: ## No e-value cut-off. Just take what is annotated by Pfam

                if uni_ac in prot_dict["seq"]: ## Keep annotation for primary UniProt accessions
                    pfams[uni_ac][pfam_ac].add((start, end, domain_e_val))
                    pfam_sets[pfam_ac].add(uni_ac)
                    pfam_names[pfam_ac] = pfam_name
                elif uni_ac in prot_dict["AC"]: ## Save annotation for secondary accessions in a different dictionary
                    pfams_temp[uni_ac][pfam_ac].add((start, end, domain_e_val))

    ## Transfer annotation to those main accessions which were not annotated
    for alt_ac in pfams_temp:
        for uni_ac in prot_dict["AC"][alt_ac]:
            if uni_ac not in pfams: ## If the associated main accession is not annotated already:
                pfams[uni_ac] = pfams_temp[alt_ac]
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

def get_interacting_linear_motifs(elm_dom_file, elm_hits_file, prot_dict, masks,
                                  max_overlap, max_eval=1):
    """ Reads pre-generated file with ELMs annotated for each protein
        If there's any filtering for these, it has already been done. All the
        info of this file is considered equally valid (take all)
    """
    # Restrict to elms apearing in the elm_dom interaction file. PIV won't use any other
    int_elms = []
    with open_file(elm_dom_file) as f:
        for line in f:
            t = line.rstrip().split("\t")
            if line[0]!="#":
                int_elms.append(t[0])

    elms = defaultdict(lambda: defaultdict(set) )
    elms_temp = defaultdict(lambda: defaultdict(set) )
    elm_names = {}
    with open_file(elm_hits_file) as f:
        for line in f:
            if line[0] != "#" and line.strip():
                t = line.rstrip().split("\t")
                elm_name, elm_ac = t[0], t[1]
                uni_ac, uni_id = t[2].split("|")[1:]
                elm_start, elm_end = int(t[4]), int(t[5])
                prob_score = float(t[6]) # probability score based on the combined expected frequencies of the AAs in the regular expression

                if elm_name in int_elms:
                    if uni_ac in prot_dict["seq"]:
                        overlap = calculate_overlap(elm_start, elm_end, masks[uni_ac])
                        if overlap <= max_overlap:
                            masks[uni_ac] = fill_mask(elm_start, elm_end, masks[uni_ac])
                            elms[uni_ac][elm_ac].add((elm_start, elm_end, prob_score))
                            elm_names[elm_ac] = elm_name
                    elif uni_ac in prot_dict["AC"]:
                        elms_temp[uni_ac][elm_ac].add((elm_start, elm_end, prob_score))

    for alt_ac in elms_temp:
        for uni_ac in prot_dict["AC"][alt_ac]:
            if uni_ac not in elms:
                overlap = calculate_overlap(elm_start, elm_end, masks[uni_ac])
                if overlap <= max_overlap:
                    masks[uni_ac] = fill_mask(elm_start, elm_end, masks[uni_ac])
                    elms[uni_ac] = elms_temp[alt_ac]
                    elm_names[elm_ac] = elm_name
    return elms, elm_names

def get_ptms(psp_file, prot_id):
    psp = defaultdict(set)
    psp_temp = defaultdict(set)
    with open_file(psp_file) as f:
        for line in f:
            t = line.rstrip().split("\t")
            pos, mod = t[4].split("-")
            res, pos = re.search("(\w)([\d]+)", pos).group(1,2)
            lt_lit="."
            ms_lit="."
            ms_cst="."
            cst_cat="."
            uni_ac = t[2]
            if len(t)>10 and hasNumbers(t[10]):
                lt_lit = assign_Number(t[10])
            if len(t)>11 and hasNumbers(t[11]):
                ms_lit = assign_Number(t[11])
            if len(t)>12 and hasNumbers(t[12]):
                ms_cst = assign_Number(t[12])
            if len(t)>13 and hasNumbers(t[13]):
                cst_cat = assign_Number(t[13])

            if uni_ac in prot_id["seq"]:
                seq = prot_id["seq"][uni_ac]
                if len(seq)>=int(pos) and seq[int(pos)-1]==res:
                    psp[uni_ac].add( (res, pos, mod) )
            elif uni_ac in prot_id["AC"]:
                psp_temp[uni_ac].add( (res, pos, mod) )


    for alt_ac in psp_temp:
        for uni_ac in prot_dict["AC"][alt_ac]:
            if uni_ac not in psp:
                seq = prot_id["seq"][uni_ac]
                for (res,pos,mod) in psp_temp[alt_ac]:
                    if len(seq)>=int(pos) and seq[int(pos)-1]==res:
                        psp[uni_ac].add( (res, pos, mod) )

    return psp

def get_ppi_pairs(biogrid_file):
	with open_file(biogrid_file) as f:
		for line in f:
			t = line.rstrip().split("\t")
			if line[0] != "#":
				a = [t[7]] + t[9].split("|")
				b = [t[8]] + t[10].split("|")

	return

def calculate_overlap(start, end, mask):
    length = end - start + 1
    sub_mask = mask[start-1:end]
    n1 = len(sub_mask) - len(sub_mask.replace("1", ""))
    overlap = n1 / float(length)
    return overlap

def fill_mask(start, end, mask):
    length = end - start + 1
    mask = mask[:start] + ("1" * length) + mask[end-1:]
    return mask

# biogrid_file = sp_data_dir + "BIOGRID-ORGANISM-3.5.165.tab2.txt.gz"
# get_ppi_pairs(biogrid_file)
# sys.exit

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
prot_dict, masks = get_protein_data_from_uniprot_text(uniprot_text_file)

## 2. Get Pfam domains and ELMs
pfams, pfam_sets, pfam_names = get_pfam_doms(pfam_matches_file, prot_dict)
pfam_des = pfam_descriptions(pfam_data_file)
elms, elm_names = get_interacting_linear_motifs(elm_dom_file, elm_hits_file,
                                                          prot_dict, masks,
                                                          max_overlap=2)

## 3. Get PTMs
ptms = get_ptms(psp_file, prot_dict)

## 4. Output
pp = pprint.PrettyPrinter(indent=4)
outfile_name = sp_data_dir+"protein_data_"+species+"_"+mode+".json.gz"

if mode == "json":
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
            seq = prot_dict["seq"][uni_ac]
            protein_data = {
                    "uniprot_acc" :     uni_ac,
                    "uniprot_id" :      list(prot_dict["ID"][uni_ac])[0], # always, only 1
                    "gene" :            list(prot_dict["GN"][uni_ac])[0],
                    "description" :     prot_dict["des"][uni_ac],
                    "data_class" :      prot_dict["dc"][uni_ac],
                    "length" :          len(seq),
                    "sequence" :        seq,
                    "pfams" :           [],
                    "elms" :            [],
                    "phosphorylation" : [],
                    "acetylation" :     [],
                    "cosmic_muts" :     []
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
                                "acc" :   pfam_ac,
                                "des":    des,
                                "name" :  pfam_name,
                                "start" : int(start),
                                "end" :   int(end),
                                "e-val" : float(e_val)
                            })

            if uni_ac in elms:
                for elm_ac in sorted(elms[uni_ac]):
                    elm_name = elm_names[elm_ac]
                    for (start, end, score) in sorted(elms[uni_ac][elm_ac],
                                                    key=lambda x: int(x[0])):
                        protein_data["elms"].append(
                            {
                                "acc" :   elm_ac,
                                "name" :  elm_name,
                                "start" : int(start),
                                "end" :   int(end),
                                "seq" :   seq[start-1:end],
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
                            "res" : res,
                            "source": "PSP"
                        })
            else:
                for i in range(9):
                    iso_uni_ac = uni_ac+"-"+str(i)
                    if iso_uni_ac in ptms:
                        for (res, pos, mod) in sorted(ptms[iso_uni_ac],
                                                     key=lambda x: int(x[1])):
                            if mod == "p":
                                ptm_type = "phosphorylation"
                            elif mod == "ac":
                                ptm_type = "acetylation"
                            protein_data[ptm_type].append(
                                {
                                    "pos" : int(pos),
                                    "res" : res,
                                    "source": "PSP,-"+str(i)
                                })
                        break

            out.write(str(protein_data)+"\n")

    sys.exit()
