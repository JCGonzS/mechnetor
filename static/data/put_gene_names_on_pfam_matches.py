#!/usr/bin/env python

import sys, re
import json, csv, gzip, pprint
from collections import defaultdict
from Bio import SwissProt


def open_file(input_file, mode="r"):
    """ Open file Zipped or not
    """
    if re.search(".gz$", input_file):
        infile = gzip.open(input_file, mode)
    else:
        infile = open(input_file, mode)
    return infile


def get_protein_data_from_uniprot_text(uniprot_file):
    """From SwissProt"""

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


    return D

sp_data_dir = "species/Hsa/"
uniprot_text_file = sp_data_dir+"uniprot_homo_sapiens_proteome_73112prts_Aug2018_data.txt.gz"
pfam_file = sp_data_dir+"pfamA_matches_9606.tsv.gz"

prot_dict = get_protein_data_from_uniprot_text(uniprot_text_file)

with open_file(pfam_file) as f:
    for line in f:
        t = line.rstrip().split("\t")
        uni_ac = t[0]
        gene = "-"
        if uni_ac in prot_dict["GN"]:
            gene = prot_dict["GN"][uni_ac]
        print uni_ac+"\t"+gene+"\t"+"\t".join(t[1:])
