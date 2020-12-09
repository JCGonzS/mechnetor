import sys, os, re, gzip, itertools, math, ast, csv
from Bio import SwissProt, SeqIO
from collections import defaultdict

sp = sys.argv[1]
data_dir = "static/data/species/"+sp+"/"

# Read fasta
gene = {}
with gzip.open(data_dir+"uniprot_Nov2019_proteome_"+sp+".fasta.gz") as f:
    for line in f:
        if line[0]==">":
            ac = line.rstrip().split()[0].split("|")[1]
            g = "-"
            if "GN=" in line:
                g = re.search("GN=([^\s]+)", line.rstrip()).group(1)
            gene[ac] = g

# Read pfam matches file
with gzip.open(data_dir+"pfamA_r32.0_matches_"+sp+"_GN.tsv.gz" ,"w") as out:
    with gzip.open(data_dir+"pfamA_r32.0_matches_"+sp+".tsv.gz") as f:
        for line in f:
            if line[0]!="#":
                t = line.rstrip().split("\t")
                ac = t[0]
                g = "-"
                if ac in gene:
                    g = gene[ac]

                out.write( "\t".join([ac, g]+t[1:])+"\n" ) 
            else:
                out.write( line ) 