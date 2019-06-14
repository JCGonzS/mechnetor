#!/usr/bin/env python

import sys, os, re, math, gzip
from collections import defaultdict

#Usage : ~$ python find_all_elms_fasta.py fasta_file

def open_file(input_file, mode="r"):
    """ Open file Zipped or not
    """

    if re.search(".gz$", input_file):
        infile = gzip.open(input_file, mode)
    else:
        infile = open(input_file, mode)
    return infile

def main(fasta_file, elm_classes_file, print_out=True, outfile="elm_hits.tsv.gz"):
    fpm2_script = "fpm2.pl"
    to_print = []
    elms = defaultdict(lambda: defaultdict(list))
    with open_file(elm_classes_file) as f:
        for line in f:
            if line[0]!="#" and not line.startswith("Accession"):
                elm_acc, elm_ide = line.split("\t")[:2]
                elm_regex, elm_prob = line.split("\t")[4:6]

                mode = "cat "
                if ".gz" in fasta_file:
                    mode = "zcat "
                os.system( mode+fasta_file+" | ./"+fpm2_script+" \""+elm_regex+"\" -m 1 > tmp.txt")

                with open_file("tmp.txt") as f2:
                    for line2 in f2:
                        if line2[0]==">":
                            label = line2.split()[0].split("/")[0].replace(">","")
                            start, end = line2.split()[0].split("/")[1].split("-")
                            gene = "-"
                            if "GN=" in line2.split("/")[1]:
                                gene = re.search("GN=([^\s]+)", line2.split("/")[1]).group(1)
                        else:
                            seq = line2.rstrip()
                            to_print.append( "\t".join([elm_ide, elm_acc, label, gene, start, end, elm_prob, "-", "!", seq]) )

                            elms[label][elm_acc].append(
                                {"start" : int(start),
                                 "end" :   int(end),
                                 "seq" :   seq,}
                                )
        os.unlink("tmp.txt")

    if print_out:
        with open_file(outfile, "w") as out:
            for ele in to_print:
                out.write(ele+"\n")
    else:
        return elms


if __name__ == "__main__":
    main(sys.argv[1],
        "static/data/common/elm_classes_May2019.tsv.gz",
        print_out=True, outfile="test.gz")
    sys.exit()
