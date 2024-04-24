#!/usr/bin/env python

import sys, os, re, math, gzip, random, string
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

def main(fasta_file, elm_classes, print_out=True,
         outfile="elm_hits.tsv.gz", tmp_file="/tmp/lm_search.txt"):
    fpm2_script = "fpm2.pl"
    to_print = []
    elms = defaultdict(lambda: defaultdict(list))
 
    for elm_acc in elm_classes:
        elm_ide = str(elm_classes[elm_acc]["ide"])
        elm_regex = str(elm_classes[elm_acc]["regex"])
        elm_prob = str(elm_classes[elm_acc]["prob"])

        mode = "cat "
        if ".gz" in fasta_file:
            mode = "zcat "
        com = "{} {} | perl {} \"{}\" -m 1 > {}".format(mode, fasta_file, fpm2_script, elm_regex, tmp_file)
        os.system(com)

        with open_file(tmp_file, "rt") as f:
            for line in f:
                if line[0]==">":
                    label = line.split()[0].split("/")[0].replace(">","")
                    start, end = line.split()[0].split("/")[1].split("-")
                    gene = "-"
                    if "GN=" in line.split("/")[1]:
                        gene = re.search("GN=([^\s]+)", line.split("/")[1]).group(1)
                else:
                    seq = line.rstrip()
                    to_print.append( "\t".join([elm_ide, elm_acc, label, gene, start, end, elm_prob, "-", "!", seq]) )

                    elms[label][elm_acc].append(
                        {   "start" : int(start),
                            "end" :   int(end),
                            "seq" :   seq,}
                        )

    os.unlink(tmp_file)

    if print_out:
        with open_file(outfile, "wt") as out:
            for ele in to_print:
                out.write(ele+"\n")
    else:
        return elms


if __name__ == "__main__":
    main(sys.argv[1],
        "static/data/common/elm_classes_Oct2019.tsv",
        print_out=True, outfile="test.gz")
    sys.exit()
