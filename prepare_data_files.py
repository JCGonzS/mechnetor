#!/usr/bin/env python

import sys, re
import gzip
from collections import defaultdict


def open_file(input_file, mode="r"):
    """ Open file Zipped or not
    """
    if re.search(".gz$", input_file):
        infile = gzip.open(input_file, mode)
    else:
        infile = open(input_file, mode)
    return infile

def parse_fasta(proteome_file):
    """Get protein sequences from FASTA file
    and the convertion between the protein IDs (accession,gene name, uniprot ID)
    """

    D = defaultdict(dict)
    seqs, masks = {}, {}

    with open_file(proteome_file) as f:
        for line in f:
            if line[0] == ">":
                uni_ac, uni_id = line.split()[0].split("|")[1:]
                gn = uni_id
                if "GN=" in line:
                    gn = re.search(r"GN=([^\s]*)\s",line).group(1)

                D["AC"][uni_ac], D["AC"][uni_id], D["AC"][gn] = uni_ac, uni_ac, uni_ac
                D["ID"][uni_ac], D["ID"][uni_id], D["ID"][gn] = uni_id, uni_id, uni_id
                D["GN"][uni_ac], D["GN"][uni_id], D["GN"][gn] = gn, gn, gn
                D["AC"][uni_ac.upper()] = uni_ac
                D["ID"][uni_ac.upper()] = uni_id
                D["GN"][uni_ac.upper()] = gn

                seqs[uni_ac] = ""
                masks[uni_ac] = ""

            else:
                seqs[uni_ac] += line.rstrip()
                masks[uni_ac] += "0" * len(line.rstrip())

    return D, seqs, masks

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

def parse_pfam_doms(pfam_hits_file, prot_dict, masks,
                    max_eval=0.1, max_overlap=0.2):
    """ Reads hmmsearch v. Pfam database results to extract the Pfam domains
        found in each protein

    pfam_hits_file = "uniprot_v_Pfam_hmmsearch_sum.txt.gz"
    ac_dict = dictionary with conversion of any protein ID to UniProt accession
    sequences = dictionary containing the protein sequences per UniProt accession
    masks = = dictionary containing the sequence mask per UniProt accession
    max_eval = maximum E-value allowed to accept the match.
               '0.1' is the maximum value found in this file anyway.
    max_overlap = maximum overlap allowed between a new domain and previously
                    accepted ones
    """

    domains = defaultdict(lambda: defaultdict(set) )
    domain_sets = defaultdict(set)
    with open_file(pfam_hits_file) as f:
        #: The entries of this file are already sorted by E-value (?), in ascending order,
        #: this means we can use a "first come, first served" to get the domains for each protein
        #:
        #: Lines look like this:
        #: PF12895.4        ANAPC3         sp|Q9W040|OSM1_DROME    -1      1239    1308    2.9     6.9e+02 -1.8    !
        for line in f:
            if line[0] != "#" and line.strip():
                t = line.rstrip().split("\t")
                if t[-1] == "!":
                    dom_id, dom_name = t[0], t[1]
                    uni_ac = t[2].split("|")[1]
                    dom_start, dom_end = int(t[4]), int(t[5])
                    e_val, some_score = float(t[6]), float(t[7]) ## Not sure about these 2 values
                    if (e_val <= max_eval and uni_ac in prot_dict["AC"]):
                        uni_ac = prot_dict["AC"][uni_ac]
                        gene = prot_dict["GN"][uni_ac]

                        #: determine the overlap with previously accepted domains
                        overlap = calculate_overlap(dom_start, dom_end,
                                                    masks[uni_ac])

                        if overlap <= max_overlap:
                            domains[uni_ac][e_val].add("\t".join([uni_ac+"|"+gene, dom_id+"|"+dom_name, str(dom_start), str(dom_end), str(e_val)]))
                            domain_sets[dom_name].add(uni_ac)
                            #: if domain is accepted, its region in the
                            #: sequence mask is filled with "1"s
                            masks[uni_ac] = fill_mask(dom_start, dom_end,
                                                        masks[uni_ac])


    return domains, domain_sets, masks


def parse_linear_motifs(elm_hits_file, prot_dict, masks,
                        max_eval=10, max_overlap=0.2):
    """ Reads hmmsearch v. ELM database results to anotate linear motifs

    elm_hits_file = "elm_hits.tsv.gz"
    input_proteins = list of protein accessions given to the program as input
    ac_dict = dictionary with conversion of any protein ID to UniProt accession
    sequences = dictionary containing the protein sequences per UniProt accession
    masks = dictionary containing the sequence mask per UniProt accession
    max_eval = maximum E-value allowed to accept the match.
    max_overlap = maximum overlap allowed between a new domain and previously
                    accepted ones
    """
    lms = defaultdict(lambda: defaultdict(set) )
    lm_sets = defaultdict(set)
    with open_file(elm_hits_file) as f:
        #: Line looks like this:
        #: CLV_C14_Caspase3-7      ELME000321      sp|P30447|1A23_HUMAN    HLA-A   187     191     0.0030937403    -       !       TCVDG
        for line in f:
            if line[0] != "#" and line.strip():
                t = line.rstrip().split("\t")
                elm_name, elm_id = t[0], t[1]
                uni_ac, uni_id = t[2].split("|")[1:]
                elm_start, elm_end = int(t[4]), int(t[5])
                some_score = float(t[6]) # conservation score?

                if uni_ac in prot_dict["AC"]:
                    uni_ac = prot_dict["AC"][uni_ac]
                    gene = prot_dict["GN"][uni_ac]

                    #: determine the overlap with previously accepted LMs
                    # overlap = calculate_overlap(elm_start, elm_end, masks[uni_ac])

                    # if overlap <= max_overlap:
                    lms[uni_ac][elm_id].add("\t".join([ uni_ac+"|"+gene, elm_id+"|"+elm_name, str(elm_start), str(elm_end), str(some_score)]))
                    lm_sets[elm_name].add(uni_ac)
                        #: if the LM is accepted, its region in the
                        #: sequence mask is filled with "1"s
                        # masks[uni_ac] = fill_mask(elm_start, elm_end,
                        #                             masks[uni_ac])
        return lms, lm_sets

def main( data_dir="static/data/",
          species="Hsa",
          max_overlap=0.2,
          no_overlap_between_doms_and_lms="y"
        ):

    # General files:
    

    #: Species Files
    sp_data_dir = data_dir+"species/"+species+"/"
    proteome_file = sp_data_dir + "uniprot_sprot_species.fasta.gz"
    pfam_hits_file = sp_data_dir + "uniprot_v_Pfam_hmmsearch_sum.txt.gz"
    elm_hits_file = sp_data_dir + "elm_hits.tsv.gz"

    #: Output Files
    pfam_parsed_file = sp_data_dir + "pfam_parsed_info.tsv.gz"
    elm_parsed_file = sp_data_dir + "elm_parsed_info_noOverlap.tsv.gz"

    #: Get protein ID-dictionary, sequences and masks
    prot_id_dict, sequences, masks = parse_fasta(proteome_file)

    #: Get Pfam domains per protein from hmmsearch result file
    domains, domain_sets, ms = parse_pfam_doms(pfam_hits_file, prot_id_dict,
                                       masks, max_overlap=max_overlap)

    #: Print them in TSV file
    with gzip.open(pfam_parsed_file, "wb") as out:
        c = ["#UniProt_acc|Gene", "Pfam_ID|Name", "Start", "End", "E-val"]
        out.write("\t".join(c) + "\n")
        for uni_ac in sorted(domains):
            for e_val in sorted(domains[uni_ac]):
                for dom in domains[uni_ac][e_val]:
                    out.write(dom + "\n")

    if no_overlap_between_doms_and_lms == "y":
        masks = ms

    #: Get Linear Motifs per protein from hmmsearch result file
    lms, lm_sets = parse_linear_motifs(elm_hits_file, prot_id_dict, masks,
                                    max_overlap=max_overlap)
    #: Print them in TSV file
    with gzip.open(elm_parsed_file, "wb") as out:
        c = ["#UniProt_acc|Gene", "ELM_ID|Name", "Start", "End", "Some_Score"]
        out.write("\t".join(c) + "\n")
        for uni_ac in sorted(lms):
            for lm in sorted(lms[uni_ac]):
                for instance in lms[uni_ac][lm]:
                    out.write(instance + "\n")
    return

##USE dom2dom.py FOR CALCULATING DOM-DOM PROPENSITIES!!!!!!!!!

if __name__ == "__main__":
    main()
