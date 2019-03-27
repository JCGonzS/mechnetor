#!/usr/bin/env python

import sys, re, os
import gzip, math
import os.path
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

def parse_uniprot_data(data_file, db="sp,tr"):
    D = defaultdict(dict)
    sequences, masks = {}, {}
    ref_proteome = set()

    for record in SwissProt.parse(open_file(data_file)):
        uni_id = record.entry_name
        accs = record.accessions
        uni_ac = accs[0]
        gn = record.gene_name
        if gn != "":
            gn = re.search("Name[s]?=([^;]+)", record.gene_name).group(1).split()[0]
        seq = record.sequence
        dc = record.data_class
        if dc == "Reviewed":
            dc = "sp"
        elif dc == "Unreviewed":
            dc = "tr"

        if dc in db.split(","):
            D["AC"][uni_ac], D["AC"][uni_id], D["AC"][gn] = uni_ac.upper(), uni_ac.upper(), uni_ac.upper()
            D["AC"][uni_ac.upper()], D["AC"][uni_id.upper()], D["AC"][gn.upper()] = uni_ac.upper(), uni_ac.upper(), uni_ac.upper()
            D["ID"][uni_ac], D["ID"][uni_id], D["ID"][gn] = uni_id.upper(), uni_id.upper(), uni_id.upper()
            D["ID"][uni_ac.upper()], D["ID"][uni_id.upper()], D["ID"][gn.upper()] = uni_id.upper(), uni_id.upper(), uni_id.upper()
            D["GN"][uni_ac], D["GN"][uni_id], D["GN"][gn] = gn.upper(), gn.upper(), gn.upper()
            D["GN"][uni_ac.upper()], D["GN"][uni_id.upper()], D["GN"][gn.upper()] = gn.upper(), gn.upper(), gn.upper()
            for acc in accs:
                D["AC"][acc], D["ID"][acc], D["GN"][acc] = uni_ac.upper(), uni_id.upper(), gn.upper()
                D["AC"][acc.upper()], D["ID"][acc.upper()], D["GN"][acc.upper()] = uni_ac.upper(), uni_id.upper(), gn.upper()

            sequences[uni_ac.upper()] = seq
            masks[uni_ac.upper()] = "0" * len(seq)

        if dc == "sp":
            ref_proteome.add(uni_ac.upper())

    return D, sequences, masks, ref_proteome

def get_protein_data_from_uniprot_text(uniprot_file):
    """Extracts information from a UniProt text file.
    Requires Bio::SwissProt module
    """

    D = defaultdict(lambda: defaultdict(set))
    masks = {}
    ref_proteome = set()

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

        for dic, val in zip(["AC", "ID", "GN"], [uni_ac, uni_id, gn]):
            for key in [uni_ac, uni_id, gn]:
                D[dic][key].add(val.upper())
                D[dic][key.upper()].add(val.upper())

        D["des"][uni_ac] = des
        D["dc"][uni_ac]  = dc
        D["seq"][uni_ac] = record.sequence
        masks[uni_ac] = "0" * len(record.sequence)
        if dc == "Reviewed":
            ref_proteome.add(uni_ac)

        for ac in accs[1:]:
            D["AC"][ac].add(uni_ac.upper())
            D["AC"][ac.upper()].add(uni_ac.upper())
        for g in gns[1:]:
            D["AC"][g].add(uni_ac.upper())
            D["AC"][g.upper()].add(uni_ac.upper())

    return D, masks, ref_proteome

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

def parse_pfam_doms(pfam_hmm_file, prot_dict, masks,
                    max_eval=0.1, max_overlap=0.2):
    """ Reads hmmsearch v. Pfam database results to extract the Pfam domains
        found in each protein

    pfam_hmm_file = "uniprot_v_Pfam_hmmsearch_sum.txt.gz"
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
    with open_file(pfam_hmm_file) as f:
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
                    uni_ac = t[2].split("|")[1].upper()
                    dom_start, dom_end = int(t[4]), int(t[5])
                    e_val, some_score = float(t[6]), float(t[7]) ## Not sure about these 2 values
                    if (e_val <= max_eval and uni_ac in prot_dict["AC"]):
                        uni_ac = prot_dict["AC"][uni_ac]
                        gene = prot_dict["GN"][uni_ac]
                        if len(gene)>1:
                            print line
                        continue

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

                if uni_ac in prot_id:# and domain_e_val <= max_eval:
                    for main_ac in prot_id[uni_ac]:
                        pfams[main_ac][pfam_ac].add((start, end, domain_e_val))
                        pfam_sets[pfam_ac].add(main_ac)
                        pfam_names[pfam_ac] = pfam_name

    return pfams, pfam_sets, pfam_names


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
                    for main_ac in prot_dict["AC"][uni_ac]:
                        gene = list(prot_dict["GN"][main_ac])[0]

                        #: determine the overlap with previously accepted LMs
                        overlap = calculate_overlap(elm_start, elm_end, masks[main_ac])

                        if overlap <= max_overlap:
                            lms[main_ac][elm_id].add("\t".join([ main_ac+"|"+gene, elm_id+"|"+elm_name, str(elm_start), str(elm_end), str(some_score)]))
                            lm_sets[elm_name].add(main_ac)
                                #: if the LM is accepted, its region in the
                                #: sequence mask is filled with "1"s
                            masks[main_ac] = fill_mask(elm_start, elm_end,
                                                        masks[main_ac])
    return lms, lm_sets

def reduce_3did(db_file, out_file):
    """ Writes a simplified version in TSV format of the 3did flat file db
    """
    with open_file(out_file, "w") as out:
        cols = ["Pfam_Name_A", "Pfam_Acc_A", "Pfam_Name_B", "Pfam_Acc_B", "PDBs"]
        out.write( "\t".join(cols)+"\n" )

        with open_file(db_file) as f:
            for line in f:
    #=ID    1-cysPrx_C      1-cysPrx_C       (PF10417.4@Pfam       PF10417.4@Pfam)
                if line.startswith("#=ID"):
                    pfam_name_a, pfam_name_b = line.rstrip().split()[1:3]
                    pfam_acc_a = line.rstrip().split()[3].split(".")[0].split("(")[1]
                    pfam_acc_b = line.rstrip().split()[4].split(".")[0]
                    pdbs = set()

    #=3D    1n8j    E:153-185       O:153-185       0.99    1.35657 0:0
                elif line.startswith("#=3D"):
                    pdb = line.rstrip().split()[1]
                    pdbs.add(pdb)

                elif line.startswith("//"):
                    row = [pfam_name_a, pfam_acc_a, pfam_name_b, pfam_acc_b,
                           ";".join(sorted(list(pdbs)))]
                    out.write( "\t".join(row)+"\n" )
    return

def edit_interprets(interprets_file, out_file):
    """ Placeholder function - does nothing right now
        Edit this function so it creates a custom interprets-results file that
        can be imported to Mongo
    """

    with open_file(interprets_file) as f:
        for line in f:
            if line[0]=="#":
                cols = line.rstrip().split("\t")
            else:
                tab = line.rstrip().split("\t")
                gn1, pdb1, eval1, pcid1, s1, e1 = tab[0:6]
                gn2, pdb2, eval2, pcid2, s2, e2 = tab[8:14]
                eval_avg = (float(eval1)+float(eval2))/2
                eval_diff = abs(float(eval1)-float(eval2))
                Z = "-"
                if len(tab) > 16:
                    Z = tab[20]
                info = [pdb1+":"+s1+"-"+e1+":"+eval1+":"+pcid1,
                        pdb2+":"+s2+"-"+e2+":"+eval2+":"+pcid2,
                        str(eval_avg), str(eval_diff), Z]

    return

def main( data_dir="",
          species="Hsa",
          max_overlap=0.2,
          no_overlap_between_doms_and_lms="y"
        ):

    # Common files:
    com_dir = data_dir+"common/"
    db_3did_file = com_dir+"3did_flat-2018_04.gz"

    # Species Files
    sp_data_dir   = data_dir+"species/"+species+"/"
    # uniprot_file  = sp_data_dir+"uniprot_homo_sapiens_proteome_73112prts_Aug2018_data.txt.gz"
    uniprot_file  = sp_data_dir+"uniprot_homo_sapiens_proteome_73928prts_Mar2019_data.txt.gz"
    pfam_matches_file = sp_data_dir+"pfamA_matches_9606_Aug18.tsv.gz"
    # pfam_hmm_file = sp_data_dir+"uniprot_v_Pfam_hmmsearch_sum.txt.gz"
    iprets_file   = sp_data_dir+"human_aaa_biogrid_i2.txt.gz"
    elm_hits_file = sp_data_dir+"elm_hits.tsv.gz"
    hippie_file   = sp_data_dir+"hippie_v2-1_July2018.tsv.gz"
    hippie_uniprot_mapping = sp_data_dir+"hippie_v2-1_uniprot_mapping_table.tsv.gz"

    # Output Files
    # pfam_parsed_file = sp_data_dir+"pfam_parsed_info.tsv.gz"
    pfam_assoc_file  = sp_data_dir+"dom_dom_association.tsv.gz"
    elm_parsed_file  = sp_data_dir+"elm_parsed_info_Overlapping.tsv.gz"
    edited_3did_file = com_dir+"3did_flat_edited-2018_04.tsv.gz"

    #1: Get protein ID-dictionary, sequences and masks
    prot_id_dict, masks, ref_prot = get_protein_data_from_uniprot_text(uniprot_file)

    # #2: Parse Pfam domains per protein from hmmsearch result file
           ### Deprecated: using pfamA_matches from now on. ###
    # if not os.path.isfile(pfam_parsed_file):
    #     domains, domain_sets, ms = parse_pfam_doms(pfam_hmm_file, prot_id_dict,
    #                                        masks, max_overlap=max_overlap)
    #
    #     #: Print them in TSV file
    #     with gzip.open(pfam_parsed_file, "wb") as out:
    #         c = ["#UniProt_acc|Gene", "Pfam_ID|Name", "Start", "End", "E-val"]
    #         out.write("\t".join(c) + "\n")
    #         for uni_ac in sorted(domains):
    #             for e_val in sorted(domains[uni_ac]):
    #                 for dom in domains[uni_ac][e_val]:
    #                     out.write(dom + "\n")

    #3: Get Linear Motifs per protein from hmmsearch result file
    if not os.path.isfile(elm_parsed_file):
        lms, lm_sets = parse_linear_motifs(elm_hits_file, prot_id_dict, masks,
                                        max_overlap=2)
        #: Print them in TSV file
        with gzip.open(elm_parsed_file, "wb") as out:
            c = ["#UniProt_acc|Gene", "ELM_ID|Name", "Start", "End", "Some_Score"]
            out.write("\t".join(c) + "\n")
            for uni_ac in sorted(lms):
                for lm in sorted(lms[uni_ac]):
                    for instance in lms[uni_ac][lm]:
                        out.write(instance + "\n")
    #4: Editing Files
    if not os.path.isfile(edited_3did_file):
        reduce_3did(db_3did_file, edited_3did_file)

    ## Calculating Dom-Dom Interactions (Association Method)
    ## using Swissprot proteins as reference set
    ## Parameters:
    ## (from HIPPIE) Medium confidence = 0.63 / High confidence = 0.73
    hippie_score = 0.63
    ## Minimum number of protein pairs with the domain-domain signature
    min_npair = 1
    ## domain count could be another parameter

    pfams, pfam_sets, pfam_names = get_pfam_doms(pfam_matches_file,
                                                prot_id_dict["AC"])

    # Source of PPI: HIPPIE (until I find something better)
    # not all ids from the hippie db matched a uniprot_ac, thus Uniprot mapping
    # tool was used and the results need to be read before:
    hippie_map = {}
    with open_file(hippie_uniprot_mapping) as f:
        for line in f:
            t = line.rstrip().split("\t")
            for x in t[0].split(","):
                hippie_map[x] = t[1]

    ppi = defaultdict(dict)
    with open_file(hippie_file) as f:
        for line in f:
            t = line.rstrip().split("\t")
            protsA = t[0].split(",")
            protsB = t[2].split(",")
            score = float(t[4])
            pmids = 0
            ## Score threshold to keep PPI
            ## (from HIPPIE) Medium confidence = 0.63 / High confidence = 0.73
            if score < hippie_score:
                continue
            info = t[5]
            if "pmids" in info:
                pmids = len(info.split(";")[1].split(":")[1].split(","))

            for protA in protsA:
                for protB in protsB:
                    if protA!="" and protB!="":
                        a = hippie_map[protA]
                        b = hippie_map[protB]
                        if a < b:
                            a = hippie_map[protB]
                            b = hippie_map[protA]
                        if (a in prot_id_dict["AC"] and b in prot_id_dict["AC"]):
                            a2 = prot_id_dict["AC"][a]
                            b2 = prot_id_dict["AC"][b]
                            for a in a2:
                                for b in b2:
                                    if (a in ref_prot and b in ref_prot):
                                        if (a in ppi and b in ppi[a]):
                                            if pmids > ppi[a][b][1]:
                                                ppi[a][b] = (score, pmids)
                                        else:
                                            ppi[a][b] = (score, pmids)

    # Count domain individual and pair frequencies in non redundant protein pairs
    # (domains are counted only once when they are repeated in the same protein)
    pfam_pair_count = defaultdict(lambda: defaultdict(int))
    total_prts = set()
    total_pp_pairs = 0
    for a in ppi:
        total_prts.add(a)
        for b in ppi[a]:
            total_prts.add(b)
            total_pp_pairs += 1

            if a in pfams and b in pfams:
                # Count domain pair frequency
                for pfam_a in pfams[a]:
                    for pfam_b in pfams[b]:
                        if pfam_a >= pfam_b:
                            pfam_pair_count[pfam_a][pfam_b] += 1
                        else:
                            pfam_pair_count[pfam_b][pfam_a] += 1

    # Count individual domain frequency
    pfam_count = defaultdict(int)
    # for prot in total_prts:
    #     for pfam in pfams[prot]:
    #         pfam_count[pfam] += 1

    ref_pfam = set()
    for uni_ac in ref_prot:
        if uni_ac in pfams:
            ref_pfam.add(uni_ac)
            for pfam in pfams[uni_ac]:
                pfam_count[pfam] += 1


    ## Compute domain-domain statistics
    ## (only those appearing in some PPI. To the rest we can assign the lowest
    ## LO score)
    lo = set()
    with gzip.open(pfam_assoc_file, "wb") as out:
        out.write("\t".join(["dom_ac_a","dom_ac_b","dom_name_a",
            "dom_name_b","dom_n_a","dom_n_b","obs","exp","or","lo"])+"\n")
        for pfam_a in sorted(pfam_pair_count):
            for pfam_b in sorted(pfam_pair_count[pfam_a]):
                n_pair = pfam_pair_count[pfam_a][pfam_b]
                n_a = pfam_count[pfam_a]
                n_b = pfam_count[pfam_b]
                exp = (n_a * n_b) / float(len(ref_prot) * len(ref_prot)) * total_pp_pairs
                oddsratio = float(n_pair) / exp
                log2 = math.log(oddsratio, 2)
                if n_pair >= min_npair:
                    out.write("\t".join([pfam_a, pfam_b,
                                    pfam_names[pfam_a], pfam_names[pfam_b],
                                    str(n_a), str(n_b), str(n_pair), str(exp),
                                    str(oddsratio), str(log2)])+"\n")

                    lo.add(log2)

    # print "#Reference proteome =",len(ref_prot),"with pfams =",len(ref_pfam)
    # print "#Total PPI pairs =", total_pp_pairs, "- involving", len(total_prts)
    # print "#Lowest / highest LO =",sorted(list(lo))[0],"/", sorted(list(lo))[-1]


if __name__ == "__main__":
    main()
