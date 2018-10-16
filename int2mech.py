#!/usr/bin/env python

""" INT2MECH.PY v2.0

Juan Carlos Gonzalez-Sanchez, Jan 2018
Python remake of Rob's script "int2mech.pl"
Description: ...

This script requires a series of pre-generated data files to run.
If these are missing, see "prepare_data_for_int2mech.py"

"""


import sys, re
import gzip, itertools
from collections import defaultdict
from flask_debugtoolbar_lineprofilerpanel.profile import line_profile
from pymongo import MongoClient

def open_file(input_file, mode="r"):
    """ Opens file Zipped or not
    """

    if re.search(".gz$", input_file):
        infile = gzip.open(input_file, mode)
    else:
        infile = open(input_file, mode)
    return infile

def count_pfam_doms(protein_data, e_val=0.1):
    """ Gets the set of unique proteins that contain every Pfam domain

    """
    ## Requires calibration... how many proteins am I counting?
    ## should be only those from proteome? or from swissprot?

    pfam_sets = defaultdict(list)
    for uni_ac in protein_data:
        for pfam_name in protein_data[uni_ac]["pfams"]:
            for domain in protein_data[uni_ac]["pfams"][pfam_name]:
                if domain["e-val"] < e_val:
                    pfam_sets[pfam_name].append(uni_ac)

    return pfam_sets

def count_elms(protein_data, e_val=0.1):
    """ Gets the set of unique proteins that contain every ELM

    """
    ## Requires calibration... how many proteins am I counting?
    ## should be only those from proteome? or from swissprot?

    elm_sets = defaultdict(list)
    for uni_ac in protein_data:
        for elm_name in protein_data[uni_ac]["elms"]:
            for elm in protein_data[uni_ac]["elms"][elm_name]:
                if elm["score"] < e_val:
                    elm_sets[elm_name].append(uni_ac)

    return elm_sets

def get_biogrid_from_mongo(client, gene_a, gene_b):
            # modify for single protein_input
    db = client['interactions_Hsa']
    data = db['biogrid_Hsa']
    info = set()
    for d in data.find( {"$or":
            [{"Official Symbol Interactor A": gene_a,
             "Official Symbol Interactor B": gene_b},
             {"Official Symbol Interactor A": gene_b,
              "Official Symbol Interactor B": gene_a}]},
            {"#BioGRID Interaction ID": 1}):

        bio_id = d["#BioGRID Interaction ID"]
        info.add(str(bio_id))

    return info

def get_3did_from_mongo(client, pfam_a, pfam_b):
    db = client['interactions_common']
    data = db['db3did']

    d = data.find_one( {"$or": [{"Pfam_Name_A": pfam_a, "Pfam_Name_B": pfam_b},
                 {"Pfam_Name_A": pfam_b, "Pfam_Name_B": pfam_a}]}, {"PDBs": 1})
    if d:
        return d["PDBs"]
    else:
        return ""

def get_dd_prop_int(dd_prop_file, obs_min, lo_min, ndom_min):
    """Gets domain-domain interactions from the pre-computed propensities
    """

    dom_int = defaultdict(dict)
    with open_file(dd_prop_file) as f:
        for line in f:
            if line[0]!="#":
                t = line.rstrip().split("\t")
                pfam_a, pfam_b = t[0:2]
                obs = int(t[2])
                lo = float(t[4])
                n_a, n_b = int(t[5]), int(t[6])
                if (obs >= obs_min and lo >= lo_min and
                    n_a >= ndom_min and n_b >= ndom_min):
                    dom_int[pfam_a][pfam_b] = "{} {} {} {:.3f}".format(obs,
                                                                 n_a, n_b, lo)
                    dom_int[pfam_b][pfam_a] = "{} {} {} {:.3f}".format(obs,
                                                                 n_b, n_a, lo)
    return dom_int

def get_elm_ints(elm_int_dom_file):
    """Gets ELM-domain interactions from the ELM database file
    """

    interactions = defaultdict(dict)
    with open_file(elm_int_dom_file) as f:
        for line in f:
            t = line.rstrip().split("\t")
            elm_id, pfam_ac, pfam_name, pfam_desc = t[:4]
            only_prts = []
            if len(t)>4:
                only_prts = t[4].split(",")
            interactions[elm_id][pfam_name] = only_prts
            interactions[pfam_name][elm_id] = only_prts

    return interactions

def get_interprets_from_mongo(client, gene_a, gene_b):
    db = client['interactions_Hsa']
    data = db['iprets_Hsa']
    


def get_interprets(interprets_file, target_prots, ac_dict):
    """Get InterPreTS predicted region-region interactions from a previously
    precomputed file from interacting proteins in BioGRID
    """

    iprets = defaultdict(lambda: defaultdict(set))
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

                if gn1.upper() in ac_dict and gn2.upper() in ac_dict:
                    ac1 = ac_dict[gn1.upper()]
                    ac2 = ac_dict[gn2.upper()]

                    if ac1 in target_prots and ac2 in target_prots:
                        iprets[ac1][ac2].add(";".join(info))
                        iprets[ac2][ac1].add(";".join(info))
                        # Like this, every interaction is taken. If there is a
                        # lot of overlapping regions at the end, I need to add
                        # a methods that searches that there aren't similar
                        # interactions already taking (compare starts & ends)
    return iprets

def calculate_p(na, nb, N):
    """Given 2 domains/elms calculates their frequencies and the probability
    of them being found in 2 random interacting proteins
    """

    pa = na / float(N)
    pb = nb / float(N)
    p = pa * pb
    pmax = pa
    if pb > pa:
        pmax = pb

    return pa, pb, p, pmax

@line_profile
def main(target_prots, protein_ids, protein_data, output_file="",
        data_dir="data/", species="Hsa", max_prots=""):

    # Parameters
    max_pmax = 1.0
    homo_int = "n"

    client = MongoClient('localhost', 27017)


    # Common files
    elm_int_dom_file = data_dir + "common/elm_interaction_domains_edited_Jan18.tsv"

    # Species files
    sp_data_dir = data_dir+"species/"+species+"/"
    dd_prop_file = sp_data_dir + "dom_dom_lo.txt"
    interprets_file = sp_data_dir + "human_aaa_biogrid_i2.txt.gz"

    pfam_sets = count_pfam_doms(protein_data)
    elm_sets = count_elms(protein_data)

    # Get Interaction Data

    ## from interaction propensities [domain-domain]
        ## needs to be adjusted !!
    dom_prop_int = get_dd_prop_int(dd_prop_file,
                                    obs_min=4,
                                    lo_min=2.0,
                                    ndom_min=4)
    ## from ELM [domain-linear motif]
    elm_int = get_elm_ints(elm_int_dom_file)
    ### Get InterpreTS interactions
    iprets = get_interprets(interprets_file, target_prots, protein_ids["AC"])
    #FIX: add here a function to remove redundancy and overlapping
        #within InterpreTS interactions.

    ## Print all interactions
    lines = []
    score = {}
    n = 0
    pairs = []

    for pair in itertools.combinations(target_prots, 2):
        ac_a, ac_b = pair
        # if homo_int == "n":
        #     if ac_a == ac_b:
        #         continue
    # for ac_1 in sorted(list(target_prots)):
    #     for ac_2 in sorted(list(target_prots)):
    #         if ac_1 < ac_2:
    #             ac_a = ac_1
    #             ac_b = ac_2
    #         else:
    #             ac_a = ac_2
    #             ac_b = ac_1
            #
            # if ac_a+"-"+ac_b in pairs:
            #     continue
            # pairs.append(ac_a+"-"+ac_b)

        gene_a = protein_ids["GN"][ac_a]
        gene_b = protein_ids["GN"][ac_b]

        ## BioGRID interactions
        info = get_biogrid_from_mongo(client, gene_a, gene_b)

        if len(info)>0:
            line = "\t".join([gene_a, ac_a, gene_b, ac_b,
                              "PROT::PROT", "", "",
                              "; ".join(list(info)), "BioGRID"
                             ])
            lines.append(line)
            score[n] = 1
            n += 1


        for pfam_pair in itertools.product(protein_data[ac_a]["pfams"], protein_data[ac_b]["pfams"]):
            pfam_a, pfam_b = pfam_pair

            pa, pb, p, pmax = calculate_p(len(pfam_sets[pfam_a]),
                                          len(pfam_sets[pfam_b]),
                                          len(protein_data))
            # if pmax <= max_pmax:

            ## 3did interactions
            pdbs = get_3did_from_mongo(client, pfam_a, pfam_b)

            if pdbs != "":
                line = "\t".join([gene_a, ac_a, gene_b, ac_b,
                                  "DOM::DOM", pfam_a, pfam_b,
                                  pdbs, "3did"
                                 ])
                lines.append(line)
                score[n] = p
                n += 1

            # if (pfam_a in dom_3did_int and
            #     pfam_b in dom_3did_int[pfam_a]):
            #     info = dom_3did_int[pfam_a][pfam_b]
            #     line = "\t".join([gene_a, ac_a, gene_b, ac_b,
            #                       "DOM::DOM", pfam_a, pfam_b,
            #                       info[0], "3did"
            #                      ])
            #                         # "; ".join(list(info))])
            #     lines.append(line)
            #     score[n] = p
            #     n += 1

            ## domain propensities*
            if (pfam_a in dom_prop_int and
                pfam_b in dom_prop_int[pfam_a]):
                info = dom_prop_int[pfam_a][pfam_b]
                line = "\t".join([gene_a, ac_a, gene_b, ac_b,
                                  "iDOM::iDOM", pfam_a, pfam_b,
                                  "; ".join(list(info)),
                                  "Statistical Prediction"
                                 ])
                lines.append(line)
                score[n] = p
                n += 1

        ## ELM-domain interactions
        for pfam_a in protein_data[ac_a]["pfams"]:
            for elm_b in protein_data[ac_b]["elms"]:
                if pfam_a in elm_int and elm_b in elm_int[pfam_a]:
                    only_prts = elm_int[pfam_a][elm_b]
                    if len(only_prts) > 0:
                        if gene_a not in only_prts:
                            continue
                    pa, pb, p, pmax = calculate_p(len(pfam_sets[pfam_a]),
                                                  len(elm_sets[elm_b]),
                                                  len(protein_data))
                    line = "\t".join([gene_a, ac_a, gene_b, ac_b,
                                      "DOM::ELM", pfam_a, elm_b,
                                      "", "ELM DB"
                                     ])
                    lines.append(line)
                    score[n] = p
                    n += 1

        for elm_a in protein_data[ac_a]["elms"]:
            for pfam_b in protein_data[ac_b]["pfams"]:
                if elm_a in elm_int and pfam_b in elm_int[elm_a]:
                    only_prts = elm_int[elm_a][pfam_b]
                    if len(only_prts) > 0:
                        if gene_b not in only_prts:
                            continue
                    pa, pb, p, pmax = calculate_p(len(elm_sets[elm_a]),
                                                  len(pfam_sets[pfam_b]),
                                                  len(protein_data))
                    line = "\t".join([gene_a, ac_a, gene_b, ac_b,
                                      "ELM::DOM", elm_a, pfam_b,
                                      "", "ELM DB"
                                     ])
                    lines.append(line)
                    score[n] = p
                    n += 1

        ## InterPreTS interactions
        if ac_a in iprets and ac_b in iprets[ac_a]:
            for info in iprets[ac_a][ac_b]:
                info = info.split(";")
                line = "\t".join([gene_a, ac_a, gene_b, ac_b,
                                  "InterPreTS", info[0], info[1],
                                  "; ".join(info[2:]), "InterPreTS prediction"
                                 ])
                lines.append(line)
                score[n] = float(info[-3])
                n += 1


    cols = ["#Gene(A)","Accession(A)","Gene(B)","Accession(B)","Type",
            "F(A)","F(B)","Info", "Source"]
    if output_file:
        with open(output_file,"w") as out:
            #out.write("#File generated with \"int2mech.py\"\n")
            #out.write("# Interprets file:"+interprets_file+"\n")
            out.write("\t".join(cols) + "\n")
            for line in lines:
                out.write(line+"\n")
            # for n in sorted(score, key=score.get, reverse=True):
            #     out.write(lines[n]+"\t"+str(score[n])+"\n")

    # else:
    #     print "\t".join(cols)
    #     for line in lines:
    #         print line
    #     for n in sorted(score, key=score.get, reverse=True):
    #         print lines[n]+"\t"+str(score[n])
