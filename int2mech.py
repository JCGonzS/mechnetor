#!/usr/bin/env python

""" INT2MECH.PY v2.0

Juan Carlos Gonzalez-Sanchez, Jan 2018
Python remake of Rob's script "int2mech.pl"
Description: ...

This script requires a series of pre-generated data files to run.
If these are missing, see "prepare_data_for_int2mech.py"

"""


import sys, re
import gzip, itertools, pymongo
from collections import defaultdict
from flask_debugtoolbar_lineprofilerpanel.profile import line_profile


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

    # ##From Mongodb
    # db = client['protein_data']
    # data = db['Hsa']


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
            {"_id": 0, "#BioGRID Interaction ID": 1}):

        bio_id = d["#BioGRID Interaction ID"]
        info.add(str(bio_id))

    return info

def get_3did_from_mongo(client, pfam_a, pfam_b):
    db = client['interactions_common']
    data = db['db3did']

    d = data.find_one( {"$or": [{"Pfam_Name_A": pfam_a, "Pfam_Name_B": pfam_b},
                 {"Pfam_Name_A": pfam_b, "Pfam_Name_B": pfam_a}]},
                 {"_id": 0, "PDBs": 1})
    if d:
        return d["PDBs"]
    else:
        return ""

def get_interprets_from_mongo(client, gene_a, gene_b):
    db = client['interactions_Hsa']
    data = db['iprets_Hsa']
    info = set()
    for d in data.find( {"$or": [{"#Gene1": gene_a, "Gene2": gene_b},
                 {"#Gene1": gene_b, "Gene2": gene_a} ]},
                 {"_id": 0, "PDB1":1, "qstart1": 1, "qend1": 1, "Blast-E1": 1, "Blast-PCID1": 1,
                 "PDB2":1, "qstart2": 1, "qend2": 1, "Blast-E2": 1, "Blast-PCID2": 1,
                 "Z": 1}):

        eval_avg = (float(d["Blast-E1"])+float(d["Blast-E2"]))/2
        eval_diff = abs(float(d["Blast-E1"])-float(d["Blast-E2"]))
        info.add( ";".join([str(d["PDB1"])+":"+str(d["qstart1"])+"-"+str(d["qend1"])+":"+str(d["Blast-E1"])+":"+str(d["Blast-PCID1"]),
                str(d["PDB2"])+":"+str(d["qstart2"])+"-"+str(d["qend2"])+":"+str(d["Blast-E2"])+":"+str(d["Blast-PCID2"]),
                str(eval_avg), str(eval_diff), str(d["Z"]) ]))

    return info

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

def pfams_from_acc(client, acc):
    db = client['protein_data']
    data = db['Hsa']
    cursor = data.find_one( { "uniprot_acc": acc },
                            { "_id": 0, "pfams.name": 1 } )
    pfams = [c["name"] for c in cursor["pfams"]]

    return sorted(list(set(pfams)))

def elms_from_acc(client, acc):
    db = client['protein_data']
    data = db['Hsa']
    cursor = data.find_one( { "uniprot_acc": acc },
                            { "_id": 0, "elms.name": 1 } )
    elms = [c["name"] for c in cursor["elms"]]

    return sorted(list(set(elms)))

@line_profile
def main(client, target_prots, protein_ids, protein_data, output_file="",
        data_dir="data/", species="Hsa", max_prots=""):

    # Parameters
    max_pmax = 1.0
    homo_int = "n"

    # Common files
    elm_int_dom_file = data_dir + "common/elm_interaction_domains_edited_Jan18.tsv"

    # Species files
    sp_data_dir = data_dir+"species/"+species+"/"
    dd_prop_file = sp_data_dir + "dom_dom_lo.txt"

    # pfam_sets = count_pfam_doms(protein_data)
    # elm_sets = count_elms(protein_data)

    # Get Interaction Data
    ## from interaction propensities [domain-domain]
        ## needs to be adjusted !!
    dom_prop_int = get_dd_prop_int(dd_prop_file,
                                    obs_min=4,
                                    lo_min=2.0,
                                    ndom_min=4)
    ## from ELM [domain-linear motif]
    elm_int = get_elm_ints(elm_int_dom_file)

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
            # score[n] = 1
            # n += 1

            ## can same pfam be repeated (i think so)
        pfams_a = pfams_from_acc(client, ac_a)
        pfams_b = pfams_from_acc(client, ac_b)
        for pfam_pair in itertools.product(pfams_a, pfams_b):
            pfam_a, pfam_b = pfam_pair

            # pa, pb, p, pmax = calculate_p(len(pfam_sets[pfam_a]),
            #                               len(pfam_sets[pfam_b]),
            #                               len(protein_data))
            # if pmax <= max_pmax:

            ## 3did interactions
            pdbs = get_3did_from_mongo(client, pfam_a, pfam_b)

            if pdbs != "":
                line = "\t".join([gene_a, ac_a, gene_b, ac_b,
                                  "DOM::DOM", pfam_a, pfam_b,
                                  pdbs, "3did"
                                 ])
                lines.append(line)
                # score[n] = p
                # n += 1


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
                # score[n] = p
                n += 1

        ## ELM-domain interactions
        elms_a = elms_from_acc(client, ac_a)
        elms_b = elms_from_acc(client, ac_b)
        for pfam_a in pfams_a:
            for elm_b in elms_b:
                if pfam_a in elm_int and elm_b in elm_int[pfam_a]:
                    only_prts = elm_int[pfam_a][elm_b]
                    if len(only_prts) > 0:
                        if gene_a not in only_prts:
                            continue
                    # pa, pb, p, pmax = calculate_p(len(pfam_sets[pfam_a]),
                    #                               len(elm_sets[elm_b]),
                    #                               len(protein_data))
                    line = "\t".join([gene_a, ac_a, gene_b, ac_b,
                                      "DOM::ELM", pfam_a, elm_b,
                                      "", "ELM DB"
                                     ])
                    lines.append(line)
                    # score[n] = p
                    n += 1

        for elm_a in elms_a:
            for pfam_b in pfams_b:
                if elm_a in elm_int and pfam_b in elm_int[elm_a]:
                    only_prts = elm_int[elm_a][pfam_b]
                    if len(only_prts) > 0:
                        if gene_b not in only_prts:
                            continue
                    # pa, pb, p, pmax = calculate_p(len(elm_sets[elm_a]),
                    #                               len(pfam_sets[pfam_b]),
                    #                               len(protein_data))
                    line = "\t".join([gene_a, ac_a, gene_b, ac_b,
                                      "ELM::DOM", elm_a, pfam_b,
                                      "", "ELM DB"
                                     ])
                    lines.append(line)
                    # score[n] = p
                    n += 1

        ## InterPreTS interactions
        hits = get_interprets_from_mongo(client, gene_a, gene_b)
        # if ac_a in iprets and ac_b in iprets[ac_a]:
            # for info in iprets[ac_a][ac_b]:
        for info in hits:
            info = info.split(";")
            line = "\t".join([gene_a, ac_a, gene_b, ac_b,
                              "InterPreTS", info[0], info[1],
                              "; ".join(info[2:]), "InterPreTS prediction"
                             ])
            lines.append(line)
            # score[n] = float(info[-3])
            n += 1


    cols = ["#Gene(A)","Accession(A)","Gene(B)","Accession(B)","Type",
            "F(A)","F(B)","Info", "Source"]
    if output_file:
        with open_file(output_file,"w") as out:
            out.write("\t".join(cols) + "\n")
            for line in lines:
                out.write(line+"\n")
            # for n in sorted(score, key=score.get, reverse=True):
            #     out.write(lines[n]+"\t"+str(score[n])+"\n")
