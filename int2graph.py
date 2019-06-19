#!/usr/bin/env python

import sys, re, os
import math, random, copy, pprint, json
import gzip, itertools, pymongo
from collections import defaultdict
from flask_debugtoolbar_lineprofilerpanel.profile import line_profile

def central_positions_layout(proteins):
    """Creates x/y coordinates of proteins to be displayed

    This function uses the parameterized archimedean spiral equation to
    generate the x/y coordinates of the central positions of the proteins in
    the network.
    """

    cpos = {}
    i = 0
    for uni_ac in sorted(list(proteins)):
        angle = i
        x = (1+100*angle)*math.cos(angle)
        y = (1+100*angle)*math.sin(angle)
        x += 650
        y += 100
        cpos[uni_ac] = (x, y)
        i = i+2

    return cpos

def muts_within_coords(ac, mutations, start, end):
    l = []
    for mut_pos in sorted(mutations):
        if mut_pos >= start and mut_pos <= end:
            for mut in mutations[mut_pos]:
                if mut not in l:
                    l.append(ac+"/"+mut)
    return l

def add_protein_main(nodes, id_n, id_dict, uni_ac, data, biogrid_data):

    id_n += 1
    prot_id = copy.deepcopy(id_n)
    id_dict[uni_ac] = defaultdict(list)
    id_dict[uni_ac]["main"] = prot_id
    biogrid_id = get_interactor_biogrid_id(biogrid_data, data["gene"])

    nodes.append({
        "group" : "nodes",
        "data" : {
            "id" : prot_id,
            "parent" : prot_id,
            "role" : "whole",
            "label" : data["gene"],
            "uni_id" : data["uni_id"],
            "biogrid_id" : biogrid_id,
            "des": data["description"],
            "length": data["length"],
            "protein" : uni_ac
            }
        })

    return nodes, id_n, id_dict, prot_id

def add_fasta_main(nodes, id_n, id_dict, header, seq, link, blast):

    id_n += 1
    prot_id = copy.deepcopy(id_n)
    id_dict[header] = defaultdict(list)
    id_dict[header]["main"] = prot_id

    nodes.append({
        "group" : "nodes",
        "data" : {
            "id" : prot_id,
            "parent" : prot_id,
            "role" : "user_seq",
            "label" : header,
            "link": link,
            "length": len(seq),
            "blast": blast
            }
        })

    return nodes, id_n, id_dict, prot_id

def add_start_end(nodes, edges, id_n, prot_id, uni_ac, len_seq, cp):

    start_x = cp[0] - len_seq/2
    start_y = cp[1]
    end_x = cp[0] + len_seq/2
    end_y = start_y

    ## node: protein start
    id_n += 1
    start_id = id_n
    nodes.append({
        "group" : "nodes",
        "data" : {
            "id" : start_id,
            "parent" : prot_id,
            "role" : "start-end",
            "label" : "0",
            "protein" : uni_ac
        },
        "position" : {
            "x" : start_x,
            "y" : start_y
        }
    })

    ## node: protein end
    id_n += 1
    end_id = id_n
    nodes.append({
        "group" : "nodes",
        "data" : {
            "id" : end_id,
            "parent" : prot_id,
            "role" : "start-end",
            "label" : str(len_seq),
            "protein" : uni_ac
        },
        "position" : {
            "x" : end_x,
            "y" : end_y
        }
    })

    ## edge: protein sequence
    id_n += 1
    edges.append({
        "group" : "edges",
        "data" : {
            "id" : id_n,
            "source" : start_id,
            "target" : end_id,
            "role" : "protein_sequence"
        }
    })

    return nodes, edges, id_n, (start_x, start_y)

def add_domains(nodes, id_n, id_dict, id_coords, id_muts, prot_id, uni_ac,
                pfams, pfam_info, start_xy, muts):

    """Add protein's Pfam domains

    Adds domain as a single node whose width is proportional to domain's length
    """

    pfam_set = set()
    start_x, start_y = start_xy
    for domain in pfams:
        start = domain["start"]
        end = domain["end"]
        length = end-start
        pfam_ac = domain["acc"].split(".")[0]
        pfam_set.add(pfam_ac)

        id_n += 1
        nodes.append({
            "group" : "nodes",
            "data" : {
                "id" : id_n,
                "role" : "domain",
                "label" : pfam_info[pfam_ac][0],
                "acc" : pfam_ac,
                "des" : pfam_info[pfam_ac][1],
                "start": str(start),
                "end": str(end),
                "length": str(length),
                "parent" : prot_id,
                "protein" : uni_ac
            },
            "position" : {
                "x" : start_x+start+(length/2),
                "y" : start_y
            }
        })
        id_dict[uni_ac][pfam_ac].append(id_n)
        id_coords[id_n] = (str(start), str(end))
        id_muts[id_n] = muts_within_coords(uni_ac, muts, start, end)

    return nodes, id_n, id_dict, id_coords, id_muts, pfam_set

def add_ptms(nodes, id_n, prot_id, uni_ac, data, start_xy):

    start_x, start_y = start_xy
    for ptm_type, x in zip(["phosphorylation", "acetylation"], ["p","ac"]):
        for ptm in data[ptm_type]:
            pos = ptm["pos"]
            res = ptm["res"]

            id_n += 1
            nodes.append({
                "group" : "nodes",
                "data" : {
                    "id" : id_n,
                    "parent" : prot_id,
                    "role" : ptm_type,
                    "label" : res+x+str(pos),
                    "protein" : uni_ac
                },
                "position" : {
                    "x" : start_x + int(pos),
                    "y": start_y
                }
            })

    return nodes, id_n

def add_mutagen(nodes, id_n, prot_id, uni_ac, mtgn, start_xy):

    start_x, start_y = start_xy
    for mtg in mtgn:
        pos = mtg["pos"]
        info = mtg["info"]
        id_n += 1
        nodes.append({
            "group" : "nodes",
            "data" : {
                "id" : id_n,
                "parent" : prot_id,
                "role" : "mutagen",
                "info" : info,
                "protein" : uni_ac
            },
            "position" : {
                "x" : start_x + int(pos),
                "y": start_y
            }
        })

    return nodes, id_n

def add_custom_mutations(nodes, id_n, prot_id, uni_ac, mutations, start_xy):
    """
    Adds user-input mutations as nodes within the proteins
    """

    start_x, start_y = start_xy
    for pos in mutations:
        label = ";".join(list(mutations[pos]))

        id_n += 1
        nodes.append({
            "group" : "nodes",
            "data" : {
                "id" : id_n,
                "parent" : prot_id,
                "role" : "input_mut",
                "label" : label,
                "protein" : uni_ac
            },
            "position" : {
                "x" : start_x + int(pos),
                "y" : start_y
            }
        })

    return nodes, id_n

def add_cosmic_mutations(nodes, id_n, prot_id, uni_ac, cosmic_data, start_xy):

    """
    Adds COSMIC mutations as nodes within the proteins
    Data extracted from internal MongoDB
    """

    muts = defaultdict(dict)
    for c in cosmic_data.find( {"uni_ac": uni_ac},
                               { "_id": 0, "enst": 0 }):
        pos = re.search("(\d+)",c["aa_mut"]).group(1)
        if int(c["samples"]) > 1:
            muts[pos][c["aa_mut"]] = (c["cosmic_id"], c["cds_mut"], str(c["samples"]))

    start_x, start_y = start_xy
    for pos in muts:
        cosmic_ids, aa_muts, cds_muts, count = [], [], [], []

        for aa_mut, (cosmic_id, cds_mut, samples) in muts[pos].iteritems():
            cosmic_ids.append("COSM"+str(cosmic_id))
            aa_muts.append(aa_mut)
            cds_muts.append(cds_mut)
            count.append(samples)

        id_n += 1
        nodes.append({
            "group": "nodes",
            "data": {
                "id": id_n,
                "parent": prot_id,
                "role": "cosmic_mut",
                "cos_id": "; ".join(cosmic_ids),
                "aa_mut": "; ".join(aa_muts),
                "cds": "; ".join(cds_muts),
                "count": "; ".join(count),
                "protein": uni_ac
            },
            "position": {
                "x": start_x + int(pos),
                "y": start_y
            }
        })

    return nodes, id_n

def get_Biogrid_from_MongoDB(data, gene_a, gene_b):
    evidence = defaultdict(set)
    for cursor in data.find( {"$or":
                                [{"Official Symbol Interactor A": gene_a,
                                 "Official Symbol Interactor B": gene_b},
                                 {"Official Symbol Interactor A": gene_b,
                                  "Official Symbol Interactor B": gene_a}]},
                                {"_id": 0, "#BioGRID Interaction ID": 1,
                                 "Throughput": 1}):

        thrput = cursor["Throughput"].split()[0]
        bio_id = cursor["#BioGRID Interaction ID"]
        evidence[thrput].add(str(bio_id))

    return evidence

def get_interactor_biogrid_id(data, gene):
    bio_id = "NA"
    cursor = data.find_one( {"$or":
                            [{"Official Symbol Interactor A": gene},
                             {"Official Symbol Interactor B": gene}]},
                             {"_id": 0,
                              "Official Symbol Interactor A": 1,
                              "Official Symbol Interactor B": 1,
                              "BioGRID ID Interactor A": 1,
                              "BioGRID ID Interactor B": 1 })
    if cursor:
        if cursor["Official Symbol Interactor A"] == gene:
            bio_id = cursor["BioGRID ID Interactor A"]
        else:
            bio_id = cursor["BioGRID ID Interactor B"]

    return bio_id

# def get_3did_from_MongoDB(data, pfam_a, pfam_b):
#
#     cursor = data.find_one( {"$or": [{"Pfam_Ide_A": pfam_a, "Pfam_Ide_B": pfam_b},
#                  {"Pfam_Ide_A": pfam_b, "Pfam_Ide_B": pfam_a}]},
#                  {"_id": 0, "PDBs": 1})
#     if cursor:
#         return cursor["PDBs"]

def get_3did_from_MongoDB(data, all_pfams):
    d = defaultdict(dict)
    cursor = data.find()
    for c in cursor:
        pfam_a = c["Pfam_Acc_A"]
        pfam_b = c["Pfam_Acc_B"]
        pdbs = c["PDBs"]
        if pfam_a in all_pfams and pfam_b in all_pfams:
            d[pfam_a][pfam_b]=str(pdbs)
    return d

def dom_dom_association_from_MongoDB(data, pfam_a, pfam_b,
                                     ndom_min, obs_min, lo_min):

    doc = data.find_one({"$or": [{"dom_ac_a": pfam_a, "dom_ac_b": pfam_b},
    					   {"dom_ac_a": pfam_b, "dom_ac_b": pfam_a}]},
    					    { "_id": 0, "dom_n_a":1, "dom_n_b":1,
                                "obs":1, "lo": 1 })
    if doc:
        if (int(doc["dom_n_a"])>=ndom_min and int(doc["dom_n_b"])>=ndom_min and
            int(doc["obs"])>=obs_min and float(doc["lo"])>=lo_min):
            return doc["lo"]

# def dom_dom_association_from_MongoDB(data, all_pfams, obs_min, lo_min, ndom_min):
#     d = defaultdict(dict)
#     cursor = data.find()
#     for c in cursor:
#         dom_a = c["dom_name_a"]
#         dom_b = c["dom_name_b"]
#         if (c["obs"]>=obs_min and c["lo"]>=lo_min
#         and c["dom_n_a"]>=ndom_min and c["dom_n_b"]>=ndom_min):
#             if (dom_a in all_pfams and dom_b in all_pfams):
#                 d[dom_a][dom_b] = c["lo"]
#     return d

def get_elm_dom_from_MongoDB(data):
    d = defaultdict(dict)
    for c in data.find():
        elm = c["#ELM identifier"]
        pfam = c["Domain identifier"]
        only_genes = c["Only in these genes"]
        d[elm][pfam] = only_genes
    return d

def get_Interprets_from_MongoDB(data, gene_a, ac_a, gene_b, ac_b):
    hits = set()
    # for d in data.find_( {"$or": [{"#Gene1": gene_a, "Gene2": gene_b},
    #              {"#Gene1": gene_b, "Gene2": gene_a} ]},
    #              {"_id": 0, "i2-raw": 0, "rand": 0, "rand-mean": 0, "rand-sd": 0,
    #               "p-value": 0,	"not-sure1": 0,	"not-sure2": 0}):

    d = data.find_one( {"$or": [{"#Gene1": gene_a, "Gene2": gene_b},
                                {"#Gene1": gene_b, "Gene2": gene_a},
                                {"#Gene1": ac_a, "Gene2": ac_b},
                                {"#Gene1": ac_b, "Gene2": ac_a} ]},
                 {"_id": 0, "i2-raw": 0, "rand": 0, "rand-mean": 0, "rand-sd": 0,
                  "not-sure1": 0,	"not-sure2": 0})
    if d:
        z = 0
        if "Z" in d:
            z = d["Z"]

        if d["#Gene1"] == gene_a:
            hit = (d["PDB1"], d["Blast-E1"], d["Blast-PCID1"],
                   d["qstart1"], d["qend1"], d["pdbstart1"], d["pdbend1"],
                   d["PDB2"], d["Blast-E2"], d["Blast-PCID2"],
                   d["qstart2"], d["qend2"], d["pdbstart2"], d["pdbend2"],
                   z, d["p-value"])
        else:
            hit = (d["PDB2"], d["Blast-E2"], d["Blast-PCID2"],
                   d["qstart2"], d["qend2"], d["pdbstart2"], d["pdbend2"],
                   d["PDB1"], d["Blast-E1"], d["Blast-PCID1"],
                   d["qstart1"], d["qend1"], d["pdbstart1"], d["pdbend1"],
                   z, d["p-value"])
        hits.add(hit)
    return hits

# def pfams_of_acc_from_MongoDB(protein_data, acc):
    #
    # cursor = protein_data.find_one( { "uni_ac": acc },
    #                         { "_id": 0, "pfams.acc": 1 } )
    # pfams = [c["acc"].split(".")[0] for c in cursor["pfams"]]
    #
    # return sorted(list(set(pfams)))

# def elms_of_acc_from_MongoDB(protein_data, acc):
#
#     cursor = protein_data.find_one( { "uni_ac": acc },
#                             { "_id": 0, "elms": 1 } )
#
#     # elm_names = sorted(list(set([c["name"] for c in cursor["elms"]])))
#     elms = defaultdict(list)
#     for c in cursor["elms"]:
#         elms[c["acc"]].append(c)
#
#     return elms

def get_random_color(already_used):
    """Get random color in hex-value
    """

    values = "123456789ABCDE"
    color = "#"
    for _ in range(100):
        for _ in range(6):
            color += values[int(round(random.random()*13))]
        if color not in already_used:
            break

    return color

def get_palette(palette):
    c = []
    if palette == "autumn":
        c = ["#CDBE70","#CD8500","#8B7500","#EE7600","#8B4513","#EE4000",
            "#528B8B","#4F94CD","#00688B", "#68838B", "#20B2AA","#8B7355"]
    elif palette == "custom1":
        c = ["#978897", "#E49AB0", "#494850", "#904C77", "#957D95"]
    return c

def color_regions(nodes, palette=""):
    """Gets colours for domains and ELMs graph elements

    Function arguments:
     nodes -- list of nodes
     palette -- if specified a custom palette of colours will be used. See
        "get_palette" for options. If empty, random colours will be used.
    """

    custom_colors = get_palette(palette)
    color_map = {}
    counter = 0
    for node in nodes:
        if node["group"] == "nodes":
            role = node["data"]["role"]
            if role in ["domain", "elm"]:
                label = node["data"]["label"]
                if label not in color_map:
                    if palette != "" and counter < len(custom_colors):
                        color_map[label] = custom_colors[counter]
                        counter += 1
                    else:
                        color_map[label] = get_random_color(color_map.values())

                node["data"]["color"] = color_map[label]

    return nodes

def color_from_zvalue(z_score):
    color = "grey"
    if z_score < 0 and z_score > -999999:
        color = "LightSkyBlue"
    elif z_score >= 0 and z_score < 1.65:
        color = "blue"
    elif z_score >= 1.65 and z_score < 2.33:
        color = "yellow"
    elif z_score >= 2.33 and z_score < 3.0:
        color = "orange"
    elif z_score >= 3.0:
        color = "red"

    return color

@line_profile
def main(make_graph, sps,
        target_prots, input_prots, custom_pairs, input_seqs, mutations,
        fasta_data, fasta_link, protein_data, cosmic_data, biogrid_data,
        iprets_data, fasta_iprets, db3did_data, dd_ass_data, elm_int_data,
        pfam_info, elm_info):

    cpos = central_positions_layout(target_prots | set(input_seqs.keys()))
    start_pos = {}
    nodes, edges = [], []
    id_n = 0
    id_dict = {}
    id_coords = {}
    id_muts = defaultdict(list)
    genes = {}
    pfams = {}
    elms = defaultdict(lambda: defaultdict(list))
    all_pfams = set()
    for uni_ac in sorted(list(target_prots)):

        data = protein_data.find_one( { "uni_ac": uni_ac },
                            { "_id": 0, "data_class": 0, "sequence": 0 })
        nodes, id_n, id_dict, prot_id = add_protein_main(nodes, id_n, id_dict,
                                                     uni_ac, data, biogrid_data)

        nodes, edges, id_n, start_pos[uni_ac] = add_start_end(nodes, edges, id_n,
                                                    prot_id, uni_ac,
                                                    data["length"], cpos[uni_ac])

        nodes, id_n, id_dict, id_coords, id_muts, pfams[uni_ac] = add_domains(nodes, id_n,
                                                    id_dict, id_coords, id_muts,
                                                    prot_id, uni_ac,
                                                    data["pfams"], pfam_info,
                                                    start_pos[uni_ac], mutations[uni_ac])

        nodes, id_n = add_ptms(nodes, id_n, prot_id, uni_ac, data,
                               start_pos[uni_ac])

        # nodes, id_n = add_mutagen(nodes, id_n, prot_id, uni_ac,
        #                           data["mutagen"], start_pos[uni_ac])

        nodes, id_n = add_custom_mutations(nodes, id_n, prot_id, uni_ac,
                                           mutations[uni_ac], start_pos[uni_ac])

        if cosmic_data != "no":
            nodes, id_n = add_cosmic_mutations(nodes, id_n, prot_id, uni_ac,
                                                cosmic_data, start_pos[uni_ac])

        genes[uni_ac] = data["gene"]
        all_pfams.update(pfams[uni_ac])
        for elm in data["elms"]:
            elms[uni_ac][elm["acc"]].append(elm)


    for header, seq in input_seqs.iteritems():

        data = fasta_data[header]
        link = fasta_link[header]
        nodes, id_n, id_dict, prot_id = add_fasta_main(nodes, id_n, id_dict,
                                                       header, seq, link, data["blast"])

        nodes, edges, id_n, start_pos[header] = add_start_end(nodes, edges, id_n,
                                                              prot_id, header,
                                                              len(seq), cpos[header])

        nodes, id_n, id_dict, id_coords, id_muts, pfams[header] = add_domains(nodes, id_n,
                                                    id_dict, id_coords, id_muts,
                                                    prot_id, header,
                                                    data["pfams"], pfam_info,
                                                    start_pos[header], mutations[header])

        nodes, id_n = add_custom_mutations(nodes, id_n, prot_id, header,
                                           mutations[header], start_pos[header])

        genes[header]=header
        all_pfams.update(pfams[header])
        elms[header]=data["elms"]

    dd_3did = get_3did_from_MongoDB(db3did_data, all_pfams)
    # dd_ass = dom_dom_association_from_MongoDB(dd_ass_data, all_pfams,
    #                                         obs_min=4, lo_min=2.0, ndom_min=4)
    elm_dom = get_elm_dom_from_MongoDB(elm_int_data)

    ass_dict = defaultdict(dict)
    elm_nodes = defaultdict(list)
    lines = []
    target_prots = target_prots | set(input_seqs.keys())
    for pair in itertools.combinations(target_prots, 2):

        ac_a, ac_b = pair
        # if make_graph==False: ## this is wrong. This is the option when you only want to see interactions between your input proteins (only to print)
        #     if ac_a not in input_prots and ac_b not in input_prots:
        #         continue
        gene_a = genes[ac_a]
        gene_b = genes[ac_b]

        if [ac_a, ac_b] in custom_pairs or [ac_b, ac_a] in custom_pairs:
            id_n += 1
            edges.append(   {"group" : "edges",
                             "data" : {  "id" : id_n,
                                         "source" : id_dict[ac_a]["main"],
                                         "target" : id_dict[ac_b]["main"],
                                         "role" : "user_interaction",
                                         "ds": "User input"
                                      }
                            })

        ## BioGRID interactions
        evidence = get_Biogrid_from_MongoDB(biogrid_data, gene_a, gene_b)
        biogrid_ints = list(evidence["Low"])+list(evidence["High"])

        if len(biogrid_ints) > 0:
            id_n += 1
            edges.append(   {"group" : "edges",
                             "data" : {  "id" : id_n,
                                         "source" : id_dict[ac_a]["main"],
                                         "target" : id_dict[ac_b]["main"],
                                         "role" : "prot_prot_interaction",
                                         "ds": "BioGRID",
                                         "low": list(evidence["Low"]),
                                         "high": list(evidence["High"])
                                      }
                            })

            line = [gene_a, ac_a, gene_b, ac_b, "PROT::PROT", "", "", "",
                    "", "", "", "; ".join(list(biogrid_ints)), "BioGRID"]
            lines.append(line)


        for pfam_pair in itertools.product(pfams[ac_a], pfams[ac_b]):
            pfam_a, pfam_b = pfam_pair
            ## 3did interactions
            # pdbs = get_3did_from_MongoDB(db3did_data, pfam_a, pfam_b)
            pdbs = ""
            if pfam_a in dd_3did and pfam_b in dd_3did[pfam_a]:
                pdbs = dd_3did[pfam_a][pfam_b]
            elif pfam_b in dd_3did and pfam_a in dd_3did[pfam_b]:
                pdbs = dd_3did[pfam_b][pfam_a]

            if pdbs:
                for i, source in enumerate(id_dict[ac_a][pfam_a], 1):
                    for j, target in enumerate(id_dict[ac_b][pfam_b], 1):
                        if source == target:
                            continue
                        id_n += 1
                        edges.append({ "group" : "edges",
                                       "data" : { "id" : id_n,
                                                  "source" : source,
                                                  "target" : target,
                                                  "role" : "DOM_interaction",
                                                  "ds": "3did",
                                                  "links": pdbs
                                                 }
                                    })

                        domA = pfam_a+":"+pfam_info[pfam_a][0]
                        domB = pfam_b+":"+pfam_info[pfam_b][0]
                        if len(id_dict[ac_a][pfam_a]) > 1:
                            domA += " ("+str(i)+")"
                        if len(id_dict[ac_b][pfam_b]) > 1:
                            domB += " ("+str(j)+")"
                        line = [gene_a, ac_a, gene_b, ac_b, "DOM::DOM",
                                domA, "-".join(id_coords[source]),
                                "; ".join(id_muts[source]),
                				domB, "-".join(id_coords[target]),
                                "; ".join(id_muts[target]),
                                pdbs, "3did"]
                        lines.append(line)

            ## Inferred Domain-Domain interactions
            if dd_ass_data != "no":
                if pfam_a in ass_dict and pfam_b in ass_dict[pfam_a]:
                    ass_lo = ass_dict[pfam_a][pfam_b]
                else:
                    ass_lo = dom_dom_association_from_MongoDB(dd_ass_data, pfam_a, pfam_b,
                                                     ndom_min=4, obs_min=5, lo_min=2.0)
                    ass_dict[pfam_a][pfam_b] = ass_lo
                    ass_dict[pfam_b][pfam_a] = ass_lo

                if ass_lo:
                    for i, source in enumerate(id_dict[ac_a][pfam_a], 1):
                        for j, target in enumerate(id_dict[ac_b][pfam_b], 1):
                            if source == target:
                                continue
                            id_n += 1
                            edges.append({ "group" : "edges",
                                           "data" : { "id" : id_n,
                                                      "source" : source,
                                                      "target" : target,
                                                      "lo": "{:.3f}".format(ass_lo),
                                                      "role" : "iDOM_interaction",
                                                      "ds": "Association Method"
                                                     }
                                        })

                            domA = pfam_a+":"+pfam_info[pfam_a][0]
                            domB = pfam_b+":"+pfam_info[pfam_b][0]
                            if len(id_dict[ac_a][pfam_a]) > 1:
                                domA += " ("+str(i)+")"
                            if len(id_dict[ac_b][pfam_b]) > 1:
                                domB += " ("+str(j)+")"
                            line = [gene_a, ac_a, gene_b, ac_b, "iDOM::iDOM",
                                    domA, "-".join(id_coords[source]),
                                    "; ".join(id_muts[source]),
                                    domB, "-".join(id_coords[target]),
                                    "; ".join(id_muts[target]),
                                    str(ass_lo), "Association Method"]
                            lines.append(line)

        ## ELM-domain interactions
        for ac1, gene1, ac2, gene2 in zip([ac_a, ac_b], [gene_a, gene_b],
                                          [ac_b, ac_a], [gene_b, gene_a]):

            for elm_acc in elms[ac1]:
                elm_ide = elm_info[elm_acc]["ide"]
                for pfam_acc in pfams[ac2]:
                    # only_prts = get_elm_dom_from_MongoDB(elm_int_data, elm_name, pfam)
                    # if only_prts:
                    #     if only_prts[0].strip() and gene2 not in only_prts:
                    #         continue
                    if elm_ide in elm_dom and pfam_acc in elm_dom[elm_ide]:
                        only_prts = elm_dom[elm_ide][pfam_acc]
                        if only_prts:
                            only_prts = only_prts.split(",")
                            if (sps=="Hsa"
                            and only_prts[0].strip() and gene2 not in only_prts):
                                    continue

                        ## Add ELM nodes if they don't exist
                        if elm_ide not in id_dict[ac1]:

                            elm_name  = elm_info[elm_acc]["name"]
                            elm_des   = elm_info[elm_acc]["des"]
                            elm_regex = elm_info[elm_acc]["regex"]
                            elm_prob  = elm_info[elm_acc]["prob"]

                            for elm in elms[ac1][elm_acc]:
                                start   = elm["start"]
                                end     = elm["end"]
                                length  = end-start
                                elm_seq = elm["seq"]

                                id_n += 1
                                nodes.append(
                                    { "group": "nodes",
                                      "data":
                                          { "id": id_n,
                                            "parent": id_dict[ac1]["main"],
                                            "role": "elm",
                                            "label": elm_ide,
                                            "acc": elm_acc,
                                            "name": elm_name,
                                            "des": elm_des,
                                            "regex": elm_regex,
                                            "start": str(start),
                                            "end": str(end),
                                            "seq": elm_seq,
                                            "length": str(length),
                                            "protein": ac1
                                          },
                                            "position": {
                                                "x": start_pos[ac1][0]+start+(length/2),
                                                "y": start_pos[ac1][1]
                                          }
                                })
                                id_dict[ac1][elm_ide].append(id_n)
                                id_coords[id_n] = (str(start), str(end))
                                id_muts[id_n] = muts_within_coords(ac1, mutations[ac1], start, end)

                        ## Add interaction edge
                        for i, source in enumerate(id_dict[ac1][elm_ide], 1):
                            for j, target in enumerate(id_dict[ac2][pfam_acc], 1):
                                id_n += 1
                                edges.append(
                                        { "group" : "edges",
                                          "data" :
                                            { "id" : id_n,
                                              "source" : source,
                                              "target" : target,
                                              "role" : "ELM_interaction",
                                              "ds": "ELM"
                                            }
                                        })

                                elmA = elm_acc+":"+elm_ide
                                domB = pfam_acc+":"+pfam_info[pfam_acc][0]
                                if len(id_dict[ac1][elm_ide]) > 1:
                                    elmA += " ("+str(i)+")"
                                if len(id_dict[ac2][pfam_acc]) > 1:
                                    domB += " ("+str(j)+")"
                                line = [gene1, ac1, gene2, ac2, "ELM::DOM",
                                         elmA, "-".join(id_coords[source]),
                                         "; ".join(id_muts[source]),
                                         domB, "-".join(id_coords[target]),
                                         "; ".join(id_muts[target]),
                                         "", "ELM"]
                                lines.append(line)

        if iprets_data == "no":
            continue

        ## InterPreTS interactions
        hits = []
        if ac_a in input_seqs or ac_b in input_seqs:
            if (ac_a, ac_b) in fasta_iprets:
                hits.append( fasta_iprets[(ac_a, ac_b)] )
            elif (ac_b, ac_a) in fasta_iprets:
                hits.append( fasta_iprets[(ac_b, ac_a)] )
        else:
            hits = get_Interprets_from_MongoDB(iprets_data, gene_a, ac_a, gene_b, ac_b)

        for hit in hits:
            pdb_a, eval_a, pcid_a, start_a, end_a, pdb_start_a, pdb_end_a = hit[:7]
            pdb_b, eval_b, pcid_b, start_b, end_b, pdb_start_b, pdb_end_b = hit[7:-2]
            z, pvalue = hit[-2:]

            label_a = pdb_a+":"+str(start_a)+"-"+str(end_a)
            label_b = pdb_b+":"+str(start_b)+"-"+str(end_b)
            eval_avg = (float(eval_a)+float(eval_b)) / 2
            eval_diff = abs(float(eval_a)-float(eval_b))
            color = color_from_zvalue(z)

            ## Add homology region node
            for ac, label, pdb, pdb_start, pdb_end, start, end, e_val, pcid in zip(pair,
                                    [label_a, label_b], [pdb_a, pdb_b],
                                    [pdb_start_a, pdb_start_b], [pdb_end_a, pdb_end_b],
                                    [start_a, start_b], [end_a, end_b],
                                    [eval_a, eval_b], [pcid_a, pcid_b]):

                length = int(end)-int(start)
                id_n += 1
                nodes.append(
                    { "group" : "nodes",
                      "data" :
                          { "id": id_n,
                            "parent": id_dict[ac]["main"],
                            "role": "iprets",
                            "label": label,
                            "pdb": pdb,
                            "pdb_start": int(pdb_start),
                            "pdb_end": int(pdb_end),
                            "start": int(start),
                            "end": int(end),
                            "length": length,
                            "eval": "{:1.0e}".format(float(e_val)),
                            "pcid": pcid,
                            "color": color,
                            "protein": ac
                          },
                            "position" : {
                                "x" : start_pos[ac][0]+int(start)+(length/2),
                                "y" : start_pos[ac][1]
                          }
                })
                id_dict[ac][label] = id_n
                id_muts[id_n] = muts_within_coords(ac, mutations[ac],
                                                   int(start), int(end))

            ## Add interaction edge
            source = id_dict[ac_a][label_a]
            target = id_dict[ac_b][label_b]
            id_n += 1
            edges.append(
                { "group" : "edges",
                  "data" :
                    { "id" : id_n,
                      "source" : source,
                      "target" : target,
                      "role" : "INT_interaction",
                      "pdb": pdb,
                      "color": color,
                      "z-score": z,
                      "p-value": "{:1.0e}".format(float(pvalue))
                    }
                })

            line = [gene_a, ac_a, gene_b, ac_b, "InterPreTS",
                   label_a, str(start_a)+"-"+str(end_a), "; ".join(id_muts[source]),
                   label_b, str(start_b)+"-"+str(end_b), "; ".join(id_muts[target]),
                   str(z), "InterPreTS prediction"]

            lines.append(line)

    ## Color domain & LMs nodes
    nodes = color_regions(nodes)
    # nodes = color_regions(nodes, palette="custom1")

    graph_elements = nodes + edges

    return  graph_elements, lines

    # ## Print graph as JSON file
    # graph_elements = nodes + edges
    # with open(graph_out, "w") as output:
    #     json.dump(graph_elements, output, indent=4)
    #
    # ## Print interactions as JSON file
    # for l in lines:
    #     print "\t".join(l)
    # int_table = {}
    # int_table["columns"] = ["#Gene(A)","Accession(A)","Gene(B)","Accession(B)",
    #                         "Type", "F(A)","Start-End(A)", "Mutations(A)",
    #                         "F(B)", "Start-End(B)", "Mutations(B)",
    #                         "Info", "Source"]
    # int_table["index"] = range(len(lines))
    # int_table["data"] = lines
    #
    # with open(ints_out,"w") as output:
    #     json.dump(int_table, output)
    #     # output.write(str(int_table))
