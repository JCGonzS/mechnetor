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

    central_pos = {}
    i = 0
    for prot_acc in sorted(list(proteins)):
        angle = i
        x = (1+100*angle)*math.cos(angle)
        y = (1+100*angle)*math.sin(angle)
        x += 650
        y += 100
        central_pos[prot_acc] = (x, y)
        i = i+2

    return central_pos


def add_protein(protein_data, cosmic_data, biogrid_data, prot_acc, central_pos,
                mutations, nodes, id_counter, id_dict, start_pos):

    """Adds a protein to the graph elements dictionary in "cytoscape-format"

    This function first adds three nodes per protein to the graph:
    the first node is the parent node and represents the whole protein, the
    other two nodes are the children of the first node and mark the start and
    end of the protein.

    Function also calls others to add the domains, ELMs, mutations, PTMs of
    the protein, which are checked to see if they are located at the start or
    at the end.
    """
    id_dict[prot_acc] = defaultdict(list)
    protein_id = copy.deepcopy(id_counter)
    cursor = protein_data.find_one( { "uniprot_acc": prot_acc },
                        { "_id": 0, "cosmic_muts": 0, "data_class": 0 })
    gene = cursor["gene"]
    biogrid_id = get_interactor_biogrid_id(biogrid_data, gene)
    ## Add protein central node
    nodes.append({
        "group" : "nodes",
        "data" : {
            "id" : protein_id,
            "parent" : protein_id,
            "role" : "whole",
            "label" : cursor["gene"],
            "uni_id" : cursor["uniprot_id"],
            "biogrid_id" : biogrid_id,
            "des": cursor["description"],
            "length": cursor["length"],
            "protein" : prot_acc
            }
        })
    id_dict[prot_acc]["main"] = protein_id
    id_counter += 1

    ## Add protein start & end nodes (protein sequence will be drawn between the two)
    start_x = central_pos[0] - cursor["length"]/2
    start_y = central_pos[1]
    end_x = central_pos[0] + cursor["length"]/2
    end_y = start_y
    start_pos[prot_acc] = (start_x, start_y)

    ## protein start
    nodes.append({
        "group" : "nodes",
        "data" : {
            "id" : id_counter,
            "parent" : protein_id,
            "role" : "start-end",
            "label" : "0",
            "protein" : prot_acc
        },
        "position" : {
            "x" : start_x,
            "y" : start_y
        }
    })
    id_counter += 1

    ## protein end
    nodes.append({
        "group" : "nodes",
        "data" : {
            "id" : id_counter,
            "parent" : protein_id,
            "role" : "start-end",
            "label" : str(cursor["length"]),
            "protein" : prot_acc
        },
        "position" : {
            "x" : end_x,
            "y" : end_y
        }
    })
    id_counter += 1

    ## Add other protein elements:
    nodes, id_counter, id_dict = add_domains(prot_acc, protein_id, cursor,
                                    start_x, start_y, nodes, id_counter, id_dict)

    nodes, id_counter = add_ptms(prot_acc, protein_id, cursor,
                                 start_x, start_y, nodes, id_counter)

    nodes, id_counter = add_custom_mutations(prot_acc, protein_id,
                                mutations, start_x, start_y, nodes, id_counter)

    nodes, id_counter = add_cosmic_mutations(prot_acc, protein_id,cosmic_data,
                                            start_x, start_y, nodes, id_counter)


    return nodes, id_counter, id_dict, start_pos


def add_domains(prot_acc, parent_id, cursor, start_x, start_y,
                nodes, id_counter, id_dict):

    """Add protein's Pfam domains

    Adds domain as a single node whose width is proportional to domain's length
    """

    for domain in cursor["pfams"]:
        start = domain["start"]
        end = domain["end"]
        length = end-start

        nodes.append({
            "group" : "nodes",
            "data" : {
                "id" : id_counter,
                "parent" : parent_id,
                "role" : "domain",
                "label" : domain["name"],
                "acc" : domain["acc"].split(".")[0],
                "des" : domain["des"],
                "start": str(start),
                "end": str(end),
                "length": str(length),
                "protein" : prot_acc
            },
            "position" : {
                "x" : start_x+start+(length/2),
                "y" : start_y
            }
        })
        id_dict[prot_acc][domain["name"]].append(id_counter)
        id_counter += 1


    return nodes, id_counter, id_dict


def add_ptms(prot_acc, parent_id, cursor, start_x, start_y,
             nodes, id_counter):

    for ptm_type, x in zip(["phosphorylation", "acetylation"], ["p","ac"]):
        for ptm in cursor[ptm_type]:
            pos = ptm["pos"]
            res = ptm["res"]

            nodes.append({
                "group" : "nodes",
                "data" : {
                    "id" : id_counter,
                    "parent" : parent_id,
                    "role" : ptm_type,
                    "label" : res+x+str(pos),
                    "protein" : prot_acc
                },
                "position" : {
                    "x" : start_x + int(pos),
                    "y": start_y
                }
            })
            id_counter += 1

    return nodes, id_counter


def add_custom_mutations(prot_acc, parent_id, mutations,
                         start_x, start_y, nodes, id_counter):
    """
    Adds user-input mutations as nodes within the proteins
    """

    for pos in mutations:
        label = ";".join(list(mutations[pos]))

        nodes.append({
            "group" : "nodes",
            "data" : {
                "id" : id_counter,
                "parent" : parent_id,
                "role" : "input_mut",
                "label" : label,
                "protein" : prot_acc
            },
            "position" : {
                "x" : start_x + int(pos),
                "y" : start_y
            }
        })
        id_counter += 1

    return nodes, id_counter

def add_cosmic_mutations(prot_acc, parent_id, cosmic_data,
                         start_x, start_y, nodes, id_counter):

    """
    Adds COSMIC mutations as nodes within the proteins
    Data extracted from internal MongoDB
    """
    muts = defaultdict(dict)
    for c in cosmic_data.find( {"uni_ac": prot_acc},
                                    { "_id": 0, "enst": 0 }):
        pos = re.search("(\d+)",c["aa_mut"]).group(1)
        if int(c["samples"]) > 1:
            muts[pos][c["aa_mut"]] = (c["cosmic_id"], c["cds_mut"], str(c["samples"]))

    for pos in muts:
        cosmic_ids, aa_muts, cds_muts, count = [], [], [], []

        for aa_mut, (cosmic_id, cds_mut, samples) in muts[pos].iteritems():
            cosmic_ids.append("COSM"+str(cosmic_id))
            aa_muts.append(aa_mut)
            cds_muts.append(cds_mut)
            count.append(samples)

        nodes.append({
            "group": "nodes",
            "data": {
                "id": id_counter,
                "parent": parent_id,
                "role": "cosmic_mut",
                "cos_id": "; ".join(cosmic_ids),
                "aa_mut": "; ".join(aa_muts),
                "cds": "; ".join(cds_muts),
                "count": "; ".join(count),
                "protein": prot_acc
            },
            "position": {
                "x": start_x + int(pos),
                "y": start_y
            }
        })
        id_counter += 1

    return nodes, id_counter


def connect_protein_sequence(prot_acc, nodes, edges, id_counter):

    # 1. Decide which nodes form the protein sequence
    prot_positions = []
    for ele in nodes:
        role = ele["data"]["role"]
        parent_prot = ele["data"]["protein"]
        if (role == "start-end" and parent_prot == prot_acc):
            position = ele["data"]["label"]
            ele_id = ele["data"]["id"]
            prot_positions.append((ele_id, int(position)))

    # 2. Join those nodes by edges
    prot_positions.sort(key=lambda x: x[1])
    for i in range(len(prot_positions)-1):
        edges.append({
            "group" : "edges",
            "data" : {
                "id" : id_counter,
                "source" : prot_positions[i][0],
                "target" : prot_positions[i+1][0],
                "role" : "protein_sequence"
            }
        })
        id_counter += 1

    return edges, id_counter


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
    bio_id = ""
    cursor = data.find_one( {"$or":
                            [{"Official Symbol Interactor A": gene},
                             {"Official Symbol Interactor B": gene}]},
                             {"_id": 0,
                              "Official Symbol Interactor A": 1,
                              "Official Symbol Interactor B": 1,
                              "BioGRID ID Interactor A": 1,
                              "BioGRID ID Interactor B": 1 })

    if cursor["Official Symbol Interactor A"] == gene:
        bio_id = cursor["BioGRID ID Interactor A"]
    else:
        bio_id = cursor["BioGRID ID Interactor B"]

    return bio_id


def get_3did_from_MongoDB(data, pfam_a, pfam_b):

    cursor = data.find_one( {"$or": [{"Pfam_Name_A": pfam_a, "Pfam_Name_B": pfam_b},
                 {"Pfam_Name_A": pfam_b, "Pfam_Name_B": pfam_a}]},
                 {"_id": 0, "PDBs": 1})
    if cursor:
        return cursor["PDBs"]


def domain_propensities_from_MongoDB(data, pfam_a, pfam_b,
                                     obs_min, lo_min, ndom_min):

    cursor = data.find_one({"$or": [{"#DOM1": pfam_a, "DOM2": pfam_b},
    					   {"#DOM1": pfam_b, "DOM2": pfam_a}],
    				  	   "OBS": {"$gte": obs_min}, "LO": {"$gte": lo_min},
    				  	   "N_DOM1": {"$gte": ndom_min},
                           "N_DOM2": {"$gte": ndom_min}},
    					    { "_id": 0, "LO": 1 })
    if cursor:
        return cursor["LO"]


def get_elm_dom_from_MongoDB(data):

    d = defaultdict(dict)
    cursor = data.find()
    for c in cursor:
        elm = c["#ELM identifier"]
        pfam = c["Domain Name"]
        only_genes = c["Only in these genes"]
        d[elm][pfam] = only_genes
    return d

def get_elm_info(data):
## Collection fields:
#    Accession       ELMIdentifier   FunctionalSiteName      Description
#    Regex   Probability     #Instances      #Instances_in_PDB
    d = {}
    cursor = data.find()
    for c in cursor:
        d[c["ELMIdentifier"]] = (c["Accession"], c["FunctionalSiteName"],
                                 c["Regex"])

    return d

def get_Interprets_from_MongoDB(data, gene_a, gene_b):
    hits = set()
    for d in data.find( {"$or": [{"#Gene1": gene_a, "Gene2": gene_b},
                 {"#Gene1": gene_b, "Gene2": gene_a} ]},
                 {"_id": 0, "i2-raw": 0, "rand": 0, "rand-mean": 0, "rand-sd": 0,
                  "p-value": 0,	"not-sure1": 0,	"not-sure2": 0}):

        if d["#Gene1"] == gene_a:
            hit = (d["PDB1"], d["Blast-E1"], d["Blast-PCID1"],
                   d["qstart1"], d["qend1"], d["pdbstart1"], d["pdbend1"],
                   d["PDB2"], d["Blast-E2"], d["Blast-PCID2"],
                   d["qstart2"], d["qend2"], d["pdbstart2"], d["pdbend2"],
                   d["Z"])
        else:
            hit = (d["PDB2"], d["Blast-E2"], d["Blast-PCID2"],
                   d["qstart2"], d["qend2"], d["pdbstart2"], d["pdbend2"],
                   d["PDB1"], d["Blast-E1"], d["Blast-PCID1"],
                   d["qstart1"], d["qend1"], d["pdbstart1"], d["pdbend1"],
                   d["Z"])
        hits.add(hit)
    return hits


def pfams_of_acc_from_MongoDB(protein_data, acc):

    cursor = protein_data.find_one( { "uniprot_acc": acc },
                            { "_id": 0, "pfams.name": 1 } )
    pfams = [c["name"] for c in cursor["pfams"]]

    return sorted(list(set(pfams)))


def elms_of_acc_from_MongoDB(protein_data, acc):

    cursor = protein_data.find_one( { "uniprot_acc": acc },
                            { "_id": 0, "elms": 1 } )

    # elm_names = sorted(list(set([c["name"] for c in cursor["elms"]])))
    elms = defaultdict(list)
    for c in cursor["elms"]:
        elms[c["name"]].append(c)

    return elms

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
def main(target_prots, custom_pairs, protein_data, cosmic_data, mutations,
        biogrid_data, iprets_data, db3did_data, dom_prop_data, elm_int_data,
        elm_classes,
        max_prots, graph_out, ints_out):

    central_pos = central_positions_layout(target_prots)
    start_pos = {}
    nodes, edges = [], []
    id_counter = 0
    id_dict = {}
    pfams, elms = {}, {}
    for prot_acc in sorted(list(target_prots)):
        nodes, id_counter, id_dict, start_pos = add_protein( protein_data,
                                                  cosmic_data, biogrid_data,
                                                  prot_acc, central_pos[prot_acc],
                                                  mutations[prot_acc],
                                                  nodes, id_counter, id_dict,
                                                  start_pos)

        edges, id_counter = connect_protein_sequence(prot_acc, nodes, edges,
                                                     id_counter)

        pfams[prot_acc] = pfams_of_acc_from_MongoDB(protein_data, prot_acc)
        elms[prot_acc] = elms_of_acc_from_MongoDB(protein_data, prot_acc)


    elm_dom = get_elm_dom_from_MongoDB(elm_int_data)
    elm_info = get_elm_info(elm_classes)

    elm_nodes = defaultdict(list)
    lines = []
    for pair in itertools.combinations(target_prots, 2):
        ac_a, ac_b = pair
        gene_a = protein_data.find_one({"uniprot_acc": ac_a},
                                       {"_id": 0, "gene": 1})["gene"]
        gene_b = protein_data.find_one({"uniprot_acc": ac_b},
                                       {"_id": 0, "gene": 1})["gene"]

        if [ac_a, ac_b] in custom_pairs or [ac_b, ac_a] in custom_pairs:
            edges.append(   {"group" : "edges",
                             "data" : {  "id" : id_counter,
                                         "source" : id_dict[ac_a]["main"],
                                         "target" : id_dict[ac_b]["main"],
                                         "role" : "user_interaction",
                                         "ds": "User input"
                                      }
                            })
            id_counter += 1

        ## BioGRID interactions
        evidence = get_Biogrid_from_MongoDB(biogrid_data, gene_a, gene_b)
        biogrid_ints = list(evidence["Low"])+list(evidence["High"])

        if len(biogrid_ints) > 0:
            edges.append(   {"group" : "edges",
                             "data" : {  "id" : id_counter,
                                         "source" : id_dict[ac_a]["main"],
                                         "target" : id_dict[ac_b]["main"],
                                         "role" : "prot_prot_interaction",
                                         "ds": "BioGRID",
                                         "low": list(evidence["Low"]),
                                         "high": list(evidence["High"])
                                      }
                            })
            id_counter += 1

            line = [gene_a, ac_a, gene_b, ac_b, "PROT::PROT", "", "",
                    "; ".join(list(biogrid_ints)), "BioGRID"]
            lines.append(line)


        for pfam_pair in itertools.product(pfams[ac_a], pfams[ac_b]):
            pfam_a, pfam_b = pfam_pair

            ## 3did interactions
            pdbs = get_3did_from_MongoDB(db3did_data, pfam_a, pfam_b)

            if pdbs:
                for source in id_dict[ac_a][pfam_a]:
                    for target in id_dict[ac_b][pfam_b]:
                        if source == target:
                            continue
                        edges.append({ "group" : "edges",
                                       "data" : { "id" : id_counter,
                                                  "source" : source,
                                                  "target" : target,
                                                  "role" : "DOM_interaction",
                                                  "ds": "3did",
                                                  "links": pdbs
                                                 }
                                    })
                        id_counter += 1

                line = [gene_a, ac_a, gene_b, ac_b, "DOM::DOM",
                                  pfam_a, pfam_b, pdbs, "3did"]
                lines.append(line)

            ## domain propensities
            prop_lo = domain_propensities_from_MongoDB(dom_prop_data, pfam_a, pfam_b,
                                             obs_min=4, lo_min=2.0, ndom_min=4)
            if prop_lo:
                for source in id_dict[ac_a][pfam_a]:
                    for target in id_dict[ac_b][pfam_b]:
                        if source == target:
                            continue
                        edges.append({ "group" : "edges",
                                       "data" : { "id" : id_counter,
                                                  "source" : source,
                                                  "target" : target,
                                                  "role" : "iDOM_interaction",
                                                  "ds": "Statistical Prediction"
                                                 }
                                    })
                        id_counter += 1

                line = [gene_a, ac_a, gene_b, ac_b, "iDOM::iDOM",
                                  pfam_a, pfam_b,
                                  str(prop_lo), "Statistical Prediction"]
                lines.append(line)

        ## ELM-domain interactions
        for ac1, gene1, ac2, gene2 in zip([ac_a, ac_b], [gene_a, gene_b],
                                          [ac_b, ac_a], [gene_b, gene_a]):

            for elm_name in elms[ac1]:
                for pfam in pfams[ac2]:
                    # only_prts = get_elm_dom_from_MongoDB(elm_int_data, elm_name, pfam)
                    # if only_prts:
                    #     if only_prts[0].strip() and gene2 not in only_prts:
                    #         continue
                    if elm_name in elm_dom and pfam in elm_dom[elm_name]:
                        only_prts = elm_dom[elm_name][pfam]
                        if only_prts:
                            only_prts = only_prts.split(",")
                            if only_prts[0].strip() and gene2 not in only_prts:
                                    continue

                        ## Add ELM nodes if they don't exist
                        if elm_name not in id_dict[ac1]:
                            elm_des, elm_regex = "", ""
                            if elm_name in elm_info:
                                elm_des, elm_regex = elm_info[elm_name][1:]

                            for elm in elms[ac1][elm_name]:
                                elm_acc = elm["acc"]
                                start = elm["start"]
                                end = elm["end"]
                                length = end-start
                                elm_seq = elm["seq"]

                                nodes.append(
                                    { "group" : "nodes",
                                      "data" :
                                          { "id" : id_counter,
                                            "parent" : id_dict[ac1]["main"],
                                            "role" : "elm",
                                            "label" : elm_name,
                                            "acc" : elm_acc,
                                            "des" : elm_des,
                                            "regex" : elm_regex,
                                            "start" : str(start),
                                            "end" : str(end),
                                            "seq" : elm_seq,
                                            "length": str(length),
                                            "protein" : ac1
                                          },
                                            "position" : {
                                                "x" : start_pos[ac1][0]+start+(length/2),
                                                "y" : start_pos[ac1][1]
                                          }
                                })
                                id_dict[ac1][elm_name].append(id_counter)
                                id_counter += 1

                            ## Add interaction edge
                            for source in id_dict[ac1][elm_name]:
                                for target in id_dict[ac2][pfam]:
                                    edges.append(
                                            { "group" : "edges",
                                              "data" :
                                                { "id" : id_counter,
                                                  "source" : source,
                                                  "target" : target,
                                                  "role" : "ELM_interaction",
                                                  "ds": "ELM"
                                                }
                                            })
                                    id_counter += 1

                            line = [gene1, ac1, gene2, ac2, "ELM::DOM",
                                              elm_name, pfam, "", "ELM"]
                            lines.append(line)

        ## InterPreTS interactions
        for hit in get_Interprets_from_MongoDB(iprets_data, gene_a, gene_b):
            pdb_a, eval_a, pcid_a, start_a, end_a, pdb_start_a, pdb_end_a = hit[:7]
            pdb_b, eval_b, pcid_b, start_b, end_b, pdb_start_b, pdb_end_b, z = hit[7:]
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
                length = end-start
                nodes.append(
                    { "group" : "nodes",
                      "data" :
                          { "id": id_counter,
                            "parent": id_dict[ac]["main"],
                            "role": "iprets",
                            "label": label,
                            "pdb": pdb,
                            "pdb_start": pdb_start,
                            "pdb_end": pdb_end,
                            "start": start,
                            "end": end,
                            "length": length,
                            "eval": str(e_val),
                            "pcid": pcid,
                            "color": color,
                            "protein": ac
                          },
                            "position" : {
                                "x" : start_pos[ac][0]+start+(length/2),
                                "y" : start_pos[ac][1]
                          }
                })
                id_dict[ac][label] = id_counter
                id_counter += 1

            ## Add interaction edge
            edges.append(
                { "group" : "edges",
                  "data" :
                    { "id" : id_counter,
                      "source" : id_dict[ac_a][label_a],
                      "target" : id_dict[ac_b][label_b],
                      "role" : "INT_interaction",
                      "pdb": pdb,
                      "color": color,
                      "z-score": z
                    }
                })
            id_counter += 1

            line = [gene_a, ac_a, gene_b, ac_b, "InterPreTS",
                   label_a, label_b, str(z), "InterPreTS prediction"]

            lines.append(line)

    ## Color domain & LMs nodes
    nodes = color_regions(nodes, palette="custom1")

    ## Print graph as JSON file
    graph_elements = nodes + edges
    with open(graph_out, "w") as output:
        json.dump(graph_elements, output, indent=4)

    ## Print interactions as JSON file
    # for l in lines:
    #     print l
    int_table = {}
    int_table["columns"] = ["#Gene(A)","Accession(A)","Gene(B)","Accession(B)","Type",
                            "F(A)","F(B)","Info", "Source"]
    int_table["index"] = range(len(lines))
    int_table["data"] = lines

    with open(ints_out,"w") as output:
        json.dump(int_table, output)
        # output.write(str(int_table))
