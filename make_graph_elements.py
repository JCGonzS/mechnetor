#!/usr/bin/env python

""" Creates cytoscape.js network's elements v2.0

Needs a set of proteins and a file with their interactions (from int2mech.py)
"""


import sys, os, re
import gzip
import json
import csv
import copy
import math
import random
from pprint import pprint


def open_file(input_file, mode="r"):
    """ Opens file Zipped or not
    """

    if re.search(".gz$", input_file):
        infile = gzip.open(input_file, mode)
    else:
        infile = open(input_file, mode)
    return infile

def get_interaction_data(proteins, protein_data, interaction_file,
                        max_pval):
    """Gets interaction data from int2mech output file

    Puts info into a dictionary that will later be used to draw edges between
    the graph elements.
    It also "tags" ELM that have interactions so they are plotted later on.

    Function arguments:
    proteins -- list of proteins input by the user
    protein_data -- protein data dictionary (db)
    interaction_file  --  TSV file with interaction info (from int2mech.py)
    interaction_pairs  --  the list of all interacting protein pairs
    max_pval -- p-value threshold to select interactions
    """

    interaction_data = []
    with open_file(interaction_file) as f:
        for line in f:
            if line[0] != "#":
                t = line.rstrip().split("\t")
                prot_acc_a = t[1]
                prot_acc_b = t[3]
                interaction_type = t[4]
                interacting_proteins = prot_acc_a + "::" + prot_acc_b
                interacting_regions = "-"
                if interaction_type != "PROT::PROT":
                    region_a, region_b = t[5:7]
                    interacting_regions = region_a + "::" + region_b

                # Checks if proteins are target
                if (prot_acc_a not in proteins or
                    prot_acc_b not in proteins):
                    continue
                # Homodimeric interactions?
                if prot_acc_a == prot_acc_b:
                    continue
                # Checks P-value threshold.
                # FIX: right now this value is not consistent between the
                # different interaction types, so can't set same threshold 4all
                # pval = ???
                # if pval > max_pval:
                #     continue


                # Decides wether to skip this interaction and not add the info to "interaction_data" variable.
                # Will skip if...
                skip = False
                for i in range(len(interaction_data)):
                    type_check = interaction_data[i]["interaction_type"]
                    proteins_check = interaction_data[i]["interacting_proteins"]
                    protein_a_check, protein_b_check = proteins_check.split("::")
                    proteins_check_inv = protein_b_check+"::"+protein_a_check
                    regions_check = interaction_data[i]["interacting_regions"]
                    regions_check_inv = "-"
                    if type_check != "PROT::PROT":
                        region_a_check,region_b_check = regions_check.split("::")
                        regions_check_inv = region_b_check+"::"+region_a_check

                    # Same interaction already added.
                         #FIX: use sets.
                         #FIX2: actually, same domains and elms can be repeated within the same protein
                    if ((interacting_proteins == proteins_check and interacting_regions == regions_check) or
                        (interacting_proteins == proteins_check_inv and interacting_regions == regions_check_inv)):

                        if interaction_type == type_check:
                            skip = True
                            # print "1",

                        #FIX: the following conditions are commented out bc at the moment we take every type of interaction:

                        # # ... the same interaction was already added as DOM::DOM, then it has preference over iDOM::iDOM.
                        # elif interaction_type == "iDOM::iDOM" and type_check == "DOM::DOM":
                        #     skip = True
                        #     print "2",
                        #
                        # # ... the same interaction was already added as iDOM::iDOM, this one will substitute it
                        # elif interaction_type == "DOM::DOM" and type_check == "iDOM::iDOM":
                        #     interaction_data[i] = {
                        #         "interaction_type" : interaction_type,
                        #         "interacting_proteins" : interacting_proteins,
                        #         "interacting_regions" : interacting_regions
                        #     }
                        #     skip = True
                        #     print "3",
                        #
                        # # ... it's a PROT::PROT interaction, but a (i)DOM::(i)DOM for those proteins was already saved
                        # elif interaction_type == "PROT::PROT" and (type_check == "iDOM::iDOM" or type_check == "DOM:DOM"):
                        #     skip = True
                        #     print "4",
                        #
                        # # interaction was saved as PROT::PROT, but a new (i)DOM::(i)DOM would substitute it
                        # elif (interaction_type == "iDOM::iDOM" or interaction_type == "DOM::DOM") and type_check == "PROT::PROT":
                        #     interaction_data[i] = {
                        #         "interaction_type" : interaction_type,
                        #         "interacting_proteins" : interacting_proteins,
                        #         "interacting_regions" : interacting_regions
                        #     }
                        #     skip = True
                        #     print "5",

                if skip == False:

                    if interaction_type == "ELM::DOM":
                        for elm in protein_data[prot_acc_a]["elms"][region_a]:
                            elm["int"] = "yes"
                    elif interaction_type == "DOM::ELM":
                        for elm in protein_data[prot_acc_b]["elms"][region_b]:
                            elm["int"] = "yes"

                    if interaction_type == "InterPreTS":
                        try:
                            z_value = float(t[-2].split("; ")[-1])
                        except ValueError:
                            z_value = -999999

                        protein_data, interaction_data = add_interprets_data(
                                                prot_acc_a, prot_acc_b,
                                                region_a, region_b, z_value,
                                                protein_data, interaction_data)

                    else:
                        interaction_data.append({
                            "interaction_type" : interaction_type,
                            "interacting_proteins" : interacting_proteins,
                            "interacting_regions" : interacting_regions
                        })

    return protein_data, interaction_data

def add_interprets_data(prot_acc_a, prot_acc_b, region_a, region_b, z_value,
                        protein_data, interaction_data):

    pdb_a = "-".join(region_a.split(":")[0].split("|")[1:])
    pdb_b = "-".join(region_b.split(":")[0].split("|")[1:])
    start_a, end_a = region_a.split(":")[1].split("-")
    start_b, end_b = region_b.split(":")[1].split("-")

    for acc, pdb, start, end in zip([prot_acc_a, prot_acc_b], [pdb_a, pdb_b],
                                    [start_a, start_b], [end_a, end_b]):
        if "iprets" not in protein_data[acc]:
            protein_data[acc]["iprets"] = []

        protein_data[acc]["iprets"].append(
            {
                "pdb" : pdb,
                "start" : int(start),
                "end" : int(end),
                "z_value" : z_value
        })

    interaction_data.append(
        {
            "interaction_type" : "InterPreTS",
            "interacting_proteins" : prot_acc_a+"::"+prot_acc_b,
            "interacting_regions" : "::".join([pdb_a+"|"+start_a+"-"+end_a,
                                               pdb_b+"|"+start_b+"-"+end_b]),
            "z_value" : z_value
    })

    return protein_data, interaction_data

# def get_lmd2_data(lmd2_file, interaction_pairs, interaction_data, protein_data):
#     with open_file(lmd2_file) as f:
#         for line in csv.reader(f, delimiter='\t'):
#             if not line[0].startswith("#"):
#                 prot_acc_a, prot_acc_b = line[0], line[3]
#                 gene_a, gene_b = line[2], line[5]
#                 interaction_type = line[6]
#                 interacting_proteins = prot_acc_a + "::" + prot_acc_b
#                 lm_name, lm_start, lm_end = line[7], line[8], line[9]
#                 lm = lm_name+"|"+lm_start+"-"+lm_end
#
#                 if ((gene_a,gene_b) not in interaction_pairs
#                     and (gene_b,gene_a) not in interaction_pairs
#                     and (prot_acc_a,prot_acc_b) not in interaction_pairs
#                     and (prot_acc_b,prot_acc_a) not in interaction_pairs):
#                     continue
#
#                 if interaction_type == "LMD2":
#
#                     interaction_data.append({
#                         "interaction_type" : interaction_type,
#                         "interacting_proteins" : interacting_proteins,
#                         "interacting_regions" : prot_acc_a+"::"+lm
#                     })
#
#                     if prot_acc_b in protein_data:
#                         if "newLM" in protein_data[prot_acc_b]:
#                             protein_data[prot_acc_b]["newLM"][lm] = {
#                                 "name" : lm_name,
#                                 "start" : lm_start,
#                                 "end" : lm_end
#                             }
#                         else:
#                             protein_data[prot_acc_b]["newLM"] = {
#                                 lm : {
#                                     "name" : lm_name,
#                                     "start" : lm_start,
#                                     "end" : lm_end
#                                 }
#                             }
#                     else:
#                         protein_data[prot_acc_b] = {
#                             "newLM" : {
#                                 lm : {
#                                     "name" : lm_name,
#                                     "start" : lm_start,
#                                     "end" : lm_end
#                                 }
#                             }
#                         }
#     return interaction_data, protein_data

# def get_mutation_data(input_file, protein_data):
#     if os.path.isfile(input_file):
#         with open_file(input_file) as f:
#             for line in f:
#                 if line.strip() and "/" in line:
#                     info, gene = line.rstrip().split()
#                     prot_acc, mutation = info.split("/")
#                     location = int(re.search(r'\d+', mutation).group())
#                     residue_old = mutation[0]
#                     residue_new = mutation[len(mutation)-1]
#                     mutation = residue_old + "->" +  residue_new
#                     if prot_acc in protein_data:
#                         if "mutations" in protein_data[prot_acc]:
#                             if location in protein_data[prot_acc]["mutations"]:
#                                 old = protein_data[prot_acc]["mutations"][location]
#                                 info = old.split("->")
#                                 last_residues = info[1].split(",")
#                                 if residue_new in last_residues:
#                                     continue
#                                 protein_data[prot_acc]["mutations"][location] = old + "," + residue_new
#                             else:
#                                 protein_data[prot_acc]["mutations"][location] = mutation
#                         else:
#                             protein_data[prot_acc]["mutations"] = {
#                                 location: mutation
#                             }
#                     else:
#                         protein_data[prot_acc] = {
#                             "gene": gene,
#                             "mutations": {
#                                 location: mutation
#                             }
#                         }
#     return protein_data

def create_positions(proteins):
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
        x += 800
        y += 225
        central_pos[prot_acc] = (x, y)
        i = i+2

    return central_pos

def add_protein(prot_acc, protein_data, mutations, central_pos,
                graph_elements, id_counter):
    """Adds a protein to the graph elements dictionary in "cytoscape-format"

    This function first adds three nodes per protein to the graph:
    the first node is the parent node and represents the whole protein, the
    other two nodes are the children of the first node and mark the start and
    end of the protein.

    Function also calls others to add the domains, ELMs, mutations, PTMs of
    the protein, which are checked to see if they are located at the start or
    at the end.

    Function arguments:
    prot_acc -- the protein accession/identifier
    protein_data -- protein data dictionary (db) for this protein
    mutations -- mutation list for this protein
    central_pos -- x, y coordinates for this protein
    graph_elements --
    id_counter  --
    """

    protein_id = copy.deepcopy(id_counter)

    ## Add protein central node
    graph_elements.append({
    "group" : "nodes",
    "data" : {
    "id" : protein_id,
    "label" : protein_data["gene"],
    "role" : "whole",
    "parent" : protein_id,
    "protein" : prot_acc,

    }
    })
    id_counter += 1

    start_x = central_pos[0] - protein_data["length"]/2
    start_y = central_pos[1]
    end_x = central_pos[0] + protein_data["length"]/2
    end_y = start_y
    start_found = False
    end_found = False

    start_found, end_found, graph_elements, id_counter = add_domains(
                                    prot_acc, protein_data, protein_id,
                                    start_x, start_y, start_found, end_found,
                                    graph_elements, id_counter)

    ## adds ELMs and LMD2-LMs
    start_found, end_found, graph_elements, id_counter = add_elms(
                                    prot_acc, protein_data, protein_id,
                                    start_x, start_y, start_found, end_found,
                                    graph_elements, id_counter)

    start_found, end_found, graph_elements, id_counter = add_interprets(
                                    prot_acc, protein_data, protein_id,
                                    start_x, start_y, start_found, end_found,
                                    graph_elements, id_counter)


    graph_elements, id_counter = add_mutations(prot_acc, protein_data,
                                               mutations, protein_id,
                                               start_x, start_y,
                                               graph_elements, id_counter)

    graph_elements, id_counter = add_ptms(prot_acc, protein_data, protein_id,
                                          start_x, start_y,
                                          graph_elements, id_counter)


    ## Add protein start node (if not occupied by another element)
    if not start_found:
        graph_elements.append({
            "group" : "nodes",
            "data" : {
                "id" : id_counter,
                "label" : "0",
                "role" : "start-end",
                "parent" : protein_id,
                "protein" : prot_acc
            },
            "position" : {
                "x" : start_x,
                "y" : start_y
            }
        })
        id_counter += 1

    ## Add protein end node (if not occupied by another element)
    if not end_found:
        graph_elements.append({
            "group" : "nodes",
            "data" : {
                "id" : id_counter,
                "label" : str(protein_data["length"]),
                "role" : "start-end",
                "parent" : protein_id,
                "protein" : prot_acc
            },
            "position" : {
                "x" : end_x,
                "y" : end_y
            }
        })
        id_counter += 1

    return graph_elements, id_counter

def add_domains(prot_acc, protein_data, parent_id, start_x, start_y,
                start_found, end_found, graph_elements, id_counter):
    """Add protein's Pfam domains

    Adds three nodes per domain: central, start and end
    """

    for pfam_name in protein_data["pfams"]:
        for instance in protein_data["pfams"][pfam_name]:
            start = instance["start"]
            end = instance["end"]
            element_id = copy.deepcopy(id_counter)
            if start == 0:
                start_found = True
            if end == int(protein_data["length"]):
                end_found = True

            # Add central node
            graph_elements.append({
                "group" : "nodes",
                "data" : {
                    "id" : element_id,
                    "label" : pfam_name,
                    "role" : "domain",
                    "parent" : parent_id,
                    "protein" : prot_acc
                }
            })
            id_counter += 1

            # Add start node
            graph_elements.append({
                "group" : "nodes",
                "data" : {
                    "id" : id_counter,
                    "label" : str(start),
                    "role" : "dom_pos",
                    "parent" : element_id,
                    "protein" : prot_acc
                },
                "position" : {
                    "x" : start_x + start,
                    "y" : start_y
                }
            })
            id_counter += 1

            # Add end node
            graph_elements.append({
                "group" : "nodes",
                "data" : {
                    "id" : id_counter,
                    "label" : str(end),
                    "role" : "dom_pos",
                    "parent" : element_id,
                    "protein" : prot_acc
                },
                "position" : {
                    "x" : start_x + end,
                    "y" : start_y
                }
            })
            id_counter += 1

    return start_found, end_found, graph_elements, id_counter

def add_elms(prot_acc, protein_data, parent_id, start_x, start_y,
             start_found, end_found, graph_elements, id_counter):
    """Add protein's ELMs

    Adds three nodes per ELM: central, start and end
    """

    for node_role in ["elms", "newLM"]:
        if node_role in protein_data:

            for lm_name in protein_data[node_role]:
                for lm in protein_data[node_role][lm_name]:
                    if "int" in lm:
                        lm_acc = lm["acc"]
                        start = lm["start"]
                        end = lm["end"]
                        element_id = copy.deepcopy(id_counter)

                        if start == 0:
                            start_found = True
                        if end == int(protein_data["length"]):
                            end_found = True

                        # Add central node
                        graph_elements.append({
                            "group" : "nodes",
                            "data" : {
                                "id" : element_id,
                                "label" : lm_name,
                                "start" : str(start),
                                "end" : str(end),
                                "role" : node_role,
                                "parent" : parent_id,
                                "protein" : prot_acc
                            }
                        })
                        id_counter += 1

                        # Add start node
                        graph_elements.append({
                            "group" : "nodes",
                            "data" : {
                                "id" : id_counter,
                                "label" : str(start),
                                "role" : "position",
                                "parent" : element_id,
                                "protein" : prot_acc
                            },
                            "position" : {
                                "x" : start_x + start,
                                "y" : start_y
                            }
                        })
                        id_counter += 1

                        # Add end node
                        graph_elements.append({
                            "group" : "nodes",
                            "data" : {
                                "id" : id_counter,
                                "label" : str(end),
                                "role" : "position",
                                "parent" : element_id,
                                "protein" : prot_acc
                            },
                            "position" : {
                                "x" : start_x + end,
                                "y" : start_y
                            }
                        })
                        id_counter += 1

    return start_found, end_found, graph_elements, id_counter

def add_interprets(prot_acc, protein_data, parent_id, start_x, start_y,
                   start_found, end_found, graph_elements, id_counter):
    """Adds Interprets regions

    Adds three nodes per regionn: central, start and end
    """

    if "iprets" in protein_data:
        for iprets in protein_data["iprets"]:
            pdb = iprets["pdb"]
            start = iprets["start"]
            end = iprets["end"]
            z_value = iprets["z_value"]
            element_id = copy.deepcopy(id_counter)

            if start == 0:
                start_found = True
            if end == int(protein_data["length"]):
                end_found = True

            # Add central node
            color = color_from_zvalue(z_value)
            graph_elements.append({
                "group" : "nodes",
                "data" : {
                    "id" : element_id,
                    "label" : pdb+"|"+str(start)+"-"+str(end),
                    "role" : "iprets",
                    "color" : color,
                    "parent" : parent_id,
                    "protein" : prot_acc
                }
            })
            id_counter += 1

            # Add start node
            graph_elements.append({
                "group" : "nodes",
                "data" : {
                    "id" : id_counter,
                    "label" : str(start),
                    "role" : "position",
                    "parent" : element_id,
                    "protein" : prot_acc
                },
                "position" : {
                    "x" : start_x + start,
                    "y" : start_y
                }
            })
            id_counter += 1

            # Add end node
            graph_elements.append({
                "group" : "nodes",
                "data" : {
                    "id" : id_counter,
                    "label" : str(end),
                    "role" : "position",
                    "parent" : element_id,
                    "protein" : prot_acc
                },
                "position" : {
                    "x" : start_x + end,
                    "y" : start_y
                }
            })
            id_counter += 1

    return start_found, end_found, graph_elements, id_counter

def add_mutations(prot_acc, protein_data, mutations, parent_id,
                  start_x, start_y,
                  graph_elements, id_counter):

    for pos in protein_data["cosmic_muts"]:
        label = ";".join(protein_data["cosmic_muts"][pos])

        if int(pos)== 0:
            start_found = True
        if int(pos) == int(protein_data["length"]):
            end_found = True

        graph_elements.append({
            "group" : "nodes",
            "data" : {
                "id" : id_counter,
                "label" : label,
                "role" : "mutation",
                "parent" : parent_id,
                "protein" : prot_acc
            },
            "position" : {
                "x" : start_x + int(pos),
                "y" : start_y
            }
        })
        id_counter += 1

    for pos in mutations:
        label = ";".join(list(mutations[pos]))

        if int(pos)== 0:
            start_found = True
        if int(pos) == int(protein_data["length"]):
            end_found = True

        graph_elements.append({
            "group" : "nodes",
            "data" : {
                "id" : id_counter,
                "label" : label,
                "role" : "mutation",
                "parent" : parent_id,
                "protein" : prot_acc
            },
            "position" : {
                "x" : start_x + int(pos),
                "y" : start_y
            }
        })
        id_counter += 1

    return graph_elements, id_counter

def add_ptms(prot_acc, protein_data, parent_id, start_x, start_y,
             graph_elements, id_counter):

    for ptm, role in zip(["phosphorylation", "acetylation"],
                         ["pp_mod", "ac_mod"]):

        for pos in protein_data[ptm]:
            residue = protein_data[ptm][pos]

            if int(pos)== 0:
                start_found = True
            if int(pos) == int(protein_data["length"]):
                end_found = True

            graph_elements.append({
                "group" : "nodes",
                "data" : {
                    "id" : id_counter,
                    "label" : residue + str(pos),
                    "role" : role,
                    "parent" : parent_id,
                    "protein" : prot_acc
                },
                "position" : {
                    "x" : start_x + int(pos),
                    "y": start_y
                }
            })
            id_counter += 1

    return graph_elements, id_counter

def add_all_proteins(proteins, protein_data, mutations, central_pos):
    """Add all proteins to the graph elements list

    This function adds, individually, all proteins that need to be displayed
    to the graph elements dictionary in "cytoscape-format".
    """

    graph_elements = []
    id_counter = 0

    for prot_acc in sorted(list(proteins)):
        try:
            graph_elements, id_counter = add_protein(
                                    prot_acc, protein_data[prot_acc],
                                    mutations[prot_acc], central_pos[prot_acc],
                                    graph_elements, id_counter)
        except KeyError as err:
            print "ERROR:", prot_acc + ": MISSING DATA!"

    return graph_elements, id_counter

def color_from_zvalue(z_value):
    if z_value == -999999:
        color = "black"
    elif z_value < 0 and z_value > -999999:
        color = "blue"
    elif z_value >= 0 and z_value < 1.6:
        color = "yellow"
    elif z_value >= 1.6 and z_value < 2.1:
        color = "orange"
    elif z_value >= 2.1:
        color = "red"

    return color

def add_interaction(graph_elements, id_counter, source, target,
                    int_type, z_value=0):
    """Add edge between source and target elements

    Function arguments:
    elements  --  the list of elements(all nodes and edges) in
                  "cytoscape-format"
    id_counter  --  the counter that is used as the id when a new elements is added
                 to the elements list
    source  --  the elements list id of the first interacting protein
    target  --  the elements list id of the second interacting protein
    int_type  --  type of interaction (DOM::DOM, iDOM::iDOM, DOM::ELM, PROT::PROT, INT::INT)
    """

    if int_type == "InterPreTS":
        color = color_from_zvalue(z_value)

        for so in source:
            for ta in target:
                if so != ta:
                    graph_elements.append({
                        "group" : "edges",
                        "data" : {
                            "id" : id_counter,
                            "source" : so,
                            "target" : ta,
                            "role" : "INT_interaction",
                            "color" : color
                        }
                    })
                    id_counter += 1
    else:
        if int_type == "DOM::DOM":
            role = "DOM_interaction"
        elif int_type == "iDOM::iDOM":
            role = "iDOM_interaction"
        elif int_type == "DOM::ELM" or int_type == "ELM::DOM":
            role = "ELM_interaction"
        elif int_type == "PROT::PROT":
            role = "prot_prot_interaction"
        elif int_type == "LMD2":
            role = "LMD2_interaction"
        else:
            role = int_type

        for so in source:
            for ta in target:
                if so != ta:
                    graph_elements.append({
                        "group" : "edges",
                        "data" : {
                            "id" : id_counter,
                            "source" : so,
                            "target" : ta,
                            "role" : role
                        }
                    })
                    id_counter += 1

    return graph_elements, id_counter

def add_all_interactions(interaction_data, graph_elements, id_counter):
    """Add all interactions that need to be displayed to the elements list

    This function adds all interactions specified in the input file to the
    elements list in "cytoscape-format".

    Function arguments:
    interaction_data  --  list of dicts that contain information about each
                          interaction
    elements  --  the list of elements(all nodes and edges) in
                  "cytoscape-format"
    id_counter  --  the counter that is used as the id when a new elements is added
                 to the elements list
    """

    for interaction in interaction_data:

        # 1. Find source & target elements id's
        interaction_type = interaction["interaction_type"]
        protein_a, protein_b = interaction["interacting_proteins"].split("::")
        interacting_regions = interaction["interacting_regions"]
        if interacting_regions != "-":
            region_a, region_b = interacting_regions.split("::")

        source = []
        target = []
        z_value = 0
        for ele in graph_elements:
            if ele["group"] == "nodes":
                role = ele["data"]["role"]
                parent_protein = ele["data"]["protein"]
                label = ele["data"]["label"]
                ele_id = ele["data"]["id"]

                if interaction_type == "PROT::PROT":
                    if role == "whole":
                        if parent_protein == protein_a:
                            source.append(ele_id)
                        if parent_protein == protein_b:
                            target.append(ele_id)

                elif interaction_type == "LMD2":
                    elm = region_b.split("|")[0]
                    elm_start, elm_end = region_b.split("|")[1].split("-")
                    if role == "whole":
                        if region_a == parent_protein:
                            source.append(ele_id)
                    elif role == "newLM":
                        start = ele["data"]["start"]
                        end = ele["data"]["end"]
                        if (label == elm and parent_protein == protein_b
                            and start == elm_start and end == elm_end):
                            target.append(ele_id)

                else:
                    if interaction_type == "InterPreTS":
                        z_value = interaction["z_value"]
                    if role in ["domain", "elms", "newLM", "iprets"]:
                        if label == region_a and parent_protein == protein_a:
                            source.append(ele_id)
                        if label == region_b and parent_protein == protein_b:
                            target.append(ele_id)

        # 2. Add interaction to graph elements as an edge
        graph_elements, id_counter = add_interaction(graph_elements,
                                                     id_counter,
                                                     source, target,
                                                     interaction_type,
                                                     z_value)

    return graph_elements, id_counter

def connect_protein_sequence(proteins, graph_elements, id_counter):
    """Add edges representing the amino acid chain of the proteins

    This function adds all the edges to the elements list that will later
    represent the amino acid chain of the proteins. First alls nodes belonging
    to the same protein are grouped together and then these are sorted according
    to their position on the protein to ensure the right order of connections
    (going from beginning to the end of the protein). After that, the edges are
    added connecting the nodes on the order that they are in after sorting the
    lists.

    Function arguments:
    proteins  --  the list containing all the proteins that need to be displayed
    elements  --  the list of elements(all nodes and edges) in
                  "cytoscape-format"
    id_counter  --  the counter that is used as the id when a new elements is added
                 to the elements list
    """

    # 1. Decide which nodes form the protein sequence
    for prot_acc in sorted(list(proteins)):
        prot_positions = []
        for ele in graph_elements:
            if ele["group"] == "nodes":
                role = ele["data"]["role"]
                parent_prot = ele["data"]["protein"]
                position = ele["data"]["label"]
                ele_id = ele["data"]["id"]
                if (role in ["position", "dom_pos", "start-end"] and
                    parent_prot == prot_acc):
                    prot_positions.append((ele_id, int(position)))

    # 2. Join those nodes by edges
        prot_positions.sort(key=lambda x: x[1])
        for i in range(len(prot_positions)-1):

            graph_elements.append({
                "group" : "edges",
                "data" : {
                    "id" : id_counter,
                    "source" : prot_positions[i][0],
                    "target" : prot_positions[i+1][0],
                    "role" : "protein_sequence"
                }
            })
            id_counter += 1

    return graph_elements, id_counter

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

    return c

def color_regions(graph_elements, palette=""):
    """Gets colours for domains and ELMs graph elements

    Function arguments:
     graph_elements --
     palette -- if specified a custom palette of colours will be used. See
        "get_palette" or options. If let empty random colours will be used.
    """

    custom_colors = get_palette(palette)
    color_map = {}
    counter = 0
    for ele in graph_elements:
        if ele["group"] == "nodes":
            role = ele["data"]["role"]
            label = ele["data"]["label"]
            if role in ["domain", "elms", "newLM"]:
                if label not in color_map:
                    if palette != "" and counter < len(custom_colors):
                        color_map[label] = custom_colors[counter]
                        counter += 1
                    else:
                        color_map[label] = get_random_color(color_map.values())

                ele["data"]["color"] = color_map[label]

    return graph_elements

def main(proteins, protein_data, mutations,
         interaction_file, lmd2_file="",  output_file="graph_elements.json",
         max_pval=999):

    protein_data, interaction_data = get_interaction_data(proteins,
                                    protein_data, interaction_file, max_pval)

    if lmd2_file != "":
        interaction_data, protein_data = get_lmd2_data(lmd2_file,
                                                       interaction_pairs,
                                                       interaction_data,
                                                       protein_data)

    ## Make graph elements
    central_pos = create_positions(proteins)

    graph_elements, id_counter = add_all_proteins(proteins, protein_data,
                                                  mutations, central_pos)
    graph_elements, id_counter = connect_protein_sequence(proteins,
                                                    graph_elements, id_counter)
    graph_elements, id_counter = add_all_interactions(interaction_data,
                                                    graph_elements, id_counter)
    graph_elements = color_regions(graph_elements, palette="autumn")

    ## Print graph as JSON file
    with open(output_file, "w") as output:
        json.dump(graph_elements, output, indent=4)
