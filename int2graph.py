#!/usr/bin/env python


import sys, re
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


def add_protein(protein_data, prot_acc, central_pos, mutations,
                graph_elements, id_counter):

    """Adds a protein to the graph elements dictionary in "cytoscape-format"

    This function first adds three nodes per protein to the graph:
    the first node is the parent node and represents the whole protein, the
    other two nodes are the children of the first node and mark the start and
    end of the protein.

    Function also calls others to add the domains, ELMs, mutations, PTMs of
    the protein, which are checked to see if they are located at the start or
    at the end.
    """

    protein_id = copy.deepcopy(id_counter)

    cursor = protein_data.find_one( { "uniprot_acc": prot_acc },
                        { "_id": 0, "cosmic_muts": 0, "data_class": 0 })

    ## Add protein central node
    graph_elements.append({
        "group" : "nodes",
        "data" : {
            "id" : protein_id,
            "parent" : protein_id,
            "role" : "whole",
            "label" : cursor["gene"],
            "uni_id" : cursor["uniprot_id"],
            "des": cursor["description"],
            "length": cursor["length"],
            "protein" : prot_acc
            }
        })
    id_counter += 1

    start_x = central_pos[0] - cursor["length"]/2
    start_y = central_pos[1]
    end_x = central_pos[0] + cursor["length"]/2
    end_y = start_y
    start_found = False
    end_found = False

    start_found, end_found, graph_elements, id_counter = add_domains(
                                    prot_acc, protein_id, cursor,
                                    start_x, start_y, start_found, end_found,
                                    graph_elements, id_counter)

    ## adds ELMs and LMD2-LMs
    start_found, end_found, graph_elements, id_counter = add_elms(
                                    prot_acc, protein_id, cursor, elms_to_show,
                                    start_x, start_y, start_found, end_found,
                                    graph_elements, id_counter)

    # start_found, end_found, graph_elements, id_counter = add_interprets(
    #                                 prot_acc, protein_data, protein_id,
    #                                 start_x, start_y, start_found, end_found,
    #                                 graph_elements, id_counter)

    graph_elements, id_counter = add_ptms( prot_acc, protein_id, cursor,
                                           start_x, start_y,
                                           graph_elements, id_counter)

    graph_elements, id_counter = add_mutations(prot_acc, protein_id, cursor,
                                               mutations, start_x, start_y,
                                               graph_elements, id_counter)



    ## Add protein start node (if not occupied by another element)
    if not start_found:
        graph_elements.append({
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

    ## Add protein end node (if not occupied by another element)
    if not end_found:
        graph_elements.append({
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

    return graph_elements, id_counter

def add_domains(prot_acc, parent_id, cursor, start_x, start_y,
                start_found, end_found, graph_elements, id_counter):

    """Add protein's Pfam domains

    Adds three nodes per domain: central, start and end
    """

    for domain in cursor["pfams"]:
        print prot_acc, domain["name"]
        element_id = copy.deepcopy(id_counter)
        start = domain["start"]
        end = domain["end"]
        if start == 0:
            start_found = True
        if end == int(cursor["length"]):
            end_found = True

        # Add central node
        graph_elements.append({
            "group" : "nodes",
            "data" : {
                "id" : element_id,
                "parent" : parent_id,
                "role" : "domain",
                "label" : domain["name"],
                "acc" : domain["acc"].split(".")[0],
                "start": str(start),
                "end": str(end),
                "protein" : prot_acc
            }
        })
        id_counter += 1

        # Add start node
        graph_elements.append({
            "group" : "nodes",
            "data" : {
                "id" : id_counter,
                "parent" : element_id,
                "role" : "dom_pos",
                "label" : str(start),
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
                "parent" : element_id,
                "role" : "dom_pos",
                "label" : str(end),
                "protein" : prot_acc
            },
            "position" : {
                "x" : start_x + end,
                "y" : start_y
            }
        })
        id_counter += 1

    return start_found, end_found, graph_elements, id_counter


@line_profile
def main(target_prots, protein_data, mutations,
        biogrid_data, iprets_data, db3did_data, dom_prop_data, elm_int_dom_file,
        output_file="", max_prots=""):


    central_pos = central_positions_layout(target_proteins)
    nodes, edges = [], []
    id_counter = 0

    for prot_acc in sorted(list(target_prots)):
        nodes, id_counter = add_protein( protein_data,
                                                  prot_acc, central_pos[prot_acc],
                                                  mutations[prot_acc],
                                                #   elms_to_show[prot_acc],
                                                  nodes, id_counter)


    nodes, id_counter = add_all_proteins(protein_data, proteins, central_pos,
                                            mutations, elms_to_show)
    edges, id_counter = connect_protein_sequence(proteins, nodes, id_counter)

    for pair in itertools.combinations(target_prots, 2):
        ac_a, ac_b = pair
        gene_a = protein_data.find_one({"uniprot_acc": ac_a},
                                       {"_id": 0, "gene": 1})["gene"]
        gene_b = protein_data.find_one({"uniprot_acc": ac_b},
                                       {"_id": 0, "gene": 1})["gene"]


        ## BioGRID interactions
        info = get_Biogrid_from_MongoDB(biogrid_data, gene_a, gene_b)
