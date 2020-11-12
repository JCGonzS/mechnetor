#!/usr/bin/env python

import sys, re, os
import math, random, copy, pprint, json
import gzip, itertools, pymongo
from collections import defaultdict
from flask_debugtoolbar_lineprofilerpanel.profile import line_profile

def get_layout_positions(proteins):
    """Creates x/y coordinates of proteins to be displayed

    This function uses the parameterized archimedean spiral equation to
    generate the x/y coordinates of the initial positions of the proteins in
    the network.
    """

    pos = {}
    i = 0
    for uni_ac in sorted(list(proteins)):
        angle = i
        x = (1+100*angle) * math.cos(angle)
        y = (1+100*angle) * math.sin(angle)
        x += 650
        y += 100
        pos[uni_ac] = (x, y)
        i = i+2

    return pos

def add_protein_main(nodes, id_n, id_dict, uni_ac, data, ini_pos):

    id_n += 1
    prot_id = copy.deepcopy(id_n)
    id_dict[uni_ac] = defaultdict(list)
    id_dict[uni_ac]["main"] = prot_id
    gene_name = data["genes"][0]
    if len(data["genes"]) > 1:
        gene_name += "(& more)"
    nodes.append({
        "group": "nodes",
        "data": {
            "id": prot_id,
            # "parent": prot_id,
            "role": "protein_main",
            "label": gene_name,
            "uni_id": data["uni_id"],
            "biogrid_id": data["biogrid_id"],
            "des": data["description"],
            "length": data["length"],
            "protein": uni_ac,
            "display": "element"
        }
    })

    id_n += 1
    nodes.append({
        "group": "nodes",
        "data": {
            "id": id_n,
            "parent": prot_id,
            "role": "protein_seq",
            "length": data["length"]
        },
        "position": {
            "x": ini_pos[0],
            "y": ini_pos[1]
        }
    })

    return nodes, id_n, id_dict, prot_id

def add_fasta_main(nodes, id_n, id_dict, header, seq, link, blast, ini_pos):

    id_n += 1
    prot_id = copy.deepcopy(id_n)
    id_dict[header] = defaultdict(list)
    id_dict[header]["main"] = prot_id

    nodes.append({
        "group": "nodes",
        "data": {
            "id": prot_id,
            # "parent": prot_id,
            "role": "protein_main",
            "label": header,
            "link": link,
            "length": len(seq),
            "blast": blast,
            "protein": header,
            "display": "element"
        }
    })

    id_n += 1
    nodes.append({
        "group": "nodes",
        "data": {
            "id": id_n,
            "parent": prot_id,
            "role": "protein_seq",
            "length": len(seq)
        },
        "position": {
            "x": ini_pos[0],
            "y": ini_pos[1]
        }
    })

    return nodes, id_n, id_dict, prot_id

def add_domains(nodes, id_n, id_dict, id_coords, id_muts, prot_id, uni_ac,
                prot_center, ini_pos, pfams, pfam_info, muts):

    """Add protein's Pfam domains

    Adds domain as a single node whose width is proportional to domain's length
    """

    pfams2 = defaultdict(list)
    for domain in pfams:
        pfams2[domain["e-val"]].append(domain) # in case there are more
                                               # than one match with same
                                               # e-value (very unlikely)
    pfam_set = set()
    for e_val in sorted(pfams2, reverse=True):
        for domain in pfams2[e_val]:
            pfam_ac = domain["acc"].split(".")[0]
            if pfam_ac not in pfam_info: # most likely old deleted pfam entries
                continue
            start = domain["start"]
            end = domain["end"]
            length = end - start
            pfam_set.add(pfam_ac)
            id_n += 1
            nodes.append({
                "group": "nodes",
                "data": {
                    "id": id_n,
                    "role": "domain",
                    "label": pfam_info[pfam_ac]["ide"],
                    "acc": pfam_ac,
                    "type": pfam_info[pfam_ac]["type"],
                    "des": pfam_info[pfam_ac]["des"],
                    "start": start,
                    "end": end,
                    "length": length,
                    "e_val": str(domain["e-val"]),
                    "parent": prot_id,
                    "protein": uni_ac
                },
                "position": {
                    "x": ini_pos[0] + start+(float(length)/2) - prot_center - 0.5,
                    "y": ini_pos[1]
                }
            })

            id_dict[uni_ac][pfam_ac].append(id_n)
            id_coords[id_n] = (str(start), str(end))
            id_muts[id_n] = muts_within_coords(uni_ac, muts, start, end)

    return nodes, id_n, id_dict, id_coords, id_muts, pfam_set

def add_ptms(nodes, id_n, prot_id, uni_ac, prot_center, ini_pos, data):
    phos = []
    for ptm_type, role, x in zip(["phosphorylation", "acetylation"],
                                 ["mod_phos", "mod_acet"],
                                 ["p","ac"]):
        for ptm in data[ptm_type]:
            pos = ptm["pos"]
            res = ptm["res"]
            if x=="p":
                phos.append(pos)
            id_n += 1
            nodes.append({
                "group": "nodes",
                "data": {
                    "id": id_n,
                    "parent": prot_id,
                    "role": role,
                    "label": res+x+str(pos),
                    "protein": uni_ac
                },
                "position": {
                    "x": ini_pos[0] + float(pos) - prot_center - 0.5,
                    "y": ini_pos[1]
                }
            })

    return nodes, id_n, phos

def parse_uni_interactors(info, prot_dict):
    interactors = set()
    for word in info.split():
        word = re.sub("[^\w]","", word)
        if len(word)>1 and word in prot_dict["AC"]:
            uni_ac = prot_dict["AC"][word][0]
            interactors.add(uni_ac)
    return list(interactors)

def add_uni_features(nodes, id_n, prot_id, prot_center, ini_pos,
                     uni_feats, uni_regions, prot_dict):

    feat_ints, feat_info = {}, {}
    for region in uni_regions:
        start, end = region["start"], region["end"]
        try:
            length = end - start
        except:
            length = 0
            if str(start).isdigit():
                end = start
            elif str(end).isdigit():
                start = end
            else:
                continue
        id_n += 1
        nodes.append({
            "group": "nodes",
            "data": {
                "id": id_n,
                "parent": prot_id,
                "role": "uni_region",
                "label": "("+str(region["start"])+"-"+str(region["end"])+") "+region["info"],
                "start": start,
                "end": end,
                "length": length
            },
            "position": {
                "x": ini_pos[0] + start+(float(length)/2) - prot_center - 0.5,
                "y": ini_pos[1]
            }
        })
        interactors = parse_uni_interactors(region["info"], prot_dict)
        if len(interactors)>0:
            feat_ints[id_n] = interactors
            feat_info[id_n] = ("uni_region", str(start), str(end), region["info"])

    uni_roles = {"VARIANT": "uni_var",
                 "MUTAGEN": "uni_mtg",
                 "METAL": "uni_metal",
                 "BINDING": "uni_binding"}
    for feat in uni_feats:
        start = int(feat["pos"].split("-")[0])
        end = int(feat["pos"].split("-")[1])
        role = uni_roles[feat["feat"]]
        if start == end:
            length = 0
            width = 2
            label = "("+str(start)+") "+feat["info"]
        else:
            length = end - start
            width = length
            label = "("+str(start)+"-"+str(end)+") "+feat["info"]
        id_n += 1
        nodes.append({
            "group": "nodes",
            "data": {
                "id": id_n,
                "parent": prot_id,
                "role": role,
                "label": label,
                "start": start,
                "end": end,
                "length": width
            },
            "position": {
                "x": ini_pos[0] + start+(float(length)/2) - prot_center - 0.5,
                "y": ini_pos[1]
            }
        })
        interactors = parse_uni_interactors(feat["info"], prot_dict)
        if len(interactors)>0:
            feat_ints[id_n] = interactors
            feat_info[id_n] = (role, str(start), str(end), feat["info"])

    return nodes, id_n, feat_ints, feat_info

def add_custom_mutations(nodes, id_n, prot_id, uni_ac, prot_len,
                         prot_center, ini_pos, mutations):
    """
    Adds user-input mutations as nodes within the proteins
    """

    for pos in mutations:
        if pos <= prot_len:
            label = ";".join(list(mutations[pos]))

            id_n += 1
            nodes.append({
                "group": "nodes",
                "data": {
                    "id": id_n,
                    "parent": prot_id,
                    "role": "mod_input",
                    "label": label,
                    "protein": uni_ac
                },
                "position": {
                    "x": ini_pos[0] + float(pos) - prot_center - 0.5,
                    "y": ini_pos[1]
                }
            })

    return nodes, id_n

def add_cosmic_mutations(nodes, id_n, prot_id, uni_ac, prot_center, ini_pos,
                         cosmic_data, min_sample_size=2):

    """
    Adds COSMIC mutations as nodes within the proteins
    Data extracted from internal MongoDB
    """

    muts = defaultdict(lambda: defaultdict(list))
    for c in cosmic_data.find( {"uni_ac": uni_ac},
                               { "_id": 0, "enst": 0 }):
        pos = re.search("(\d+)",c["aa_mut"]).group(1)
        if int(c["samples"]) >= min_sample_size:
            muts[pos][c["samples"]].append(c)

    for pos in muts:
        aa_muts, cosmic_ids, cds_muts, count, tot_count = [], [], [], [], 0
        for samples in sorted(muts[pos], reverse=True):
            for mut in muts[pos][samples]:
                aa_muts.append(mut["aa_mut"])
                cosmic_ids.append("COSM"+str(mut["cosmic_id"]))
                cds_muts.append(mut["cds_mut"])
                count.append(str(samples))
                tot_count += int(samples)

        id_n += 1
        nodes.append({
            "group": "nodes",
            "data": {
                "id": id_n,
                "parent": prot_id,
                "role": "mod_cosmic",
                "cos_id": "; ".join(cosmic_ids),
                "aa_mut": "; ".join(aa_muts),
                "cds": "; ".join(cds_muts),
                "count": "; ".join(count),
                "tot_count": tot_count
            },
            "position": {
                "x": ini_pos[0] + float(pos) - prot_center - 0.5,
                "y": ini_pos[1]
            }
        })

    return nodes, id_n

def add_elm_node(nodes, id_n, id_dict, id_coords, id_muts, prot_acc,
                 prot_center, ini_pos, elm_acc, elms, elm_info, mutations,
                 phospho_hits):

    elm_ide   = elm_info[elm_acc]["ide"]
    elm_name  = elm_info[elm_acc]["name"]
    elm_des   = elm_info[elm_acc]["des"]
    elm_regex = elm_info[elm_acc]["regex"]
    elm_prob  = elm_info[elm_acc]["prob"]

    for hit in elms[prot_acc][elm_acc]:
        hit_start   = hit["start"]
        hit_end     = hit["end"]
        hit_seq     = hit["seq"]
        hit_length  = hit_end - hit_start
        phos        = ""
        if phospho_hits != "none":
            if (hit_start, hit_end) in phospho_hits:
                phos = ["yes"]+hit["phos_pos"]
            else:
                phos = ["no"]+hit["phos_pos"]
        hit_status = hit.get("status", "")
        if hit_status == "":
            hit_status = "Unknown"
        id_n += 1
        nodes.append({
            "group": "nodes",
            "data": {
                "id": id_n,
                "parent": id_dict[prot_acc]["main"],
                "protein": prot_acc,
                "role": "elm",
                "label": elm_ide,
                "acc": elm_acc,
                "name": elm_name,
                "des": elm_des,
                "regex": elm_regex,
                "start": hit_start,
                "end": hit_end,
                "seq": hit_seq,
                "length": hit_length,
                "status": hit_status,
                "phos": phos
            },
            "position": {
                "x": ini_pos[0]+hit_start+(float(hit_length)/2)-prot_center-0.5,
                "y": ini_pos[1]
            }
        })

        id_dict[prot_acc][elm_ide].append(id_n)
        id_coords[id_n] = (str(hit_start), str(hit_end))
        id_muts[id_n] = muts_within_coords(prot_acc, mutations[prot_acc],
                                           hit_start, hit_end)

    return nodes, id_n, id_dict, id_coords, id_muts

def add_lm_node(nodes, id_n, id_dict, id_coords, id_muts, prot_acc,
                 prot_center, ini_pos, lm_acc, lms, mutations):

    for hit in lms[prot_acc][lm_acc]:
        hit_start   = hit["start"]
        hit_end     = hit["end"]
        hit_seq     = hit["seq"]
        hit_length  = hit_end - hit_start
        hit_status = hit.get("status", "")
        if hit_status == "":
            hit_status = "Unknown"
        id_n += 1
        nodes.append({
            "group": "nodes",
            "data": {
                "id": id_n,
                "parent": id_dict[prot_acc]["main"],
                "protein": prot_acc,
                "role": "3dlm",
                "label": lm_acc,
                "acc": lm_acc,
                "start": hit_start,
                "end": hit_end,
                "seq": hit_seq,
                "length": hit_length,
                "status": hit_status
            },
            "position": {
                "x": ini_pos[0]+hit_start+(float(hit_length)/2)-prot_center-0.5,
                "y": ini_pos[1]
            }
        })

        id_dict[prot_acc][lm_acc].append(id_n)
        id_coords[id_n] = (str(hit_start), str(hit_end))
        id_muts[id_n] = muts_within_coords(prot_acc, mutations[prot_acc],
                                           hit_start, hit_end)

    return nodes, id_n, id_dict, id_coords, id_muts

def muts_within_coords(ac, mutations, start, end):
    l = []
    for mut_pos in sorted(mutations):
        if mut_pos >= start and mut_pos <= end:
            for mut in mutations[mut_pos]:
                if mut not in l:
                    l.append(mut)
    return l

def get_PPI_from_MongoDB(data, acc_a, acc_b):
    evidence = defaultdict(set)
    for doc in data.find( {"$or":
                                [{"Acc_A": acc_a, "Acc_B": acc_b},
                                 {"Acc_A": acc_b, "Acc_B": acc_a}]},
                             {"_id": 0,
                              "Source:ID:PubMedID:Throughput": 1}):

        for info in doc["Source:ID:PubMedID:Throughput"].split(";"):
            int_id, throughput = info.split(":")[1], info.split(":")[-1]
            for th in ["High", "Low"]:
                if th in throughput:
                    evidence[th].add(int_id)

    return evidence

def extract_DDI_from_MongoDB(data, all_pfams):
    d = {}
    collection = data.find()
    for c in collection:
        pfam_a = c["Pfam_Acc_A"]
        pfam_b = c["Pfam_Acc_B"]
        source = c["Source"]
        pdbs   = [x for x in str(c["PDBs"]).split(";") if x!=""]
        if pfam_a in all_pfams and pfam_b in all_pfams:
            d[(pfam_a, pfam_b)] = {"dbs":source, "pdbs": pdbs}
            d[(pfam_b, pfam_a)] = {"dbs":source, "pdbs": pdbs}
    return d

def get_3did_dmi_from_MongoDB(data):
    d = {}
    for doc in data.find():
        d[(doc["MOTIF"], doc["DOMAIN"])] = (doc["REGEX"], doc["PDBS"])
    return d

def get_pair_association_from_MongoDB(data, dic, acc_a, acc_b,
                                     n_min, obs_min, lo_min):

    p_value, lo = dic.get(acc_a, {}).get(acc_b, (False, False))

    if not (p_value and lo):
        p_value, log_odds = 0, 0
        if data:
            doc = data.find_one({"$or": [{"Acc_A": acc_a, "Acc_B": acc_b},
                                        {"Acc_A": acc_b, "Acc_B": acc_a}]},
                                { "_id": 0 })
            if doc:
                p_value = doc["P-val"]
                if (doc["n_A"] >= n_min
                and doc["n_B"] >= n_min
                and doc["Obs"] >= obs_min
                and doc["LO"] >= lo_min):
                    log_odds = doc["LO"]

        dic[acc_a][acc_b] = (p_value, lo)
        dic[acc_b][acc_a] = (p_value, lo)

    return p_value, log_odds, dic

def string2list_fix(string):

    return [str(x) for x in str(string).upper().replace(" ","").split(",") if x!=""]

def get_elm_dom_from_MongoDB(data):
    d = defaultdict(list)
    for doc in data.find():
        elm = doc["ELM identifier"].upper()
        doms = string2list_fix(doc["Interaction Domain Id"].upper())
        for dom in doms:
            d[(elm, dom)].append( {
                "in_taxon"        : string2list_fix(doc["Present in Taxon"]),
                "not_in_taxon"    : string2list_fix(doc["Not Present in Taxon"]),
                "elm_genes"       : string2list_fix(doc["ELM-containing Genes"]),
                "human_dom_genes" : string2list_fix(doc["Human Domain-containing Genes"]),
                "dom_genes"       : string2list_fix(doc["Domain-containing Genes"]),
                "phos"            : string2list_fix(doc["Phosphosites"]),
                "obs"             : doc["Observations"],
                "req_elms"        : string2list_fix(doc["Requires other ELM"])
            })

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
        p_val = "-"
        z = "-"
        if "Z" in d:
            z = d["Z"]
        if "p-value" in d:
            p_val = "{:1.0e}".format(float(d["p-value"]))

        if (d["#Gene1"] in [gene_a, ac_a]):
            hit = (d["PDB1"], d["Blast-E1"], d["Blast-PCID1"],
                   d["qstart1"], d["qend1"], d["pdbstart1"], d["pdbend1"],
                   d["PDB2"], d["Blast-E2"], d["Blast-PCID2"],
                   d["qstart2"], d["qend2"], d["pdbstart2"], d["pdbend2"],
                   z, p_val)
        elif (d["Gene2"] in [gene_a, ac_a]):
            hit = (d["PDB2"], d["Blast-E2"], d["Blast-PCID2"],
                   d["qstart2"], d["qend2"], d["pdbstart2"], d["pdbend2"],
                   d["PDB1"], d["Blast-E1"], d["Blast-PCID1"],
                   d["qstart1"], d["qend1"], d["pdbstart1"], d["pdbend1"],
                   z, p_val)
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
            # if role == "domain":
            #     node["data"]["color-bg"] = "white magenta"# + color_map[label]


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

def check_gene_match(target_gene, gene_list, sps_tax):
    target_gene = target_gene.upper()
    allow = False

    gene_list2 = []
    for gene in gene_list:
        if "[" in gene:
            gene_tax = re.search("\[(\w+)\]", gene).group(1)
            gene = gene.split("[")[0]
            if gene_tax not in sps_tax:
                continue
        gene_list2.append(gene)

    if len(gene_list2) == 0:
        allow = True
    else:
        for gene in gene_list2:
            if target_gene == gene:
                allow = True
                break
            else:
                if gene == "[transmembrane]":
                    continue
                # if "[" in gene:
                #     gene_tax = re.search("\[(\w+)\]", gene).group(1)
                #     gene = gene.split("[")[0]
                #     if gene_tax not in sps_tax:
                #         continue
                if "?" in gene:
                    if target_gene.startswith(gene.split("?")[0]):
                        allow = True
                        break
                if gene[0]=="^":
                    if target_gene == gene.split("^")[1]:
                        allow = False
                        break
    return allow

def check_ELM_domain_interaction(elm_dom_info, elm_ide, elm_hits, dom_acc,
                                 elm_gene, dom_gene, sps, sp_map, sp_tax, prot_phos, input_seqs):

    elm_ide = elm_ide.upper()
    dom_acc = dom_acc.upper()
    phospho_hits = "none"
    # 'elm_hits' is returned but nothing is done to it. This function should returns only the hits
    # of this motif that pass the filter. Right now all of them do.
    try:
        int_infos = elm_dom_info[(elm_ide, dom_acc)]
    except:
        return False, phospho_hits, elm_hits

    if len(sps)==0: # this happens when both proteins are "input sequences"
        return True, phospho_hits, elm_hits

    sps_tax = []
    for sp in sps:
        sps_tax += sp_tax.get(sp, [])
    sps_tax = set(sps_tax)

    # 1. Check allowed/not allowed taxons
    int_info = "-"
    for info in int_infos:
        if len(sps_tax.intersection(set(info["not_in_taxon"]))) > 0:
            continue
        elif (len(info["in_taxon"]) > 0
        and len(sps_tax.intersection(set(info["in_taxon"]))) == 0):
            continue
        else:
            int_info = info
            break
    if int_info == "-":
        return False, phospho_hits, elm_hits

    # 2. Check genes with ELM
    if len(sps) == 1 and elm_gene not in input_seqs:
        allow = check_gene_match(elm_gene, int_info["elm_genes"], sps_tax)
        if not allow:     
            return False, phospho_hits, elm_hits
    # 3. Check genes with Domain
    if len(sps) == 1 and dom_gene not in input_seqs:
        if "Hsa" in sps:
            allow = check_gene_match(dom_gene, int_info["human_dom_genes"], sps_tax)
        else:
            allow = check_gene_match(dom_gene, int_info["dom_genes"], sps_tax)
        if not allow:     
            return False, phospho_hits, elm_hits

    # 4. Supporting ELMs


    # 5. Check Phosphosites
    if len(int_info["phos"]) > 0:
        elm_phos = [int(i) for i in int_info["phos"]]
        # print "\tPHOS:", elm_phos, prot_phos
        phospho_hits = []
        for hit in elm_hits:
            hit_start = int(hit["start"])
            hit_end   = int(hit["end"])
            hit_phos  = []
            for phos in elm_phos:
                if phos < 0:
                    phos_pos = hit_end + phos + 1
                else:
                    phos_pos = hit_start + phos - 1
                hit_phos.append(phos_pos)
            if len(set(hit_phos).intersection(set(prot_phos))) == len(hit_phos):
                phospho_hits.append((hit_start, hit_end))
            hit["phos_pos"] = hit_phos
        # if len(phospho_hits) == 0:
        #     return False, phospho_hits, elm_hits

    return True, phospho_hits, elm_hits


@line_profile
def main(target_prots, custom_pairs, input_seqs, mutations, sp_map,
        fasta_data, fasta_link, prot_dict, CLIENT, protein_data, cosmic_data,
        fasta_iprets, ddi_data, elm_int_data, dmi_data,
        pfam_info, elm_info,
        make_network=True, hide_no_int=True, inferred_elmdom_ints=False):

    nw = make_network
    nodes, edges = [], []
    id_n = 0
    id_dict, id_coords  = {}, {}
    id_muts = defaultdict(list)
    genes, uni_ids, pfams, phospho_positons = {}, {}, {}, {}
    uni_feat_ints, uni_feat_info = {}, {}
    elms = defaultdict(lambda: defaultdict(list))
    lm3d = defaultdict(lambda: defaultdict(list))
    all_pfams = set()
    ini_pos = get_layout_positions(target_prots | set(input_seqs.keys()))
    p_center = {}

    ### Create Nodes for Input Proteins
    for uni_ac in sorted(list(target_prots)):
        org = sp_map[uni_ac]
        data = protein_data.find_one( { "organism": org, "uni_ac": uni_ac },
                            { "_id": 0, "data_class": 0, "sequence": 0,
                              "alt_ids": 0, "biogrid_interactors": 0})
        genes[uni_ac] = data["genes"][0]
        uni_ids[uni_ac] = data["uni_id"]
        prot_center = float(data["length"])/2
        p_center[uni_ac] = prot_center
        nodes, id_n, id_dict, prot_id = add_protein_main(nodes, id_n, id_dict,
                                                uni_ac, data, ini_pos[uni_ac])

        (nodes, id_n, id_dict, id_coords, id_muts,
        pfams[uni_ac]) = add_domains(nodes, id_n, id_dict, id_coords, id_muts,
                                     prot_id, uni_ac, prot_center,
                                     ini_pos[uni_ac], data["pfams"], pfam_info,
                                     mutations[uni_ac])

        nodes, id_n, phos = add_ptms(nodes, id_n, prot_id, uni_ac, prot_center,
                                     ini_pos[uni_ac], data)

        (nodes, id_n, uni_feat_ints[uni_ac],
            uni_feat_info[uni_ac]) = add_uni_features(nodes, id_n,
                                       prot_id, prot_center, ini_pos[uni_ac],
                                       data["uni_features"], data["regions"],
                                       prot_dict[org])

        nodes, id_n = add_custom_mutations(nodes, id_n, prot_id, uni_ac,
                                           data["length"],
                                           prot_center, ini_pos[uni_ac],
                                           mutations[uni_ac])

        if sp_map[uni_ac] == "Hsa":
            nodes, id_n = add_cosmic_mutations(nodes, id_n, prot_id, uni_ac,
                                               prot_center, ini_pos[uni_ac],
                                               cosmic_data)

        all_pfams.update(pfams[uni_ac])
        for lm in data["elms"]:
            elms[uni_ac][lm["acc"]].append(lm)
        for lm in data["3dlms"]:
            lm3d[uni_ac][lm["acc"]].append(lm)
            
        phospho_positons[uni_ac] = phos

    ### Create Nodes for Input Protein Sequences
    for header, seq in input_seqs.iteritems():
        data = fasta_data[header]
        link = fasta_link[header]
        prot_center = float(len(seq))/2
        p_center[header] = prot_center
        nodes, id_n, id_dict, prot_id = add_fasta_main(nodes, id_n, id_dict,
                                                       header, seq, link,
                                                       data["blast"],
                                                       ini_pos[header])

        (nodes, id_n, id_dict, id_coords, id_muts,
        pfams[header]) = add_domains(nodes, id_n, id_dict, id_coords, id_muts,
                                     prot_id, header, prot_center,
                                     ini_pos[header], data["pfams"], pfam_info,
                                     mutations[header])

        nodes, id_n = add_custom_mutations(nodes, id_n, prot_id, header,
                                           len(seq),
                                           prot_center, ini_pos[header],
                                           mutations[header])

        genes[header] = header
        uni_ids[header] = header
        all_pfams.update(pfams[header])
        elms[header] = data["elms"]
        phospho_positons[header] = []

    ### Load data from Mongo
    ddi = extract_DDI_from_MongoDB(ddi_data, all_pfams)
    elm_dom = get_elm_dom_from_MongoDB(elm_int_data)
    
    dmi = get_3did_dmi_from_MongoDB(dmi_data)

    sp_tax = {
        "Ath" : ["Eukaryota", "Plantae", "Viridiplantae", "Ath"],
        "Cel" : ["Eukaryota", "Opisthokonta", "Metazoa", "Bilateria", "Chordata", "Nematoda", "Cel"],
        "Dme" : ["Eukaryota", "Opisthokonta", "Metazoa", "Bilateria", "Chordata", "Arthropoda", "Insecta", "Dme"],
        "Dre" : ["Eukaryota", "Opisthokonta", "Metazoa", "Bilateria", "Chordata", "Vertebrata", "Dre"],
        "Hsa" : ["Eukaryota", "Opisthokonta", "Metazoa", "Bilateria", "Chordata", "Vertebrata", "Mammalia", "Hsa"],
        "Mmu" : ["Eukaryota", "Opisthokonta", "Metazoa", "Bilateria", "Chordata", "Vertebrata", "Mammalia", "Mmu"],
        "Sce" : ["Eukaryota", "Opisthokonta", "Fungi", "Sce"],
        "Xla" : ["Eukaryota", "Opisthokonta", "Metazoa", "Bilateria", "Chordata", "Vertebrata", "Amphibia", "Xla"],
        "SARS-CoV-2" : ["Eukaryota", "Opisthokonta", "Metazoa", "Chordata", "Vertebrata", "Mammalia", "Hsa","Viridae", "Virus"]
    }
    for s in sp_tax:
        sp_tax[s] = [x.upper() for x in sp_tax[s]]

    ## FIX: IF NUMBER OF PROTEINS ABOVE A THRESHOLD THEN LOAD INFO INTO MEMORY
    ass_dict = defaultdict(dict)
    elm_nodes = defaultdict(list)
    lines = []
    n_ints = defaultdict(int)
    mech_ints = defaultdict(int)
    target_prots = target_prots | set(input_seqs.keys())
    columns = ["UniProt ID (A)","Gene (A)", "UniProt AC (A)", "UniProt ID (B)", "Gene (B)", "UniProt AC (B)",
        "Type", "Element (A)", "Start-End (A)", "Mutations (A)",
        "Element (B)", "Start-End (B)", "Mutations (B)",
        "Scores", "Data Source"]

    if len(custom_pairs) > 0:
        protein_pairs = custom_pairs
    else:
        protein_pairs = itertools.combinations(target_prots, 2)

    for (ac_a, ac_b) in protein_pairs:
        gene_a = genes[ac_a]
        gene_b = genes[ac_b]
        sps = list(set([sp_map[ac] for ac in [ac_a, ac_b] if ac in sp_map]))
        if len(sps) == 1:
            ppi_data = CLIENT[sps[0]]["ppi_db"]
            ass_prob_data = CLIENT[sps[0]]["association_probabilities"]
            iprets_data = CLIENT[sps[0]]["interprets"]
        else:
            ppi_data = None
            ass_prob_data = None
            iprets_data = None

        ### Protein-Protein interaction
        biogrid_ints = []
        if ppi_data:
            evidence = get_PPI_from_MongoDB(ppi_data, ac_a, ac_b)
            biogrid_ints = list(evidence["Low"]) + list(evidence["High"])

        if len(biogrid_ints) > 0:
            id_n += 1
            edges.append({
                "group": "edges",
                "data": {
                    "id": id_n,
                    "source": id_dict[ac_a]["main"],
                    "target": id_dict[ac_b]["main"],
                    "role": "prot_prot_interaction",
                    "ds": "BioGRID",
                    "evidence_n": len(biogrid_ints),
                    "low": list(evidence["Low"]),
                    "high": list(evidence["High"])
                }
            })

            line = [uni_ids[ac_a], gene_a[:20], ac_a, uni_ids[ac_b] ,gene_b[:20], ac_b, "PROT::PROT", "", "",
                    "", "", "", "", len(biogrid_ints), "BioGRID"]
            lines.append(line)
            n_ints[ac_a] += 1
            n_ints[ac_b] += 1

        ### UniProt Features interactions
        for (ac1, ac2) in itertools.permutations([ac_a, ac_b], 2):

            for reg_id in uni_feat_ints.get(ac1, []):
                if ac2 in uni_feat_ints[ac1][reg_id]:
                    role, start, end, info = uni_feat_info[ac1][reg_id]
                    id_n += 1
                    edges.append({
                        "group": "edges",
                        "data": {
                            "role": role+"_interaction",
                            "id": id_n,
                            "source": reg_id,
                            "target": id_dict[ac2]["main"],
                            "label": info,
                            "ds": "UniProt"
                        }
                    })

                    line = [uni_ids[ac1], genes[ac1][:20], ac1, uni_ids[ac2], genes[ac2][:20], ac2,
                            "UniProt region",
                            "Region "+start+"-"+end, start+"-"+end, "",
                            "", "", "", info, "UniProt"]
                    lines.append(line)
                    n_ints[ac_a] += 1
                    n_ints[ac_b] += 1


        ### Domain-Domain Interactions
        for (pfam_a, pfam_b) in itertools.product(pfams[ac_a], pfams[ac_b]):

            ## Probabilities & association
            p_value, lo, ass_dict = get_pair_association_from_MongoDB(
                                        ass_prob_data, ass_dict, pfam_a, pfam_b,
                                        n_min=4, obs_min=5, lo_min=2.0)

            ## From DDI database
            ddi_info = ddi.get((pfam_a, pfam_b), False)
            if ddi_info:
                for i, source in enumerate(id_dict[ac_a][pfam_a], 1):
                    for j, target in enumerate(id_dict[ac_b][pfam_b], 1):
                        if source == target:
                            continue
                        id_n += 1
                        edges.append({
                            "group": "edges",
                            "data": {
                                "id": id_n,
                                "source": source,
                                "target": target,
                                "role": "DOM_interaction",
                                "ds": ddi_info["dbs"].replace(";", ", "),
                                "pdb_n": len(ddi_info["pdbs"]),
                                # "links": pdbs, # uncomment if neccessary
                                "p_val": p_value
                            }
                        })

                        domA = pfam_a+":"+pfam_info[pfam_a]["ide"]
                        domB = pfam_b+":"+pfam_info[pfam_b]["ide"]
                        if len(id_dict[ac_a][pfam_a]) > 1:
                            domA += " ("+str(i)+")"
                        if len(id_dict[ac_b][pfam_b]) > 1:
                            domB += " ("+str(j)+")"
                        line = [uni_ids[ac_a], gene_a[:20], ac_a, uni_ids[ac_b], gene_b[:20], ac_b,
                                "DOM::DOM",
                                domA, "-".join(id_coords[source]),
                                "; ".join(id_muts[source]),
                				domB, "-".join(id_coords[target]),
                                "; ".join(id_muts[target]),
                                "; ".join(sorted(ddi_info["pdbs"])), ddi_info["dbs"]]
                        lines.append(line)
                        n_ints[ac_a]+=1
                        n_ints[ac_b]+=1
                        mech_ints[ac_a]+=1
                        mech_ints[ac_b]+=1

            # Inferred DDI
            if lo > 0:
                for i, source in enumerate(id_dict[ac_a][pfam_a], 1):
                    for j, target in enumerate(id_dict[ac_b][pfam_b], 1):
                        if source == target:
                            continue
                        id_n += 1
                        edges.append({
                            "group": "edges",
                            "data": {
                                "id": id_n,
                                "source": source,
                                "target": target,
                                "role": "iDOM_interaction",
                                "ds": "Predicted",
                                "lo": lo,
                                "p_val": p_value
                            }
                        })

                        domA = pfam_a+":"+pfam_info[pfam_a]["ide"]
                        domB = pfam_b+":"+pfam_info[pfam_b]["ide"]
                        if len(id_dict[ac_a][pfam_a]) > 1:
                            domA += " ("+str(i)+")"
                        if len(id_dict[ac_b][pfam_b]) > 1:
                            domB += " ("+str(j)+")"
                        line = [uni_ids[ac_a], gene_a[:20], ac_a, uni_ids[ac_b], gene_b[:20], ac_b,
                                "iDOM::iDOM",
                                domA, "-".join(id_coords[source]),
                                "; ".join(id_muts[source]),
                                domB, "-".join(id_coords[target]),
                                "; ".join(id_muts[target]),
                                str(lo), "Association Method"]
                        lines.append(line)
                        n_ints[ac_a]+=1
                        n_ints[ac_b]+=1
                        mech_ints[ac_a]+=1
                        mech_ints[ac_b]+=1

        ### ELM-domain Interactions
        for ac1, gene1, ac2, gene2 in zip([ac_a, ac_b], [gene_a, gene_b],
                                          [ac_b, ac_a], [gene_b, gene_a]):
            for elm_acc in elms[ac1]:
                elm_ide = elm_info[elm_acc]["ide"]
                for pfam_acc in pfams[ac2]:

                    ## Annotated ELM-Domain Interaction
                    (cont, phospho_hits,
                    elms[ac1][elm_acc]) = check_ELM_domain_interaction(
                                                elm_dom, elm_ide,
                                                elms[ac1][elm_acc], pfam_acc,
                                                gene1, gene2, sps, sp_map, sp_tax,
                                                phospho_positons[ac1], set(input_seqs.keys()))
                    if not cont and not inferred_elmdom_ints:
                        continue

                    ## Probabilities & association
                    p_value, lo, ass_dict = get_pair_association_from_MongoDB(
                                        ass_prob_data, ass_dict, elm_acc, pfam_acc,
                                        n_min=4, obs_min=5, lo_min=2.0)

                    if cont:
                        ## Add ELM nodes if they don't exist
                        if elm_ide not in id_dict[ac1]:
                            (nodes, id_n, id_dict, id_coords,
                            id_muts) = add_elm_node(nodes, id_n, id_dict,
                                                    id_coords, id_muts,
                                                    ac1, p_center[ac1],
                                                    ini_pos[ac1], elm_acc,
                                                    elms, elm_info,
                                                    mutations, phospho_hits)

                        ## Add interaction edge
                        for i, source in enumerate(id_dict[ac1][elm_ide], 1):
                            for j, target in enumerate(id_dict[ac2][pfam_acc], 1):
                                id_n += 1
                                edges.append({
                                    "group": "edges",
                                    "data": {
                                        "id": id_n,
                                        "source": source,
                                        "target": target,
                                        "role": "ELM_interaction",
                                        "ds": "ELM",
                                        "p_val": p_value
                                    }
                                })

                                elmA = elm_acc+":"+elm_ide
                                domB = pfam_acc+":"+pfam_info[pfam_acc]["ide"]
                                if len(id_dict[ac1][elm_ide]) > 1:
                                    elmA += " ("+str(i)+")"
                                if len(id_dict[ac2][pfam_acc]) > 1:
                                    domB += " ("+str(j)+")"
                                line = [uni_ids[ac1], gene1[:20], ac1, uni_ids[ac2], gene2[:20], ac2,
                                        "ELM::DOM",
                                         elmA, "-".join(id_coords[source]),
                                         "; ".join(id_muts[source]),
                                         domB, "-".join(id_coords[target]),
                                         "; ".join(id_muts[target]),
                                         "", "ELM"]
                                lines.append(line)
                                n_ints[ac1] += 1
                                n_ints[ac2] += 1
                                mech_ints[ac1] += 1
                                mech_ints[ac2] += 1

                    ## Predicted ELM-DOM Interactions
                    if inferred_elmdom_ints and lo > 0:

                        ## Add ELM nodes if they don't exist
                        if elm_ide not in id_dict[ac1]:
                            (nodes, id_n, id_dict, id_coords,
                                id_muts) = add_elm_node(ac1, elm_acc, elms,
                                                        elm_info, id_n, nodes,
                                                        start_pos, id_dict,
                                                        id_coords, id_muts,
                                                        mutations, [])

                        ## Add interaction edge
                        for i, source in enumerate(id_dict[ac1][elm_ide], 1):
                            for j, target in enumerate(id_dict[ac2][pfam_acc], 1):
                                id_n += 1
                                edges.append(
                                        { "group" : "edges",
                                          "data" :
                                            { "id": id_n,
                                              "source": source,
                                              "target": target,
                                              "role": "iELM_interaction",
                                              "ds": "Predicted",
                                              "lo": lo,
                                              "p_val": p_value
                                            }
                                        })

                                elmA = elm_acc+":"+elm_ide
                                domB = pfam_acc+":"+pfam_info[pfam_acc]["ide"]
                                if len(id_dict[ac1][elm_ide]) > 1:
                                    elmA += " ("+str(i)+")"
                                if len(id_dict[ac2][pfam_acc]) > 1:
                                    domB += " ("+str(j)+")"
                                line = [uni_ids[ac1], gene1[:20], ac1, uni_ids[ac2], gene2[:20], ac2, "iELM::DOM",
                                         elmA, "-".join(id_coords[source]),
                                         "; ".join(id_muts[source]),
                                         domB, "-".join(id_coords[target]),
                                         "; ".join(id_muts[target]),
                                         "", "ELM"]
                                lines.append(line)
                                n_ints[ac1]+=1
                                n_ints[ac2]+=1
                                mech_ints[ac1]+=1
                                mech_ints[ac2]+=1

        ## 3DID DMI
        for ac1, gene1, ac2, gene2 in zip([ac_a, ac_b], [gene_a, gene_b],
                                          [ac_b, ac_a], [gene_b, gene_a]):
            for lm_acc in lm3d[ac1]:  
                # print ac1, lm_acc
                for pfam_acc in pfams[ac2]:
                    dmi_info = dmi.get((lm_acc, pfam_acc), False)
                    if dmi_info:
                        ## Probabilities & association
                        p_value, lo, ass_dict = get_pair_association_from_MongoDB(
                                            ass_prob_data, ass_dict, lm_acc, pfam_acc,
                                            n_min=4, obs_min=5, lo_min=2.0)

                        ## Add LM nodes if they don't exist
                        if lm_acc not in id_dict[ac1]:
                            (nodes, id_n, id_dict, id_coords,
                            id_muts) = add_lm_node(nodes, id_n, id_dict,
                                                    id_coords, id_muts,
                                                    ac1, p_center[ac1],
                                                    ini_pos[ac1], lm_acc,
                                                    lm3d, mutations)
                        ## Make edges
                        for i, source in enumerate(id_dict[ac1][lm_acc], 1):
                            for j, target in enumerate(id_dict[ac2][pfam_acc], 1):
                                id_n += 1
                                edges.append(
                                    { "group" : "edges",
                                          "data" :
                                            { "id": id_n,
                                              "source": source,
                                              "target": target,
                                              "role": "DMI_interaction",
                                              "ds": "3did",
                                              "lo": lo,
                                              "p_val": p_value
                                            }
                                    })
                           
                                eleA = lm_acc
                                eleB = pfam_acc+":"+pfam_info[pfam_acc]["ide"]
                                if len(id_dict[ac1][lm_acc]) > 1:
                                    eleA += " ("+str(i)+")"
                                if len(id_dict[ac2][pfam_acc]) > 1:
                                    eleB += " ("+str(j)+")"
                                line = [uni_ids[ac1], gene1[:20], ac1, uni_ids[ac2], gene2[:20], ac2, "LM::DOM",
                                            eleA, "-".join(id_coords[source]),
                                            "; ".join(id_muts[source]),
                                            eleB, "-".join(id_coords[target]),
                                            "; ".join(id_muts[target]),
                                            "", "3did"]
                                lines.append(line)
                                # print line
                                n_ints[ac1]+=1
                                n_ints[ac2]+=1
                                mech_ints[ac1]+=1
                                mech_ints[ac2]+=1
                    
        ## InterPreTS interactions
        hits = []
        if ac_a in input_seqs or ac_b in input_seqs:
            if (ac_a, ac_b) in fasta_iprets:
                hits.append( fasta_iprets[(ac_a, ac_b)] )
            elif (ac_b, ac_a) in fasta_iprets:
                hits.append( fasta_iprets[(ac_b, ac_a)] )
        elif iprets_data:
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
            for (ac, label, pdb, pdb_start, pdb_end,
                start, end, e_val, pcid) in zip([ac_a, ac_b],
                                    [label_a, label_b], [pdb_a, pdb_b],
                                    [pdb_start_a, pdb_start_b], [pdb_end_a, pdb_end_b],
                                    [start_a, start_b], [end_a, end_b],
                                    [eval_a, eval_b], [pcid_a, pcid_b]):

                length = int(end) - int(start)
                id_n += 1
                nodes.append({
                    "group": "nodes",
                    "data": {
                        "id": id_n,
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
                    "position": {
                        "x": ini_pos[ac][0]+int(start)+(float(length)/2)-p_center[ac]-0.5,
                        "y": ini_pos[ac][1]
                    }
                })

                id_dict[ac][label] = id_n
                id_muts[id_n] = muts_within_coords(ac, mutations[ac],
                                                   int(start), int(end))

            ## Add interaction edge
            source = id_dict[ac_a][label_a]
            target = id_dict[ac_b][label_b]
            id_n += 1
            edges.append({
                "group": "edges",
                "data": {
                    "id": id_n,
                    "source": source,
                    "target": target,
                    "role": "INT_interaction",
                    "pdb": pdb,
                    "color": color,
                    "z-score": z,
                    "p-value": pvalue
                }
            })

            line = [uni_ids[ac_a], gene_a[:20], ac_a, uni_ids[ac_b], gene_b[:20], ac_b, "InterPreTS",
                   label_a, str(start_a)+"-"+str(end_a), "; ".join(id_muts[source]),
                   label_b, str(start_b)+"-"+str(end_b), "; ".join(id_muts[target]),
                   pvalue, "InterPreTS prediction"]
            lines.append(line)
            n_ints[ac_a]+=1
            n_ints[ac_b]+=1
            mech_ints[ac_a]+=1
            mech_ints[ac_b]+=1
    
    ## Color domain & LMs nodes
    nodes = color_regions(nodes)
    # nodes = color_regions(nodes, palette="custom1")

    ## Do not display proteins without interactions
    no_int_prots = []
    if hide_no_int:
        for node in nodes:
            if node["data"]["role"] == "protein_main":
                if node["data"]["protein"] not in n_ints:
                    node["data"]["display"]="none"
                    no_int_prots.append(node["data"]["label"])

    graph_elements = nodes + edges

    if not make_network:
        graph_elements = []

    return  graph_elements, columns, lines, no_int_prots
