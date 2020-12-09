#!/usr/bin/env python

import sys, re, os
import math, random, copy, pprint, json
import gzip, itertools
from collections import defaultdict
# from flask_debugtoolbar_lineprofilerpanel.profile import line_profile

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

def merge_gene_names(genes):
    """For a list of similar gene names ("GENE1A", "GENE2B", "GENE3C"),
    returns a string following the format: "GENE(1A,2B,3C)"
    """
    if len(genes)==1:
        return genes[0]
    else:
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

def add_protein_main(prot_nodes, seq_nodes, id_n, id_dict, uni_id, data, ini_pos):

    id_n += 1
    prot_id = copy.deepcopy(id_n)
    id_dict[uni_id] = defaultdict(list)
    id_dict[uni_id]["main"] = prot_id
    gene_name = data["genes"][0]
    if len(data["genes"]) > 1:
        gene_name += "(& more)"
    prot_nodes.append({
        "group": "nodes",
        "data": {
            "id": prot_id,
            # "parent": prot_id,
            "role": "protein_main",
            "label": gene_name,
            "uni_ac": data["uni_ac"],
            "biogrid_id": data["biogrid_ids"].split(", ")[0],
            "des": data["description"],
            "length": data["length"],
            "protein": uni_id,
            "display": "element"
        }
    })

    id_n += 1
    seq_nodes.append({
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

    return prot_nodes, seq_nodes, id_n, id_dict, prot_id

def add_domains(nodes, id_n, id_dict, id_coords, id_muts, prot_id, uni_id,
                prot_center, ini_pos, pfams, pfam_info, muts):
    """Add protein's Pfam domains

    Adds domain as a single node whose width is proportional to domain's length
    """

    pfams2 = defaultdict(list)
    for domain in pfams:
        pfams2[domain["e-value"]].append(domain) # in case there are more
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
                    "e_val": str(domain["e-value"]),
                    "parent": prot_id,
                    "protein": uni_id
                },
                "position": {
                    "x": ini_pos[0] + start+(float(length)/2) - prot_center - 0.5,
                    "y": ini_pos[1]
                }
            })

            id_dict[uni_id][pfam_ac].append(id_n)
            id_coords[id_n] = (str(start), str(end))
            id_muts[id_n] = muts_within_coords(uni_id, muts, start, end)

    return nodes, id_n, id_dict, id_coords, id_muts, pfam_set

def add_ptms(nodes, id_n, prot_id, uni_ac, prot_center, ini_pos, data):

    for ptm_type, role, x in zip(["Phosphorylation", "Acetylation"],
                                 ["mod_phos", "mod_acet"],
                                 ["p","ac"]):
        for ptm in data[ptm_type]:
            pos = ptm["pos"]
            res = ptm["res"]

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

    return nodes, id_n

def add_cosmic_mutations(nodes, id_n, prot_id, uni_ac, prot_center, ini_pos,
                         cosmic_muts):
    for pos in cosmic_muts:
        aa_muts, cosmic_ids, cds_muts, count, tot_count = [], [], [], [], 0
        for mut in sorted(cosmic_muts[pos], key=lambda k: k["sample_num"], reverse=True):
            aa_muts.append(mut["aa_mut"])
            cosmic_ids.append("COSM"+str(mut["cosmic_id"]))
            cds_muts.append(mut["cds_mut"])
            count.append(str(mut["sample_num"]))
            tot_count += int(mut["sample_num"])

        h = 10+tot_count
        if h > 110:
            h = 110
        offset = h/2

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
                "tot_count": tot_count,
                "height": str(h)+"px"
            },
            "position": {
                "x": ini_pos[0] + float(pos) - prot_center - 0.5,
                "y": ini_pos[1] - offset
            }
        })
        
    return nodes, id_n

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

def muts_within_coords(ac, mutations, start, end):
    l = []
    for mut_pos in sorted(mutations):
        if mut_pos >= start and mut_pos <= end:
            for mut in mutations[mut_pos]:
                if mut not in l:
                    l.append(mut)
    return l

def get_PPI_sql(cursor, sql_table, bios_a, bios_b):
    data = defaultdict(set)
    for bio_a in bios_a:
        for bio_b in bios_b:
            cursor.execute("SELECT interaction_id, throughput"+
                            " FROM "+sql_table+ 
                            " WHERE (interactor_a=\'"+bio_a+"\' AND interactor_b=\'"+bio_b+"\')"+
                            " OR (interactor_a=\'"+bio_b+"\' AND interactor_b=\'"+bio_a+"\');")
            for row in cursor.fetchall():
                data[row[1]].add(row[0])
    return data

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

def check_ELM_domain_interaction(elm_dom_info, elm_acc, elm_hits, dom_acc,
                                 elm_gene, dom_gene, sps, org_map, org_tax, prot_phos, input_seqs):

    phospho_hits = None
    # 'elm_hits' is returned but nothing is done to it. This function should returns only the hits
    # of this motif that pass the filter. Right now all of them do.
    try:
        int_infos = elm_dom_info[(elm_acc, dom_acc)]
    except:
        return False, phospho_hits, elm_hits
    
    if len(sps)==0: # this happens when both proteins are "input sequences"
        return True, phospho_hits, elm_hits

    sps_tax = []
    for sp in sps:
        sps_tax += org_tax.get(sp, [])
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

def add_elm_node(nodes, id_n, id_dict, id_coords, id_muts, uni_id,
                 prot_center, ini_pos, elm_acc, elm_hits, elm_info, mutations,
                 phospho_hits):
    ## Add ELM nodes if they don't exist
    if elm_acc not in id_dict[uni_id]:

        elm_ide   = elm_info[elm_acc]["ide"]
        elm_name  = elm_info[elm_acc]["name"]
        elm_des   = elm_info[elm_acc]["des"]
        elm_regex = elm_info[elm_acc]["regex"]
        elm_prob  = elm_info[elm_acc]["prob"]

        for hit in elm_hits:
            hit_start   = hit["start"]
            hit_end     = hit["end"]
            hit_seq     = hit["seq"]
            hit_length  = hit_end - hit_start
            phos        = ""
            if phospho_hits:
                if (hit_start, hit_end) in phospho_hits:
                    phos = ["yes"]+hit["phos_pos"]
                else:
                    phos = ["no"]+hit["phos_pos"]

            hit_status = hit["status"]
            if not hit_status:
                hit_status = "Unknown"
   
            id_n += 1
            nodes.append({
                "group": "nodes",
                "data": {
                    "id": id_n,
                    "parent": id_dict[uni_id]["main"],
                    "protein": uni_id,
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

            id_dict[uni_id][elm_acc].append(id_n)
            id_coords[id_n] = (str(hit_start), str(hit_end))
            id_muts[id_n] = muts_within_coords(uni_id, mutations[uni_id],
                                            hit_start, hit_end)

    return nodes, id_n, id_dict, id_coords, id_muts

def add_lm_node(nodes, id_n, id_dict, id_coords, id_muts, uni_id,
                 prot_center, ini_pos, lm_acc, lm_hits, mutations):

    if lm_acc not in id_dict[uni_id]:

        for hit in lm_hits:
            hit_start   = hit["start"]
            hit_end     = hit["end"]
            hit_seq     = hit["seq"]
            hit_length  = hit_end - hit_start
            hit_status = hit["status"]
            if not hit_status:
                hit_status = "Unknown"
            id_n += 1
            nodes.append({
                "group": "nodes",
                "data": {
                    "id": id_n,
                    "parent": id_dict[uni_id]["main"],
                    "protein": uni_id,
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

            id_dict[uni_id][lm_acc].append(id_n)
            id_coords[id_n] = (str(hit_start), str(hit_end))
            id_muts[id_n] = muts_within_coords(uni_id, mutations[uni_id],
                                            hit_start, hit_end)

    return nodes, id_n, id_dict, id_coords, id_muts

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

def color_nodes(nodes, palette=""):
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

def hide_proteins_without_interactions(protein_nodes, n_ints):
    no_int_prots = []
    for node in protein_nodes:
        if node["data"]["protein"] not in n_ints:
            node["data"]["display"] = "none"
            no_int_prots.append(node["data"]["label"])
    return protein_nodes, no_int_prots
####

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

def parse_uni_interactors(info, prot_dict):
    interactors = set()
    for word in info.split():
        word = re.sub("[^\w]","", word)
        if len(word)>1 and word in prot_dict["AC"]:
            uni_ac = prot_dict["AC"][word][0]
            interactors.add(uni_ac)
    return list(interactors)

def add_uni_features(nodes, id_n, prot_id, prot_center, ini_pos,
                     uni_feats):#, uni_regions, prot_dict):

    # feat_ints, feat_info = {}, {}
    # for region in uni_regions:
    #     start, end = region["start"], region["end"]
    #     try:
    #         length = end - start
    #     except:
    #         length = 0
    #         if str(start).isdigit():
    #             end = start
    #         elif str(end).isdigit():
    #             start = end
    #         else:
    #             continue
    #     id_n += 1
    #     nodes.append({
    #         "group": "nodes",
    #         "data": {
    #             "id": id_n,
    #             "parent": prot_id,
    #             "role": "uni_region",
    #             "label": "("+str(region["start"])+"-"+str(region["end"])+") "+region["info"],
    #             "start": start,
    #             "end": end,
    #             "length": length
    #         },
    #         "position": {
    #             "x": ini_pos[0] + start+(float(length)/2) - prot_center - 0.5,
    #             "y": ini_pos[1]
    #         }
    #     })
    #     interactors = parse_uni_interactors(region["info"], prot_dict)
    #     if len(interactors)>0:
    #         feat_ints[id_n] = interactors
    #         feat_info[id_n] = ("uni_region", str(start), str(end), region["info"])

    uni_roles = {"VARIANT": "uni_var",
                 "MUTAGEN": "uni_mtg",
                 "METAL": "uni_metal",
                 "BINDING": "uni_binding"}
    for feat_type in uni_feats:
        role = uni_roles[feat_type]
        for (start, end) in uni_feats[feat_type]:
            var, info = [], []
            for feature in uni_feats[feat_type][(start, end)]:
                var.append(feature["var"])
                info.append(feature["info"])
            
            length = end - start
            if length <= 1:
                width = 2
                label = str(end)+" "+"; ".join(var)
            else:
                width = length
                label = str(start)+"-"+str(end)+" "+"; ".join(var)
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
                    "length": width,
                    "var": var,
                    "info": info
                },
                "position": {
                    "x": ini_pos[0] + start+(float(length)/2) - prot_center - 0.5,
                    "y": ini_pos[1]
                }
            })
            # interactors = parse_uni_interactors(info, prot_dict)
            # if len(interactors)>0:
            #     feat_ints[id_n] = interactors
            #     feat_info[id_n] = (role, str(start), str(end), feat["info"])

    return nodes, id_n#, feat_ints, feat_info



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


# @line_profile
def main(target_prots, protein_pairs, input_seqs, mutations, org_map,
        fasta_data, fasta_link, conn, protein_data, pfam_matches, 
        lms, ptms, uni_feats, cosmic_muts, fasta_iprets,
        ppi_table, ddi, dmi_elm, dmi_3did, pfam_info, elm_info,
        make_network=True, hide_no_int=True, inferred_elmdom_ints=False):
    nw = make_network
    id_n = 0 ##### MAKE IT GLOBAL
    id_dict, id_coords  = {}, {}
    id_muts = defaultdict(list)
    genes, uni_acs, pfams, phos_positons = {}, {}, {}, {}
    uni_feat_ints, uni_feat_info = {}, {}

    org_tax = {
        "ARATH" : ["Eukaryota", "Plantae", "Viridiplantae", "ARATH"],
        "CAEEL" : ["Eukaryota", "Opisthokonta", "Metazoa", "Bilateria", "Chordata", "Nematoda", "CAEEL"],
        "DROME" : ["Eukaryota", "Opisthokonta", "Metazoa", "Bilateria", "Chordata", "Arthropoda", "Insecta", "DROME"],
        "DANRE" : ["Eukaryota", "Opisthokonta", "Metazoa", "Bilateria", "Chordata", "Vertebrata", "DANRE"],
        "HUMAN" : ["Eukaryota", "Opisthokonta", "Metazoa", "Bilateria", "Chordata", "Vertebrata", "Mammalia", "HUMAN"],
        "MOUSE" : ["Eukaryota", "Opisthokonta", "Metazoa", "Bilateria", "Chordata", "Vertebrata", "Mammalia", "MOUSE"],
        "YEAST" : ["Eukaryota", "Opisthokonta", "Fungi", "YEAST"],
        "XENLA" : ["Eukaryota", "Opisthokonta", "Metazoa", "Bilateria", "Chordata", "Vertebrata", "Amphibia", "XENLA"],
        "SARS-CoV-2" : ["Eukaryota", "Opisthokonta", "Metazoa", "Chordata", "Vertebrata", "Mammalia", "Hsa","Viridae", "Virus"]
    }
    for org in org_tax:
        org_tax[org] = [x.upper() for x in org_tax[org]]

    ini_pos = get_layout_positions(target_prots | set(input_seqs.keys()))
    p_center = {}
    protein_nodes, seq_nodes, region_nodes, position_nodes = [], [], [], []
    edges = []

    cursor = conn.cursor()
    
    ### Create Nodes for Input Proteins
    for uni_id in sorted(list(target_prots)):
        prot_center = float(protein_data[uni_id]["length"])/2
        p_center[uni_id] = prot_center

        (protein_nodes, seq_nodes, id_n, id_dict, 
            prot_id) = add_protein_main(protein_nodes, seq_nodes, id_n, id_dict,
                                        uni_id, protein_data[uni_id], ini_pos[uni_id])
        
        (region_nodes, id_n, id_dict, id_coords, id_muts,
            pfams[uni_id]) = add_domains(region_nodes, id_n, id_dict, 
                                        id_coords, id_muts,
                                        prot_id, uni_id, prot_center,
                                        ini_pos[uni_id], pfam_matches[uni_id], pfam_info,
                                        mutations[uni_id])

        position_nodes, id_n = add_ptms(position_nodes, id_n, prot_id, uni_id, prot_center,
                                     ini_pos[uni_id], ptms[uni_id])
        
        nodes, id_n = add_uni_features(position_nodes, id_n, prot_id, prot_center, ini_pos[uni_id],
                                       uni_feats[uni_id])

        # (nodes, id_n, uni_feat_ints[uni_ac],
        #     uni_feat_info[uni_ac]) = add_uni_features(nodes, id_n,
        #                                prot_id, prot_center, ini_pos[uni_ac],
        #                                data["uni_features"], data["regions"],
        #                                prot_dict[org])

        position_nodes, id_n = add_custom_mutations(
                                    position_nodes, id_n, 
                                    prot_id, uni_id,
                                    protein_data[uni_id]["length"],
                                    prot_center, ini_pos[uni_id],
                                    mutations[uni_id])

        if uni_id in cosmic_muts:
            position_nodes, id_n = add_cosmic_mutations(
                                        position_nodes, id_n, prot_id, uni_id,
                                        prot_center, ini_pos[uni_id], 
                                        cosmic_muts[uni_id])

        phos_positons[uni_id] = [p["pos"] for p in ptms[uni_id]["Phosphorylation"]]
    
  
    ### Create Nodes for Input Protein Sequences
    for header, seq in input_seqs.items():
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

        # genes[header] = header
        # uni_ids[header] = header
        # elms[header] = data["elms"]
        # phos_positons[header] = []



    ## FIX: IF NUMBER OF PROTEINS ABOVE A THRESHOLD THEN LOAD INFO INTO MEMORY
    ass_dict = defaultdict(dict)

    rows = []
    n_ints = defaultdict(int)
    mech_ints = defaultdict(int)
    target_prots = target_prots | set(input_seqs.keys())
    columns = ["UniProt ID (A)","Gene (A)", "UniProt AC (A)", "UniProt ID (B)", "Gene (B)", "UniProt AC (B)",
        "Type", "Element (A)", "Start-End (A)", "Mutations (A)",
        "Element (B)", "Start-End (B)", "Mutations (B)",
        "Scores", "Data Source"]
    
    for (id_a, id_b) in protein_pairs:
        gene_a = merge_gene_names(protein_data[id_a]["genes"])
        gene_b = merge_gene_names(protein_data[id_b]["genes"])
        ac_a = protein_data[id_a]["uni_ac"]
        ac_b = protein_data[id_b]["uni_ac"]
        sps = list(set([org_map[idx] for idx in [id_a, id_b] if idx in org_map]))
   
        # if len(sps) == 1:
        #     ppi_data = CLIENT[sps[0]]["ppi_db"]
        #     ass_prob_data = CLIENT[sps[0]]["association_probabilities"]
        #     iprets_data = CLIENT[sps[0]]["interprets"]
        # else:
        #     ppi_data = None
        #     ass_prob_data = None
        #     iprets_data = None

        row_main = [id_a, gene_a[:20], ac_a, id_b, gene_b[:20], ac_b]

        biogrid_ints = 0

        ### Protein-Protein interaction
        if len(sps) == 1:
            bio_ids_a = protein_data[id_a]["biogrid_ids"].split(", ")
            bio_ids_b = protein_data[id_b]["biogrid_ids"].split(", ")
            evidence = get_PPI_sql(cursor, ppi_table, bio_ids_a, bio_ids_b)
            biogrid_ints = len(list(evidence["Low"]) + list(evidence["High"]))
    
        if biogrid_ints > 0:
            id_n += 1
            edges.append({
                "group": "edges",
                "data": {
                    "id": id_n,
                    "source": id_dict[id_a]["main"],
                    "target": id_dict[id_b]["main"],
                    "role": "prot_prot_interaction",
                    "ds": "BioGRID",
                    "evidence_n": biogrid_ints,
                    "low": list(evidence["Low"]),
                    "high": list(evidence["High"])
                }
            })

            row = ["PROT::PROT", "", "","", "", "", "", biogrid_ints, "BioGRID"]
            rows.append(row_main+row)
            n_ints[id_a] += 1
            n_ints[id_b] += 1

        ### UniProt Features interactions
        # for (ac1, ac2) in itertools.permutations([ac_a, ac_b], 2):

        #     for reg_id in uni_feat_ints.get(ac1, []):
        #         if ac2 in uni_feat_ints[ac1][reg_id]:
        #             role, start, end, info = uni_feat_info[ac1][reg_id]
        #             id_n += 1
        #             edges.append({
        #                 "group": "edges",
        #                 "data": {
        #                     "role": role+"_interaction",
        #                     "id": id_n,
        #                     "source": reg_id,
        #                     "target": id_dict[ac2]["main"],
        #                     "label": info,
        #                     "ds": "UniProt"
        #                 }
        #             })

        #             row = [uni_ids[ac1], genes[ac1][:20], ac1, uni_ids[ac2], genes[ac2][:20], ac2,
        #                     "UniProt region",
        #                     "Region "+start+"-"+end, start+"-"+end, "",
        #                     "", "", "", info, "UniProt"]
        #             rows.append(row)
        #             n_ints[ac_a] += 1
        #             n_ints[ac_b] += 1
 
        ### Domain-Domain Interactions
        for (pfam_a, pfam_b) in itertools.product(sorted(list(pfams[id_a])), sorted(list(pfams[id_b]))):
        
            ## Probabilities & association
            p_value, lo = 0, 0
            # p_value, lo, ass_dict = get_pair_association_from_MongoDB(
            #                             ass_prob_data, ass_dict, pfam_a, pfam_b,
            #                             n_min=4, obs_min=5, lo_min=2.0)

            ## From DDI database
            ddi_info = ddi.get((pfam_a, pfam_b), False)
            if ddi_info:
                for i, source in enumerate(id_dict[id_a][pfam_a], 1):
                    for j, target in enumerate(id_dict[id_b][pfam_b], 1):
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
                                "ds": ddi_info["dbs"],
                                "pdb_n": len(ddi_info["pdbs"]),
                                # "links": pdbs, # uncomment if neccessary
                                "p_val": p_value
                            }
                        })

                        domA = pfam_a+":"+pfam_info[pfam_a]["ide"]
                        domB = pfam_b+":"+pfam_info[pfam_b]["ide"]
                        if len(id_dict[id_a][pfam_a]) > 1:
                            domA += " ("+str(i)+")"
                        if len(id_dict[id_b][pfam_b]) > 1:
                            domB += " ("+str(j)+")"
                        pdbs = [pdb for i, pdb in enumerate(ddi_info["pdbs"]) if i<3]
                        row = ["DOM::DOM",
                                domA, "-".join(id_coords[source]),
                                "; ".join(id_muts[source]),
                				domB, "-".join(id_coords[target]),
                                "; ".join(id_muts[target]),
                                "; ".join(sorted(pdbs)), ddi_info["dbs"]]
                        rows.append(row_main+row)
                        n_ints[id_a]+=1
                        n_ints[id_b]+=1
                        mech_ints[id_a]+=1
                        mech_ints[id_b]+=1
            
            # Inferred DDI
            if lo > 0:
                for i, source in enumerate(id_dict[id_a][pfam_a], 1):
                    for j, target in enumerate(id_dict[id_b][pfam_b], 1):
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
                        if len(id_dict[id_a][pfam_a]) > 1:
                            domA += " ("+str(i)+")"
                        if len(id_dict[id_b][pfam_b]) > 1:
                            domB += " ("+str(j)+")"
                        row = ["iDOM::iDOM",
                                domA, "-".join(id_coords[source]),
                                "; ".join(id_muts[source]),
                                domB, "-".join(id_coords[target]),
                                "; ".join(id_muts[target]),
                                str(lo), "Association Method"]
                        rows.append(row_main+row)
                        n_ints[id_a]+=1
                        n_ints[id_b]+=1
                        mech_ints[id_a]+=1
                        mech_ints[id_b]+=1
        
        ### ELM-domain Interactions
        for id1, gene1, id2, gene2 in zip([id_a, id_b], [gene_a, gene_b],
                                          [id_b, id_a], [gene_b, gene_a]):
            for lm_acc in lms[id1]:
                if lm_acc in elm_info:
                    elm_acc = lm_acc
                    elm_ide = elm_info[lm_acc]["ide"]
                    elm_hits = lms[id1][lm_acc]
                    for pfam_acc in pfams[id2]:
                        ## Annotated ELM-Domain Interaction
                        (reviewed, phospho_hits,
                        elm_hits) = check_ELM_domain_interaction(
                                        dmi_elm, elm_acc, elm_hits, 
                                        pfam_acc, gene1, gene2, sps, org_map, org_tax,
                                        phos_positons[id1], set(input_seqs.keys()))
                        accept_motif = False
                        if reviewed or inferred_elmdom_ints:
                            p_value, lo = 0, 0
                            # p_value, lo, ass_dict = get_pair_association_from_MongoDB(
                            #                     ass_prob_data, ass_dict, elm_acc, pfam_acc,
                            #                     n_min=4, obs_min=5, lo_min=2.0)
                            if reviewed:
                                accept_motif = True
                                edge_role = "ELM_interaction" #"iELM_interaction",
                                edge_ds = "ELM" #Predicted
                                int_type = "ELM::DOM" #"iELM::DOM"
                            elif inferred_elmdom_ints and lo > 0:
                                accept_motif = True
                                edge_role = "iELM_interaction",
                                edge_ds = "Predicted"
                                int_type = "iELM::DOM"

                        if accept_motif:
                            (region_nodes, id_n, id_dict, id_coords,
                            id_muts) = add_elm_node(region_nodes, id_n, id_dict,
                                                    id_coords, id_muts,
                                                    id1, p_center[id1],
                                                    ini_pos[id1], elm_acc,
                                                    elm_hits, elm_info,
                                                    mutations, phospho_hits)
                         
                            ## Add interaction edge
                            for i, source in enumerate(id_dict[id1][elm_acc], 1):
                                for j, target in enumerate(id_dict[id2][pfam_acc], 1):
                                    id_n += 1
                                    edges.append({
                                        "group": "edges",
                                        "data": {
                                            "id": id_n,
                                            "source": source,
                                            "target": target,
                                            "role": edge_role,
                                            "ds": edge_ds,
                                            "p_val": p_value
                                        }
                                    })

                                    elmA = elm_acc+":"+elm_ide
                                    domB = pfam_acc+":"+pfam_info[pfam_acc]["ide"]
                                    if len(id_dict[id1][elm_acc]) > 1:
                                        elmA += " ("+str(i)+")"
                                    if len(id_dict[id2][pfam_acc]) > 1:
                                        domB += " ("+str(j)+")"
                                    row = [ int_type,
                                            elmA, "-".join(id_coords[source]),
                                            "; ".join(id_muts[source]),
                                            domB, "-".join(id_coords[target]),
                                            "; ".join(id_muts[target]),
                                            "", "ELM"]
                                    rows.append(row_main+row)
                                    n_ints[id1] += 1
                                    n_ints[id2] += 1
                                    mech_ints[id1] += 1
                                    mech_ints[id2] += 1
                                    
                else: #3dmi
                    for pfam_acc in pfams[id2]:
                        dmi_info = dmi_3did.get((lm_acc.upper(), pfam_acc), False)
                        lm_hits = lms[id1][lm_acc]
                        if dmi_info:
                            p_value, lo = 0, 0
                            ## Probabilities & association
                            # p_value, lo, ass_dict = get_pair_association_from_MongoDB(
                            #                     ass_prob_data, ass_dict, lm_acc, pfam_acc,
                            #                     n_min=4, obs_min=5, lo_min=2.0)

                            (region_nodes, id_n, id_dict, id_coords,
                            id_muts) = add_lm_node(region_nodes, id_n, id_dict,
                                                    id_coords, id_muts,
                                                    id1, p_center[id1],
                                                    ini_pos[id1], lm_acc,
                                                    lm_hits, mutations)
                   
                            ## Make edges
                            for i, source in enumerate(id_dict[id1][lm_acc], 1):
                                for j, target in enumerate(id_dict[id2][pfam_acc], 1):
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
                                    if len(id_dict[id1][lm_acc]) > 1:
                                        eleA += " ("+str(i)+")"
                                    if len(id_dict[id2][pfam_acc]) > 1:
                                        eleB += " ("+str(j)+")"
                                    row = ["LM::DOM",
                                            eleA, "-".join(id_coords[source]),
                                            "; ".join(id_muts[source]),
                                            eleB, "-".join(id_coords[target]),
                                            "; ".join(id_muts[target]),
                                            "", "3did"]
                                    rows.append(row_main+row)
                                    n_ints[id1]+=1
                                    n_ints[id2]+=1
                                    mech_ints[id1]+=1
                                    mech_ints[id2]+=1

        ## InterPreTS interactions
        # hits = []
        # if ac_a in input_seqs or ac_b in input_seqs:
        #     if (ac_a, ac_b) in fasta_iprets:
        #         hits.append( fasta_iprets[(ac_a, ac_b)] )
        #     elif (ac_b, ac_a) in fasta_iprets:
        #         hits.append( fasta_iprets[(ac_b, ac_a)] )
        # elif iprets_data:
        #     hits = get_Interprets_from_MongoDB(iprets_data, gene_a, ac_a, gene_b, ac_b)

        # for hit in hits:
        #     pdb_a, eval_a, pcid_a, start_a, end_a, pdb_start_a, pdb_end_a = hit[:7]
        #     pdb_b, eval_b, pcid_b, start_b, end_b, pdb_start_b, pdb_end_b = hit[7:-2]
        #     z, pvalue = hit[-2:]

        #     label_a = pdb_a+":"+str(start_a)+"-"+str(end_a)
        #     label_b = pdb_b+":"+str(start_b)+"-"+str(end_b)
        #     eval_avg = (float(eval_a)+float(eval_b)) / 2
        #     eval_diff = abs(float(eval_a)-float(eval_b))
        #     color = color_from_zvalue(z)

        #     ## Add homology region node
        #     for (ac, label, pdb, pdb_start, pdb_end,
        #         start, end, e_val, pcid) in zip([ac_a, ac_b],
        #                             [label_a, label_b], [pdb_a, pdb_b],
        #                             [pdb_start_a, pdb_start_b], [pdb_end_a, pdb_end_b],
        #                             [start_a, start_b], [end_a, end_b],
        #                             [eval_a, eval_b], [pcid_a, pcid_b]):

        #         length = int(end) - int(start)
        #         id_n += 1
        #         nodes.append({
        #             "group": "nodes",
        #             "data": {
        #                 "id": id_n,
        #                 "parent": id_dict[ac]["main"],
        #                 "role": "iprets",
        #                 "label": label,
        #                 "pdb": pdb,
        #                 "pdb_start": int(pdb_start),
        #                 "pdb_end": int(pdb_end),
        #                 "start": int(start),
        #                 "end": int(end),
        #                 "length": length,
        #                 "eval": "{:1.0e}".format(float(e_val)),
        #                 "pcid": pcid,
        #                 "color": color,
        #                 "protein": ac
        #             },
        #             "position": {
        #                 "x": ini_pos[ac][0]+int(start)+(float(length)/2)-p_center[ac]-0.5,
        #                 "y": ini_pos[ac][1]
        #             }
        #         })

        #         id_dict[ac][label] = id_n
        #         id_muts[id_n] = muts_within_coords(ac, mutations[ac],
        #                                            int(start), int(end))

        #     ## Add interaction edge
        #     source = id_dict[ac_a][label_a]
        #     target = id_dict[ac_b][label_b]
        #     id_n += 1
        #     edges.append({
        #         "group": "edges",
        #         "data": {
        #             "id": id_n,
        #             "source": source,
        #             "target": target,
        #             "role": "INT_interaction",
        #             "pdb": pdb,
        #             "color": color,
        #             "z-score": z,
        #             "p-value": pvalue
        #         }
        #     })

        #     row = [uni_ids[ac_a], gene_a[:20], ac_a, uni_ids[ac_b], gene_b[:20], ac_b, "InterPreTS",
        #            label_a, str(start_a)+"-"+str(end_a), "; ".join(id_muts[source]),
        #            label_b, str(start_b)+"-"+str(end_b), "; ".join(id_muts[target]),
        #            pvalue, "InterPreTS prediction"]
        #     rows.append(row)
        #     n_ints[ac_a]+=1
        #     n_ints[ac_b]+=1
        #     mech_ints[ac_a]+=1
        #     mech_ints[ac_b]+=1
    
    ## Do not display proteins without interactions
    no_int_prots = []
    if hide_no_int:
        protein_nodes, no_int_prots = hide_proteins_without_interactions(
                                                    protein_nodes, n_ints)
    graph_elements = []
    if make_network:
        ## Color domain & LMs nodes
        region_nodes = color_nodes(region_nodes) #, palette="custom1")
        graph_elements = protein_nodes+seq_nodes+region_nodes+position_nodes+edges

    return graph_elements, columns, rows, no_int_prots
