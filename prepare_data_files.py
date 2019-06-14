#!/usr/bin/env python

import sys, re, os, gzip, math, pprint, datetime
import find_all_elms
import os.path
from collections import defaultdict
from Bio import SwissProt
import run_interprets
from Bio import SearchIO, SeqIO

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

def get_protein_data_from_uniprot_text(uniprot_file):
    """Extracts information from a UniProt text file.
    Requires Bio::SwissProt module
    """

    D = defaultdict(lambda: defaultdict(list))
    masks = {}
    ref_proteome = set()
    modres = defaultdict(lambda: defaultdict(dict))
    mutagen = defaultdict(lambda: defaultdict(list))
    region = defaultdict(list)
    pdb2uni = defaultdict(set)
    for record in SwissProt.parse(open_file(uniprot_file)):
        dc  = record.data_class
        uni_id = record.entry_name
        accs = record.accessions
        uni_ac = accs[0].upper()
        des = re.search("Name: Full=([^;|{]+)", record.description).group(1)
        seq = record.sequence
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
            for key in [uni_ac, uni_id, gn,
                        uni_ac.upper(), uni_id.upper(), gn.upper()]:
                if val not in D[dic][key]:
                    if dc == "Reviewed":
                        D[dic][key].insert(0, val)
                    else:
                        D[dic][key].append(val)

        for ac in accs[1:]:
            for ac0 in [ac, ac.upper()]:
                if uni_ac not in D["AC"][ac0]:
                    if dc == "Reviewed":
                        D["AC"][ac0].insert(0, uni_ac)
                    else:
                        D["AC"][ac0].append(uni_ac)

        for g in gns[1:]:
            for g0 in [g, g.upper()]:
                if uni_ac not in D["AC"][g0]:
                    if dc == "Reviewed":
                        D["AC"][g0].insert(0, uni_ac)
                    else:
                        D["AC"][g0].append(uni_ac)

        D["des"][uni_ac] = des
        D["dc"][uni_ac]  = dc
        D["seq"][uni_ac] = seq
        masks[uni_ac] = "0" * len(record.sequence)

        if dc == "Reviewed":
            ref_proteome.add(uni_ac)

        # for ref in record.cross_references:
        #     if "PDB" in ref:
        #         chains = set()
        #         for x in ref[-1].split(", "):
        #             for y in x.split("=")[0].split("/"):
        #                 chains.add(y)
        #         for c in chains:
        #             pdb2uni["pdb|"+ref[1]+"|"+c].add(uni_ac)

        for feat in record.features:
            if feat[0]=="MOD_RES":
                pos = feat[1]
                res = seq[pos-1]
                if "Phospho" in feat[3].split(";")[0]:
                    modres[uni_ac][str(pos)+"-"+res]["p"] = ["mod_res"]

                elif "acetyl" in feat[3].split(";")[0]:
                    modres[uni_ac][str(pos)+"-"+res]["ac"] = ["mod_res"]

            elif feat[0]=="MUTAGEN":
                start = str(feat[1])
                end = str(feat[2])
                pos = start+"-"+end
                if start == end:
                    pos = start
                info = feat[3].split(" {ECO")[0]
                mutagen[uni_ac][pos].append(info)

            elif feat[0]=="REGION":
                start = str(feat[1])
                end = str(feat[2])
                info = feat[3].split(" {ECO")[0]
                region[uni_ac].append(
                    {"start": start,
                     "end": end,
                     "info": info,
                    }
                )
    return D, masks, ref_proteome, modres, mutagen, region #,pdb2uni

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

def get_pfam_doms(pfam_file, prot_dict, max_eval=999):
    """Pfam-A matches in species proteome. File downloaded from PFAM.
    """
    pfams = defaultdict(lambda: defaultdict(set) )
    pfams_temp = defaultdict(lambda: defaultdict(set) )
    pfam_names = {}
    with open_file(pfam_file) as f:
        for line in f:
            if line[0] != "#":
                t = line.rstrip().split("\t")
                if len(t)>14:
                    t.remove(t[1])
                uni_ac = t[0].upper()
                start, end = int(t[3]), int(t[4])
                pfam_ac, pfam_name, domain_e_val = t[5], t[6], float(t[12])

                # if domain_e_val <= max_eval: ## No e-value cut-off. Just take what is annotated by Pfam

                if uni_ac in prot_dict["seq"]: ## Keep annotation for primary UniProt accessions
                    pfams[uni_ac][pfam_ac].add((start, end, domain_e_val))
                    pfam_names[pfam_ac]=pfam_name
                elif uni_ac in prot_dict["AC"]: ## Save annotation for secondary accessions in a different dictionary
                    pfams_temp[uni_ac][pfam_ac].add((start, end, domain_e_val))

    ## Transfer annotation to those main accessions which were not annotated
    for alt_ac in pfams_temp:
        for uni_ac in prot_dict["AC"][alt_ac]:
            if uni_ac not in pfams: ## If the associated main accession is not annotated already:
                pfams[uni_ac] = pfams_temp[alt_ac]
                pfam_names[pfam_ac]=pfam_name
    return pfams, pfam_names

def edit_pfam_dat(pfam_dat_file, out_file):
    with open_file(out_file, "w") as out:
        cols = ["Acc", "Ide", "Des"]
        out.write("\t".join(cols)+"\n")

        with open_file(pfam_dat_file) as f:
            for line in f:
                if line.startswith("#=GF ID"):
                    ide = re.search("#=GF\s+ID\s+(.+)", line.rstrip()).group(1)
                elif line.startswith("#=GF AC"):
                    acc = re.search("#=GF\s+AC\s+(.+)", line.rstrip()).group(1)
                    acc = acc.split(".")[0]
                elif line.startswith("#=GF DE"):
                    des = re.search("#=GF\s+DE\s+(.+)", line.rstrip()).group(1)
                elif line.startswith("//"):
                    out.write("\t".join([acc, ide, des])+"\n")

def reduce_3did(db_file, out_file):
    """ Writes a simplified version in TSV format of the 3did flat file db
    """
    with open_file(out_file, "w") as out:
        cols = ["Pfam_Ide_A", "Pfam_Acc_A", "Pfam_Ide_B", "Pfam_Acc_B", "PDBs"]
        out.write("\t".join(cols)+"\n")

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

def compare_ints(t, gene1, gene2, ints):
    s1, e1 = int(t[4]), int(t[5])
    s2, e2 = int(t[12]), int(t[13])
    for info in ints[gene1][gene2]:
        i = info.split("\t")
        sA, eA = int(i[4]), int(i[5])
        sB, eB = int(i[12]), int(i[13])

        mask1 = [0 for _ in range(0,eA)]
        mask2 = [0 for _ in range(0,eB)]
        if e1 > eA:
            mask1 = [0 for _ in range(0,e1)]
        if e2 > eB:
            mask2 = [0 for _ in range(0,e2)]

        mask1 = mask1[:sA-1]+[n+1 for n in mask1[sA-1:eA]]+mask1[eA:]
        mask1 = mask1[:s1-1]+[n+1 for n in mask1[s1-1:e1]]+mask1[e1:]
        overlap1 = mask1.count(2)

        mask2 = mask2[:sB-1]+[n+1 for n in mask2[sB-1:eB]]+mask2[eB:]
        mask2 = mask2[:sB-1]+[n+1 for n in mask2[sB-1:eB]]+mask2[eB:]
        overlap2 = mask2.count(2)


        print "OLD",gene1,":",sA,"-",eA,"\t",gene2,":",sB,"-",eB
        print "NEW",gene1,":",s1,"-",e1,"\t",gene2,":",s2,"-",e2,"\n"
        # if overlap1 == 0 and overlap2 == 0:
        #
        #     sys.exit()

def check_file(somefile):
    if os.path.isfile(somefile):
        print_status(somefile, "exists")
        return True
    else:
        sys.exit("FileError: '"+somefile+"' does not exist.\n"+
            "Make sure the file is located and named correctly\n")

def get_interacting_linear_motifs(elm_dom_file, elm_hits_file, masks,
                                  max_overlap, max_eval=1):
    """ Reads pre-generated file with ELMs annotated for each protein
        Restricted to ELMs appearing in the elm_dom interaction file. PIV won't use any other
    """
    int_elms = []
    with open_file(elm_dom_file) as f:
        for line in f:
            t = line.rstrip().split("\t")
            if line[0]!="#":
                int_elms.append(t[0])

    elms = defaultdict(lambda: defaultdict(set) )
    elms_temp = defaultdict(lambda: defaultdict(set) )
    with open_file(elm_hits_file) as f:
        for line in f:
            if line[0] != "#" and line.strip():
                t = line.rstrip().split("\t")
                elm_name, elm_ac = t[0], t[1]
                uni_ac, uni_id = t[2].split("|")[1:]
                elm_start, elm_end = int(t[4]), int(t[5])
                prob_score = float(t[6]) # probability score based on the combined expected frequencies of the AAs in the regular expression

                if elm_name in int_elms:
                        overlap = calculate_overlap(elm_start, elm_end, masks[uni_ac])
                        if overlap <= max_overlap:
                            masks[uni_ac] = fill_mask(elm_start, elm_end, masks[uni_ac])
                            elms[uni_ac][elm_ac].add((elm_start, elm_end))
    return elms

def make_protein_data_json(prot_dict, pfams, elms,
                           ptms, mutagen, regions,
                           outfile_name, mode="mongo"):

    pp = pprint.PrettyPrinter(indent=4)
    json_data = {}
    with open_file(outfile_name, "w") as out:
        for uni_ac in prot_dict["seq"]:

            seq = prot_dict["seq"][uni_ac]
            protein_data = {
                    "uni_ac": uni_ac,
                    "uni_id": list(prot_dict["ID"][uni_ac])[0], # always, only 1
                    "gene": list(prot_dict["GN"][uni_ac])[0],
                    "description": prot_dict["des"][uni_ac],
                    "data_class": prot_dict["dc"][uni_ac],
                    "length" : len(seq),
                    "sequence" : seq,
                    "pfams" : [],
                    "elms" : [],
                    "phosphorylation" : [],
                    "acetylation" : [],
                    "mutagen" : [],
                    "regions" : []
            }

            if uni_ac in pfams:
                for pfam_ac in sorted(pfams[uni_ac]):
                    for (start, end, e_val) in sorted(pfams[uni_ac][pfam_ac],
                                                      key=lambda x: int(x[0])):
                        protein_data["pfams"].append(
                            {
                                "acc" :   pfam_ac,
                                "start" : int(start),
                                "end" :   int(end),
                                "e-val" : float(e_val)
                            })

            if uni_ac in elms:
                for elm_ac in sorted(elms[uni_ac]):
                    for (start, end) in sorted(elms[uni_ac][elm_ac],
                                                    key=lambda x: int(x[0])):
                        protein_data["elms"].append(
                            {
                                "acc" :   elm_ac,
                                "start" : int(start),
                                "end" :   int(end),
                                "seq" :   seq[int(start)-1:int(end)]
                            })

            if uni_ac in ptms:
                for pos_res in ptms[uni_ac]:
                    pos, res = pos_res.split("-")

                    for mod, source in ptms[uni_ac][pos_res].iteritems():
                        if mod == "p":
                            ptm_type = "phosphorylation"
                        elif mod == "ac":
                            ptm_type = "acetylation"
                        protein_data[ptm_type].append(
                            {
                                "pos" : int(pos),
                                "res" : res,
                                "source": ";".join(source)
                            })
            # else:
            #     for i in range(9):
            #         iso_uni_ac = uni_ac+"-"+str(i)
            #         if iso_uni_ac in ptms:
            #             for (res, pos, mod) in sorted(ptms[iso_uni_ac],
            #                                          key=lambda x: int(x[1])):
            #                 if mod == "p":
            #                     ptm_type = "phosphorylation"
            #                 elif mod == "ac":
            #                     ptm_type = "acetylation"
            #                 protein_data[ptm_type].append(
            #                     {
            #                         "pos" : int(pos),
            #                         "res" : res,
            #                         "source": "PSP,-"+str(i)
            #                     })
            #             break

            if uni_ac in mutagen:
                for pos in mutagen[uni_ac]:
                    protein_data["mutagen"].append(
                        {
                            "pos": pos,
                            "info": mutagen[uni_ac][pos]
                        }
                    )

            if uni_ac in regions:
                protein_data["regions"] = regions[uni_ac]

            if mode == "mongo":
                out.write(str(protein_data)+"\n")
            elif mode == "json":
                json_data[uni_ac] = protein_data

        if mode == "json":
            json.dump(json_data, out)

        print_status(outfile_name, "updated")

def hasNumbers(inputString):
	return any(char.isdigit() for char in inputString)

def assign_Number(inputString):
	if ";" in inputString:
		return inputString.split(";")[0]
	else:
		return inputString

def get_ptms(psp_file, prot_dict, modres):
    psp = defaultdict(set)
    psp_temp = defaultdict(set)
    with open_file(psp_file) as f:
        for line in f:
            t = line.rstrip().split("\t")
            pos, mod = t[4].split("-")
            res, pos = re.search("(\w)([\d]+)", pos).group(1,2)
            lt_lit="."
            ms_lit="."
            ms_cst="."
            cst_cat="."
            uni_ac = t[2]
            if len(t)>10 and hasNumbers(t[10]):
                lt_lit = assign_Number(t[10])
            if len(t)>11 and hasNumbers(t[11]):
                ms_lit = assign_Number(t[11])
            if len(t)>12 and hasNumbers(t[12]):
                ms_cst = assign_Number(t[12])
            if len(t)>13 and hasNumbers(t[13]):
                cst_cat = assign_Number(t[13])

            if uni_ac in prot_dict["seq"]:
                seq = prot_dict["seq"][uni_ac]
                if len(seq)>=int(pos) and seq[int(pos)-1]==res:
                    psp[uni_ac].add( (res, pos, mod) )
                    if (uni_ac in modres and pos+"-"+res in modres[uni_ac]
                    and mod in modres[uni_ac][pos+"-"+res]):
                        modres[uni_ac][pos+"-"+res][mod].append("PSP")
                    else:
                        modres[uni_ac][pos+"-"+res][mod] = ["PSP"]
            elif uni_ac in prot_dict["AC"]:
                psp_temp[uni_ac].add( (res, pos, mod) )

    for alt_ac in psp_temp:
        for uni_ac in prot_dict["AC"][alt_ac]:
            if uni_ac not in psp:
                seq = prot_dict["seq"][uni_ac]
                for (res,pos,mod) in psp_temp[alt_ac]:
                    if len(seq)>=int(pos) and seq[int(pos)-1]==res:
                        psp[uni_ac].add( (res, pos, mod) )
                        if (uni_ac in modres and pos+"-"+res in modres[uni_ac]
                        and mod in modres[uni_ac][pos+"-"+res]):
                            modres[uni_ac][pos+"-"+res][mod].append("PSP")
                        else:
                            modres[uni_ac][pos+"-"+res][mod] = ["PSP"]

    return modres

def print_status(thefile, status):
    l = 62-len(thefile)-len(status)
    print "'"+thefile+"'"+"."*l+status

def main( sp="Hsa",
          data_dir="",
          max_overlap=0.2,
          no_overlap_between_doms_and_lms=True,
          force=False
        ):

    uni_version  = "May2019"
    pfam_version = "Aug2018"
    elm_version  = "May2019"
    psp_version  = "Mar2019"
    biogrid_version = "3.5.172"

    ## Common files:
    com_dir = data_dir+"common/"
    pfam_dat_file    = com_dir+"Pfam-A.hmm_r32.0.dat.gz"
    flat_3did_file   = com_dir+"3did_flat-2018_04.gz"
    # pdbchain2uniprot = com_dir+"pdbsws_chain.txt.gz"
    elm_classes_file = com_dir+"elm_classes_"+elm_version+".tsv.gz"
    elm_int_dom_file = com_dir+"elm_interaction_domains_edited_Jan18.tsv"
    ## Output Files:
    edited_pfam_dat_file = com_dir+"Pfam-A.hmm_r32.0.tsv.gz"
    edited_3did_file     = com_dir+"3did_flat_edited-2018_04.tsv.gz"

    ## Species files:
    sp_dir = data_dir+"species/"+sp+"/"
    uni_text_file     = sp_dir+"uniprot_proteome_"+sp+"_"+uni_version+".txt.gz"
    uni_fasta_file    = sp_dir+"uniprot_proteome_"+sp+"_"+uni_version+".fasta.gz"
    pfam_matches_file = sp_dir+"pfamA_matches_"+sp+"_"+pfam_version+".tsv.gz"
    elm_hits_file     = sp_dir+"elm_hits_"+sp+"_"+elm_version+".tsv.gz"
    biogrid_file      = sp_dir+"BIOGRID-ORGANISM-"+sp+"-"+biogrid_version+".tab2.txt.gz"
    pfam_assoc_file   = sp_dir+"dom_dom_association_"+sp+".tsv.gz"
    prot_data_file    = sp_dir+"protein_data_"+sp+"_mongo.json.gz"
    iprets_file       = sp_dir+"i2_biogrid_"+sp+".txt.gz"
    if sp in ["Hsa", "Mmu"]:
        psp_file      = sp_dir+"PSP_ptms_"+sp+"_"+psp_version+".tsv.gz"
    if sp == "Hsa":
        iprets_file   = sp_dir+"human_aaa_biogrid_i2.txt.gz"
        hippie_file   = sp_dir+"hippie_v2-1_July2018.tsv.gz"
        hippie_uniprot_mapping = sp_dir+"hippie_v2-1_uniprot_mapping_table.tsv.gz"

    ## Check/Edit common files.
    print "Common files:"
    # Pfam.dat
    if not os.path.isfile(edited_pfam_dat_file):
        if check_file(pfam_dat_file):
            edit_pfam_dat(pfam_dat_file, edited_pfam_dat_file)
            print_status(edited_pfam_dat_file, "created")
    else:
        print_status(edited_pfam_dat_file, "exists")
    # 3did flat
    if not os.path.isfile(edited_3did_file):
        if check_file(flat_3did_file):
            reduce_3did(flat_3did_file, edited_3did_file)
            print_status(edited_3did_file, "created")
    else:
        print_status(edited_3did_file, "exists")
    # elm classes & domain-ints
    check_file(elm_classes_file)
    check_file(elm_int_dom_file)

    ## Check/Edit species files.
    print "Species files:"
    for f in [uni_text_file, uni_fasta_file, pfam_matches_file, biogrid_file]:
        check_file(f)
    # Check/Generate ELM annotation
    if not os.path.isfile(elm_hits_file):
        find_all_elms.main(uni_fasta_file, elm_classes_file,
                           print_out=True, outfile=elm_hits_file)
        print_status(elm_hits_file, "created")
    else:
        print_status(elm_hits_file, "exists")
    # Check Dom-Dom Association file
    flag=0
    if os.path.isfile(pfam_assoc_file):
        print_status(pfam_assoc_file, "exists")
        flag=1
    flag2=0
    if os.path.isfile(iprets_file):
        print_status(iprets_file, "exists")
        flag2=1

    # Get protein ID-dictionary, sequences and masks.
    prot_dict, masks, ref_prot, modres, mutagen, regions = get_protein_data_from_uniprot_text(uni_text_file)
    # Get Pfam domains.
    pfams, pfam_names = get_pfam_doms(pfam_matches_file, prot_dict)
    # Get elm info.
    elms = get_interacting_linear_motifs(elm_int_dom_file,
                                         elm_hits_file,
                                         masks,
                                         max_overlap=2)

    if os.path.isfile(prot_data_file) and force==False:
        print_status(prot_data_file, "exists")
    else:
        # Get PTMs
        if sp in ["Hsa", "Mmu"]:
            check_file(psp_file)
            ptms = get_ptms(psp_file, prot_dict, modres)
        else:
            ptms = modres
        # Create protein_data.json for Mongo
        make_protein_data_json(prot_dict, pfams, elms,
                               ptms, mutagen, regions,
                               prot_data_file, mode="mongo")


    # Get PP interaction from BioGRID file
    bio_ppi = defaultdict(lambda: defaultdict(set))
    tot_ints = 0
    with open_file(biogrid_file) as f:
        for line in f:
            t = line.rstrip().split("\t")
            if line[0]!="#":
                biogrid_int_id = t[0]
                geneA, geneB = t[7], t[8]
                if (geneA in prot_dict["AC"] and geneB in prot_dict["AC"]):
                    protA = list(prot_dict["AC"][geneA])[0]
                    protB = list(prot_dict["AC"][geneB])[0]
                    if protA < protB:
                        a = protA
                        b = protB
                    else:
                        a = protB
                        b = protA
                        if b not in bio_ppi[a]:
                            tot_ints+=1
                        bio_ppi[a][b].add(biogrid_int_id)

    if flag==0:
        ## Calculate Dom-Dom Interactions (Association Method)
        ## using Swissprot proteins as reference set

        if sp=="Hsa":
            # Source of PPI: HIPPIE (until I find something better)
            # not all ids from the hippie db matched a uniprot_ac, thus Uniprot mapping
            # tool was used and the results need to be read before:
            hippie_map = defaultdict(set)
            hippie_pfams = {}
            check_file(hippie_uniprot_mapping)
            check_file(hippie_file)
            with open_file(hippie_uniprot_mapping) as f:
                for line in f:
                    t = line.rstrip().split("\t")
                    uni_ac = t[1]
                    if uni_ac in prot_dict["AC"]:
                        uni_ac = list(prot_dict["AC"][uni_ac])[0] ## always only 1. I checked
                        for x in t[0].split(","):
                            hippie_map[x].add(uni_ac)
                            hippie_pfams[x] = pfams[uni_ac].keys()

            # Parameters:
            # (from HIPPIE) Medium confidence = 0.63 / High confidence = 0.73
            hippie_score = 0.63
            ppi = defaultdict(dict)
            with open_file(hippie_file) as f:
                for line in f:
                    t = line.rstrip().split("\t")
                    protsA = t[0].split(",")
                    protsB = t[2].split(",")
                    score = float(t[4])
                    ## Score threshold to keep PPI
                    ## (from HIPPIE) Medium confidence = 0.63 / High confidence = 0.73
                    if score < hippie_score:
                        continue
                    pmids = 0
                    if "pmids" in t[5]:
                        pmids = len(t[5].split(";")[1].split(":")[1].split(","))

                    for protA in protsA:
                        for protB in protsB:
                            if protA!="" and protB!="":
                                if protA < protB:
                                    a = protA
                                    b = protB
                                else:
                                    a = protB
                                    b = protA
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
                for b in ppi[a]:
                    if a in hippie_map and b in hippie_map:
                        total_prts.add(a)
                        total_prts.add(b)
                        total_pp_pairs += 1

                        # Count        # with open_file(iprets_file, "w") as out:
        #     cols = ["#Gene1","PDB1","Blast-E1","Blast-PCID1","qstart1","qend1","pdbstart1","pdbend1",
        #             "Gene2","PDB2","Blast-E2","Blast-PCID2","qstart2","qend2","pdbstart2","pdbend2",
        #             "i2-raw","rand","rand-mean","rand-sd","Z","p-value","not-sure1","not-sure2"]
        #     out.write("\t".join(cols)+"\n") domain pair frequency
                        for pfam_a in hippie_pfams[a]:
                            for pfam_b in hippie_pfams[b]:
                                if pfam_a >= pfam_b:
                                    pfam_pair_count[pfam_a][pfam_b] += 1
                                else:
                                    pfam_pair_count[pfam_b][pfam_a] += 1


            ## Count individual domain frequency
            pfam_count = defaultdict(int)
            ref_pfam = set()
            for prot in total_prts:
                if len(hippie_pfams[prot])>1:
                    ref_pfam.add(prot)
                for pfam in hippie_pfams[prot]:
                    pfam_count[pfam] += 1

            N_total = len(ref_prot)


        else:
            ppi = bio_ppi
            # Count domain pair frequencies in non redundant protein pairs
            # (domains are counted only once when they are repeated in the same protein)
            pfam_pair_count = defaultdict(lambda: defaultdict(int))
            total_prts = set()
            total_pp_pairs = 0
            for a in ppi:
                for b in ppi[a]:
                    total_prts.add(a)
                    total_prts.add(b)
                    total_pp_pairs += 1

                    # Count domain pair frequency
                    for pfam_a in pfams[a]:
                        for pfam_b in pfams[b]:
                            if pfam_a >= pfam_b:
                                pfam_pair_count[pfam_a][pfam_b] += 1
                            else:
                                pfam_pair_count[pfam_b][pfam_a] += 1

            # Count individual domain frequency
            pfam_count = defaultdict(int)
            for prot in total_prts:
                for pfam in pfams[prot]:
                    pfam_count[pfam] += 1

            N_total = len(total_prts)

        # Compute domain-domain statistics
        # (only those appearing in some PPI. To the rest we can assign the lowest
        # LO score)
        # Minimum number of protein pairs with the domain-domain signature
        min_npair = 1
        # domain count could be another parameter
        lo = set()
        with gzip.open(pfam_assoc_file, "wb") as out:
            out.write("\t".join(["dom_ac_a","dom_ac_b","dom_name_a",
                "dom_name_b","dom_n_a","dom_n_b","obs","exp","or","lo"])+"\n")
            for pfam_a in sorted(pfam_pair_count):
                for pfam_b in sorted(pfam_pair_count[pfam_a]):
                    n_a = pfam_count[pfam_a]
                    n_b = pfam_count[pfam_b]
                    f_exp = (n_a * n_b) / float(N_total * N_total)
                    exp = f_exp * total_pp_pairs
                    obs = float(pfam_pair_count[pfam_a][pfam_b])
                    oddsratio = obs / exp
                    log2 = math.log(oddsratio, 2)
                    if obs >= min_npair:
                        out.write("\t".join([pfam_a, pfam_b,
                                        pfam_names[pfam_a], pfam_names[pfam_b],
                                        str(n_a), str(n_b), str(obs), str(exp),
                                        str(oddsratio), str(log2)])+"\n")
                        lo.add(log2)

        print_status(pfam_assoc_file, "created")

    # print "#Reference proteome =",len(ref_prot),"with pfams =",len(ref_pfam)
    # print "#Total PPI pairs =", total_pp_pairs, "- involving", len(total_prts)
    # print "#Lowest / highest LO =",sorted(list(lo))[0],"/", sorted(list(lo))[-1]

    if flag2 != 0:
        sys.exit()
        # Run InterPreTS using BioGRID pairs

        # with open_file(iprets_file, "w") as out:
        #     cols = ["#Gene1","PDB1","Blast-E1","Blast-PCID1","qstart1","qend1","pdbstart1","pdbend1",
        #             "Gene2","PDB2","Blast-E2","Blast-PCID2","qstart2","qend2","pdbstart2","pdbend2",
        #             "i2-raw","rand","rand-mean","rand-sd","Z","p-value","not-sure1","not-sure2"]
        #     out.write("\t".join(cols)+"\n")

        # seq = {}
        # with open_file(uni_fasta_file) as f:
        #     for record in SeqIO.parse(f, "fasta"):
        #         seq[str(record.id).split("|")[1]]=str(record.seq)

        n = 0
        for a in sorted(bio_ppi):
            input_seqs = {}
            input_seqs[a] = prot_dict["seq"][a]
            pairs = []
            for b in sorted(bio_ppi[a]):
                pairs.append((a,b))
                input_seqs[b] = prot_dict["seq"][b]
                n+=1

            print "[{}] Running InterPrets for {} and {} partners: {}".format(datetime.datetime.now(), a, str(len(bio_ppi[a])), ";".join(bio_ppi[a].keys()))
            iprets_dir = sp_dir+"iprets/"
            # run_interprets.main(input_seqs, iprets_dir, iprets_file,
            #                     these_pairs=pairs)

            print "[{}] Done {} out of {}".format(datetime.datetime.now(), n, tot_ints)

        print_status(iprets_file, "created")

    ## Preprocessing InterPreTS info
    # chain2uni = {}
    ## WRONG
    ## I cant use fixed PDB-UniAcc conversion because InterPreTS uses homology
    ## to find templates!
    ## Don't do anything for the moment
    # s = set()
    # with open_file(pdbchain2uniprot) as f:
    #     for line in f:
    #         t = line.upper().rstrip().split()
    #         if len(t) == 3:
    #             pdb = t[0]
    #             chain = t[1]
    #             uni_ac = t[2]
    #             if uni_ac in prot_dict["AC"]:
    #                 for ac in prot_dict["AC"][uni_ac]:
    #                     chain2uni["pdb|"+pdb+"|"+chain].add(ac)

    # ints = defaultdict(lambda: defaultdict(list))
    # with open_file(iprets_file) as f: ## lines in this file are already sorted from highest-to-lowest Z-score!!
    #     for line in f:
    #         t = line.rstrip().split("\t")
    #         if line[0] == "#":
    #             columns = t
    #         else:
                # pdb1 = t[1]
                # pdb2 = t[9]
                # if pdb1 not in chain2uni:
                #     s.add(pdb1)
                # if pdb2 not in chain2uni:
                #     s.add(pdb2)
                # if pdb1 in chain2uni and pdb2 in chain2uni:
                #     for uni_ac1 in chain2uni[pdb1]:
                #         for uni_ac2 in chain2uni[pdb2]:
                #             print t[0]+"\t"+uni_ac1+"\t"+"\t".join(t[1:9])+"\t"+uni_ac2+"\t"+"\t".join(t[9:])

                # gene1, gene2 = t[0], t[8]
                # if len(ints[gene1][gene2])>0:
                #     compare_ints(t, gene1, gene2, ints)
                # else:
                #     ints[gene1][gene2].append("\t".join(t))


    # print  "# not found:", len(s)
    # for pdb in s:
    #     print pdb


if __name__ == "__main__":
    force=False
    if "-force" in sys.argv:
        force=True
    main(sp=sys.argv[1], data_dir="static/data/", force=force)
    sys.exit()
