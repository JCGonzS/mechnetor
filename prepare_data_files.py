#!/usr/bin/env python

import sys, re, os, gzip, math, pprint, datetime, random, string
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

def check_file_exists(somefile):
    if os.path.isfile(somefile):
        print_status(somefile, "exists")
        return True
    else:
        sys.exit("FileError: '"+somefile+"' does not exist.\n"+
            "Make sure the file is located and named correctly\n")

def print_status(thefile, status):
    l = 80 - len(thefile) - len(status)
    print "'"+thefile+"'"+"."*l+status

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

def extract_protein_data_from_uniprot_text(uniprot_file):
    """Extracts information from a UniProt text file.
    Requires Bio::SwissProt module
    """
    # all_genes = defaultdict(set)
    # for record in SwissProt.parse(open_file(uniprot_file)):
    #     if record.gene_name.strip():
    #         ## Main gene names
    #         for name in [match.split()[0] for match in re.findall("Name=([^;]+)",
    #                                                         record.gene_name)]:
    #             name = name.split("_")[0]
    #             all_genes[name].add(record.data_class)

    D = defaultdict(lambda: defaultdict(list))
    alt_ids = defaultdict(lambda: defaultdict(list))
    masks = {}
    modres = defaultdict(lambda: defaultdict(dict))
    uni_features = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    region = defaultdict(list)
    pdb2uni = defaultdict(set)
    reviewed = set()
    previous_data_class = ""
    for record in SwissProt.parse(open_file(uniprot_file)):
        if previous_data_class == "Unreviewed" and record.data_class == "Reviewed":
            print "WARNING: UniProt entries not sorted Reviewed->Unreviewed"

        uni_id = record.entry_name.upper()
        accs = record.accessions
        uni_ac = accs[0].upper()
        des = re.search("Name: Full=([^;|{]+)", record.description).group(1)
        genes, syns = [uni_id], []
        if record.gene_name.strip():

            ## Main gene names
            names = [match.split()[0] for match in re.findall("Name=([^;]+)",
                                                            record.gene_name)]
            if len(names) > 0:
                genes = []
                for name in names:
                    name = name.split("_")[0]
                    if name not in genes:
                        genes.append(name)

            ## ORF & Locus names
            orfnames = []
            for match in (re.findall("ORFNames=([^;]+)", record.gene_name)+
                    re.findall("OrderedLocusNames=([^;]+)", record.gene_name)):
                for name in match.split(", "):
                    name = re.search("([^{]+)", name).group(1).replace(" ","")
                    if not name.startswith("ECO:"):
                        for name2 in name.split("/"):
                            # if name2 not in all_genes:
                            orfnames.append(name2)
            if len(orfnames) > 0 and len(names) == 0:
                genes = [orfnames[0]]
            syns = orfnames

            ## Other Synonyms
            for match in re.findall("Synonyms=([^;]+)", record.gene_name):
                for syn in match.split(", "):
                    syn = re.search("([^{]+)", syn).group(1).replace(" ","")
                    if not syn.startswith("ECO:"):
                        # if syn not in all_genes:
                        syns.append(syn)

        main_gene = genes[0]


        ## Rememeber: all Reviewed entries come first, then all Unreviewed ones
        if record.data_class == "Reviewed":
            for dic, val in zip(["AC", "ID"], [uni_ac, uni_id]):
                for key in [uni_ac, uni_id, main_gene, main_gene.upper()]:
                    reviewed.add(key)
                    if val not in D[dic][key]:
                        D[dic][key].insert(0, val)

                for gene in genes[1:]+syns:
                    for gn in [gene, gene.upper()]:
                        reviewed.add(gn)
                        if val not in D[dic][gn]:
                            D[dic][gn].append(val)

            for key in [uni_ac, uni_id]:
                if main_gene not in D["GN"][key]:
                    D["GN"][key].insert(0, main_gene)
                for gene in genes[1:]:
                    if gene not in D["GN"][key]:
                        D["GN"][key].append(gene)

        else: # Unreviewed
            for dic, val in zip(["AC", "ID"], [uni_ac, uni_id]):
                for key in [uni_ac, uni_id, main_gene, main_gene.upper()]:
                    if val not in D[dic][key]:
                        if key not in reviewed:
                            D[dic][key].insert(0, val)
                        else:
                            D[dic][key].append(val)

                for gene in genes[1:]+syns:
                    for gn in [gene, gene.upper()]:
                        if val not in D[dic][gn]:
                            D[dic][gn].append(val)

            for key in [uni_ac, uni_id]:
                if main_gene not in D["GN"][key]:
                    if main_gene not in reviewed:
                        D["GN"][key].insert(0, main_gene)
                    else:
                        D["GN"][key].append(main_gene)
                for gene in genes[1:]:
                    if gene not in D["GN"][key]:
                        D["GN"][key].append(gene)

        alt_ids[uni_ac] = accs[1:]+syns
        D["des"][uni_ac] = des
        D["dc"][uni_ac]  = record.data_class
        D["seq"][uni_ac] = record.sequence
        masks[uni_ac] = ["0"] * len(record.sequence)

        for feat in record.features:
            start, end = feat[1], feat[2]
            info = feat[3].split(" {ECO")[0]
            if feat[0]=="MOD_RES": # start == end, always
                res = record.sequence[start-1]
                if "Phospho" in feat[3].split(";")[0]:
                    modres[uni_ac][str(start)+"-"+res]["p"] = ["mod_res"]
                elif "acetyl" in feat[3].split(";")[0]:
                    modres[uni_ac][str(start)+"-"+res]["ac"] = ["mod_res"]

            elif feat[0] in ["MUTAGEN", "VARIANT", "METAL", "BINDING"]:
                pos = str(start)+"-"+str(end)
                uni_features[feat[0]][uni_ac][pos].append(info)

            elif feat[0]=="REGION":
                region[uni_ac].append(
                    {"start": start,
                     "end": end,
                     "info": info
                    }
                )

            else:
                print feat

        # for ref in record.cross_references:
        #     if "PDB" in ref:
        #         chains = set()
        #         for x in ref[-1].split(", "):
        #             for y in x.split("=")[0].split("/"):
        #                 chains.add(y)
        #         for c in chains:
        #             pdb2uni["pdb|"+ref[1]+"|"+c].add(uni_ac)
    return D, alt_ids, masks, modres, uni_features, region

def calculate_overlap(start, end, mask):
    sub_mask = mask[start-1:end]
    overlap = sub_mask.count("1") / float(len(sub_mask))
    return overlap

def fill_mask(start, end, mask):
    length = end - start + 1
    for i in range(start-1, end):
        mask[i] = "1"
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

def extract_pfam_doms(pfam_file, prot_dict, max_eval=999):
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
                pfam_names[pfam_ac] = pfam_name

                ## Keep annotation for primary UniProt accessions
                if uni_ac in prot_dict["seq"]:
                    pfams[uni_ac][pfam_ac].add((start, end, domain_e_val))
                ## Save annotation for 2ary accessions in different dictionary
                elif uni_ac in prot_dict["AC"]:
                    pfams_temp[uni_ac][pfam_ac].add((start, end, domain_e_val))

    ## Transfer annotation to those main accessions which were not annotated
    for alt_ac in pfams_temp:
        for uni_ac in prot_dict["AC"][alt_ac]:
            ## If the associated main accession is not annotated already:
            if uni_ac not in pfams:
                pfams[uni_ac] = pfams_temp[alt_ac]

    return pfams, pfam_names

def get_elm_names(elm_classes_file):
    elm_names = {}
    with open_file(elm_classes_file) as f:
        for line in f:
            if line.startswith("ELM"):
                t = line.rstrip().split("\t")
                elm_names[t[0]] = t[1]
    return elm_names

def extract_elm_interactions(elm_intdom_file, elm_names):
    elm_acc = dict((v, k) for k, v in elm_names.iteritems())
    elm_int = defaultdict(set)
    with open_file(elm_intdom_file) as f:
        for line in f:
            if line[0] != "#":
                t = line.rstrip().split("\t")
                if t[0] in elm_acc:
                    elm, dom = elm_acc[t[0]], t[1]
                    elm_int[elm].add(dom)
                    elm_int[dom].add(elm)

    return elm_int

def extract_linear_motifs(elm_hits_file, elm_names, masks,
                                  max_overlap, max_eval=1):
    """ Reads pre-generated file with ELMs annotated for each protein
    """

    elms = defaultdict(lambda: defaultdict(set) )
    elms_temp = defaultdict(lambda: defaultdict(set) )
    with open_file(elm_hits_file) as f:
        for line in f:
            if line[0] != "#" and line.strip():
                t = line.rstrip().split("\t")
                elm_name, elm_ac = t[0], t[1]
                uni_ac, uni_id = t[2].split("|")[1:]
                elm_start, elm_end = int(t[4]), int(t[5])
                prob_score = float(t[6]) # probability score based on the
                                         # combined expected frequencies of
                                         # the AAs in the regular expression

                # overlap = calculate_overlap(elm_start, elm_end, masks[uni_ac])
                # if overlap <= max_overlap:
                #     masks[uni_ac] = fill_mask(elm_start, elm_end, masks[uni_ac])
                elms[uni_ac][elm_ac].add((elm_start, elm_end))

    return elms

def edit_pfam_dat(pfam_dat_file, out_file):
    """ Writes a simple TSV file with a few selected columns from Pfam data.
    """
    with open_file(out_file, "w") as out:
        cols = ["Accession", "Identifier", "Description"]
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

def edit_3did(db_file, out_file):
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
                    if pfam_acc_a < pfam_acc_b:
                        row = [pfam_name_a, pfam_acc_a, pfam_name_b, pfam_acc_b,
                               ";".join(sorted(list(pdbs)))]
                    else:
                        row = [pfam_name_b, pfam_acc_b, pfam_name_a, pfam_acc_a,
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


def extract_biogrid_interactions(biogrid_file, prot_dict):
    """From a BioGRID proteom file, extracts:
        - All protein-protein interactions (mapped to UniProt accessions)
        - UniProt accession - BioGRID protein identifier diccionary
        - Total number of PPI
    """
    biogrid_ppi = defaultdict(lambda: defaultdict(set))
    biogrid_ppi_per_prot = defaultdict(lambda: defaultdict(int))
    biogrid_prot_id = {}
    protein_not_found = set()
    with open_file(biogrid_file) as f:
        for line in f:
            t = line.rstrip().split("\t")
            if line[0] != "#":
                interaction_id = t[0]
                idA, idB = t[3], t[4]
                geneA,geneB = t[7].upper(), t[8].upper()
                synsA,synsB = t[9].upper().split("|"), t[10].upper().split("|")
                pubmed, throughput = t[14], t[17].replace(" Throughput","")
                accA, accB = "-", "-"
                for protA in [geneA]+synsA:
                    if protA in prot_dict["AC"]:
                        accA = prot_dict["AC"][protA][0]
                        break
                for protB in [geneB]+synsB:
                    if protB in prot_dict["AC"]:
                        accB = prot_dict["AC"][protB][0]
                        break

                if accA != "-" and accB != "-":
                    biogrid_prot_id[accA] = idA
                    biogrid_prot_id[accB] = idB
                    biogrid_ppi_per_prot[accA][accB] += 1
                    biogrid_ppi_per_prot[accB][accA] += 1
                    if accA < accB:
                        biogrid_ppi[accA][accB].add(
                                    "BioGRID:"+interaction_id+":"+pubmed+":"+throughput)
                    else:
                        biogrid_ppi[accB][accA].add(
                                    "BioGRID:"+interaction_id+":"+pubmed+":"+throughput)
                else:
                    for gene in [geneA, geneB]:
                        if gene not in prot_dict["AC"]:
                            protein_not_found.add(gene)

    print "From BioGRID, not found:", len(protein_not_found), "proteins"
    # for prot in protein_not_found:
    #     print prot
    return biogrid_ppi, biogrid_ppi_per_prot, biogrid_prot_id

def create_protein_data_json(prot_dict, alt_ids, pfams, elms, ptms,
                           uni_features, regions, bio_id, bio_set,
                           outfile_name, mode="mongo"):

    pp = pprint.PrettyPrinter(indent=4)
    json_data = {}
    with open_file(outfile_name, "w") as out:
        for uni_ac in prot_dict["seq"]:

            seq = prot_dict["seq"][uni_ac]
            protein_data = {
                    "uni_ac": uni_ac,
                    "uni_id": list(prot_dict["ID"][uni_ac])[0], # always, only 1
                    "genes": list(prot_dict["GN"][uni_ac]),
                    "alt_ids": alt_ids[uni_ac],
                    "description": prot_dict["des"][uni_ac],
                    "data_class": prot_dict["dc"][uni_ac],
                    "length" : len(seq),
                    "sequence" : seq,
                    "pfams" : [],
                    "elms" : [],
                    "phosphorylation" : [],
                    "acetylation" : [],
                    "uni_features" : [],
                    "regions" : []
            }

            if uni_ac in bio_id:
                protein_data["biogrid_id"] = bio_id[uni_ac]
                sorted_ints = sorted(bio_set[uni_ac], key=bio_set[uni_ac].get,
                                     reverse=True)
                protein_data["biogrid_interactors"] = sorted_ints
            else:
                protein_data["biogrid_id"] = "NA"
                protein_data["biogrid_interactors"] = []

            for pfam_ac in sorted(pfams.get(uni_ac, [])):
                for (start, end, e_val) in sorted(pfams[uni_ac][pfam_ac],
                                                  key=lambda x: int(x[0])):
                    protein_data["pfams"].append(
                        {
                            "acc" :   pfam_ac,
                            "start" : int(start),
                            "end" :   int(end),
                            "e-val" : float(e_val)
                        })

            for elm_ac in sorted(elms.get(uni_ac, [])):
                for (start, end) in sorted(elms[uni_ac][elm_ac],
                                                key=lambda x: int(x[0])):
                    protein_data["elms"].append(
                        {
                            "acc" :   elm_ac,
                            "start" : int(start),
                            "end" :   int(end),
                            "seq" :   seq[int(start)-1:int(end)]
                        })

            for pos_res in ptms.get(uni_ac, []):
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

            for feat in uni_features:
                for pos in uni_features[feat].get(uni_ac, []):
                    protein_data["uni_features"].append(
                        {
                            "feat": feat,
                            "pos": pos,
                            "info": ";".join(uni_features[feat][uni_ac][pos])
                        }
                    )

            protein_data["regions"] = regions.get(uni_ac, [])

            if mode == "mongo":
                out.write(str(protein_data)+"\n")
            elif mode == "json":
                json_data[uni_ac] = protein_data

        if mode == "json":
            json.dump(json_data, out)


def hasNumbers(inputString):
	return any(char.isdigit() for char in inputString)

def assign_Number(inputString):
	if ";" in inputString:
		return inputString.split(";")[0]
	else:
		return inputString

def extract_ptms(psp_file, prot_dict, modres):
    psp = defaultdict(set)
    psp_temp = defaultdict(set)
    with open_file(psp_file) as f:
        for line in f:
            t = line.rstrip().split("\t")
            pos, mod = t[4].split("-")
            res, pos = re.search("(\w)([\d]+)", pos).group(1,2)
            lt_lit = "."
            ms_lit = "."
            ms_cst = "."
            cst_cat = "."
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

def get_HIPPIE_UniProt_map(map_file, prot_dict):
    """
    # Source of PPI: HIPPIE
    # not all ids from the hippie db matched a uniprot_ac, thus Uniprot mapping
    # tool was used and the results need to be read before:
    """
    hippie_map = defaultdict(set)
    with open_file(map_file) as f:
        for line in f:
            if (line.startswith("yourlist") or "Deleted." in line
                or "Homo sapiens" not in line):
                continue
            t = line.rstrip().split("\t")
            uni_ac = t[1]
            uni_ac = list(prot_dict["AC"][uni_ac])[0] #always 1 (checked)
            for protein in t[0].split(","):
                hippie_map[protein].add(uni_ac)

    return hippie_map

def extract_HIPPIE_interactions(hippie_file, hippie_map, score_cutoff=0.63):
    """
        (from HIPPIE) Medium confidence = 0.63 / High confidence = 0.73
    """
    ppi = defaultdict(set)
    with open_file(hippie_file) as f:
        for line in f:
            t = line.rstrip().split("\t")
            protsA = t[0].split(",")
            protsB = t[2].split(",")
            score = float(t[4])
            if score >= score_cutoff:
                for protA in protsA:
                    if protA in hippie_map:
                        a = list(hippie_map[protA])[0]
                        for protB in protsB:
                            if protB in hippie_map:
                                b = list(hippie_map[protB])[0]
                                if b < a:
                                    a = list(hippie_map[protB])[0]
                                    b = list(hippie_map[protA])[0]
                                ppi[a].add(b)
    return ppi

def extract_3did_interactions(edited_3did_file):
    ints = defaultdict(set)
    with open_file(edited_3did_file) as f:
        for i, line in enumerate(f):
            if i > 0:
                t = line.rstrip().split("\t")
                ac_a, ac_b = t[1], t[3]
                ints[ac_a].add(ac_b)
                ints[ac_b].add(ac_a)

    return ints

def create_ppi_database(outfile, ppi, prot_dict):

    with open_file(outfile, "w") as out:
        out.write("\t".join(["Acc_A", "Gene_A", "Acc_B", "Gene_B",
                            "Source:ID:PubMedID:Throughput"])+"\n")

        for acc_a in ppi:
            for acc_b in ppi[acc_a]:
                gene_a = prot_dict["GN"][acc_a][0]
                gene_b = prot_dict["GN"][acc_b][0]
                cols = [acc_a, gene_a, acc_b, gene_b,
                        ";".join(ppi[acc_a][acc_b])]
                out.write("\t".join(cols)+"\n")

def count_protein_and_domain_pairs(ppi, pfams, elms):
    # From non-redundant PPI set, get:
    #  - total number of proteins
    #  - total number of interactions (or protein pairs)
    all_proteins = set()
    total_interactions = 0
    ele_pair_count = defaultdict(lambda: defaultdict(int))
    for a in ppi:
        for b in ppi[a]:
            all_proteins.add(a)
            all_proteins.add(b)
            total_interactions += 1

            # - domain/elm pair counts
            ele_pairs = set()
            for ele_a in sorted(pfams[a].keys() + elms[a].keys()):
                for ele_b in sorted(pfams[b].keys() + elms[b].keys()):
                    if ele_a <= ele_b:
                        eleA, eleB = ele_a, ele_b
                    else:
                        eleA, eleB = ele_b, ele_a

                    if (eleA,eleB) not in ele_pairs:
                        ele_pair_count[eleA][eleB] += 1
                        ele_pairs.add((eleA, eleB))

    return all_proteins, total_interactions, ele_pair_count

def count_individual_element_frequency(all_proteins, pfams, elms):
    element_count = defaultdict(int)
    for protein in all_proteins:
        for pfam in pfams[protein]:
            element_count[pfam] += 1
        for elm in elms[protein]:
            element_count[elm] += 1

    return element_count

def create_prob_file(out_file, total_proteins, total_interactions,
                     ele_pair_count, ele_count, ele_names,
                     ints_3did, ints_elm, min_npair=1):
    """ Calculate element-element* probabilities and log-odds
        *elements are Pfam domains and ELMs

        Only for those appearing in some PPI. To the rest we can assign the
        lowestLO score)

        Parameters:
            - Minimum number of protein pairs with the domain-domain signature
    """

    lo = set()
    # with gzip.open(out_file, "wb") as out:
    with open_file(out_file, "w") as out:
        line = ("{:<8} {:<10} {:<24} {:<6} {:<10}"+
                " {:<10} {:<24} {:<6} {:<10}"+
                " {:<9} {:<12} {:<8} {:<9} {:<12} {:<8}").format(
                "#Type", "Acc_A", "Name_A", "n_A", "Prob_A",
                          "Acc_B", "Name_B", "n_B", "Prob_B",
                "Prob_A&B", "Exp", "Obs", "P-val", "Ratio", "LO")
        # print line
        cols = ["#Type", "Acc_A", "Name_A", "n_A", "Prob_A",
                         "Acc_B", "Name_B", "n_B", "Prob_B",
                "Prob_A&B", "Exp", "Obs", "P-val", "Ratio", "LO"]
        out.write("\t".join(cols)+"\n")


        for ele_a in sorted(ele_pair_count):
            for ele_b in sorted(ele_pair_count[ele_a]):
                if ele_b in ints_3did.get(ele_a, []):
                    cat = "3DID"
                elif ele_b in ints_elm.get(ele_a, []):
                    cat = "ELM"
                else:
                    cat = "Pred"

                prob_a = float(ele_count[ele_a]) / total_proteins
                prob_b = float(ele_count[ele_b]) / total_proteins
                prob_ab = prob_a * prob_b
                exp = prob_ab * total_interactions
                obs = ele_pair_count[ele_a][ele_b]
                p_val = float(obs)/total_interactions
                ratio = float(obs) / exp
                log2 = math.log(ratio, 2)

                # if obs >= min_npair:
                line = ("{:<8} {:<10} {:<24} {:<6} {:<.3E} "+
                     " {:<10} {:<24} {:<6} {:<.3E}"+
                     " {:<.3E} {:<12.3f} {:<8} {:<.3E} {:<12.3f} {:<8.4f}").format(
                    cat, ele_a, ele_names[ele_a], ele_count[ele_a], prob_a,
                    ele_b, ele_names[ele_b], ele_count[ele_b], prob_b,
                    prob_ab, exp, obs, p_val, ratio, log2)
                # print line
                cols = [cat,
                        ele_a, ele_names[ele_a], str(ele_count[ele_a]), str(prob_a),
                        ele_b, ele_names[ele_b], str(ele_count[ele_b]), str(prob_b),
                        str(prob_ab), str(exp), str(obs), str(p_val),
                        str(ratio), str(log2)]
                out.write("\t".join(cols)+"\n")

def create_interprets_file(outfile, ppi, prot_dict, iprets_dir):

    with open_file(outfile, "w") as out:
        cols = ["#Gene1","PDB1","Blast-E1","Blast-PCID1","qstart1","qend1",
                "pdbstart1","pdbend1","Gene2","PDB2","Blast-E2","Blast-PCID2",
                "qstart2","qend2","pdbstart2","pdbend2","i2-raw","rand",
                "rand-mean","rand-sd","Z","p-value","not-sure1","not-sure2"]
        out.write("\t".join(cols)+"\n")

    for a in sorted(ppi):
        input_seqs = {}
        input_seqs[a] = prot_dict["seq"][a]
        pairs = []
        for b in sorted(ppi[a]):
            pairs.append((a,b))
            input_seqs[b] = prot_dict["seq"][b]

        ide = ''.join(random.choice(string.ascii_uppercase +
                        string.ascii_lowercase + string.digits) for _ in range(8))
        print "Running InterPrets for {} and {} partners: {}".format(a,
                                        str(len(ppi[a])), ";".join(ppi[a]))
        try:
            run_interprets.main(input_seqs, iprets_dir, outfile, ide,
                                these_pairs=pairs, fasta_as_input=False,
                                verbose=False)
        except:
            print "InterPreTS failed"
        print "[IntePreTs finished for {} protein pairs]".format(len(pairs))


def main( SP="Hsa",
          DATA_DIR="",
          max_overlap=0.2,
          no_overlap_between_doms_and_lms=True,
          force_prot_data=False,
          force_ppi=False,
          force_probs=False,
          force_iprets=False
        ):

    ### 1. Define required data files
    UNI_VERSION     = "Nov2019" # date of download
    PFAM_VERSION    = "r32.0"   # release 32
    ELM_VERSION     = "Oct2019"
    BIOGRID_VERSION = "3.5.178" # release
    PSP_VERSION     = "Mar2019"

    ## Common files:
    COM_DIR              = DATA_DIR+"common/"
    PFAM_DAT_FILE        = COM_DIR+"Pfam-A.hmm_"+PFAM_VERSION+".dat.gz"
    EDITED_PFAM_DAT_FILE = COM_DIR+"Pfam-A.hmm_"+PFAM_VERSION+".tsv.gz"
    FLAT_3DID_FILE       = COM_DIR+"3did_flat-2018_04.gz"
    EDITED_3DID_FILE     = COM_DIR+"3did_flat_edited-2018_04.tsv.gz"
    # pdbchain2uniprot   = COM_DIR+"pdbsws_chain.txt.gz"
    ELM_CLASSES_FILE     = COM_DIR+"elm_classes_"+ELM_VERSION+".tsv"
    ELM_INTDOM_FILE      = (COM_DIR+"elm_interaction_domains_"
                            +ELM_VERSION+"_reviewed.tsv")

    ## Species files:
    SP_DIR            = DATA_DIR+"species/"+SP+"/"
    UNI_TEXT_FILE     = SP_DIR+"uniprot_"+UNI_VERSION+"_proteome_"+SP+".txt.gz"
    UNI_FASTA_FILE    = (SP_DIR+"uniprot_"+UNI_VERSION+"_proteome_"+SP
                        +".fasta.gz")
    PFAM_MATCHES_FILE = SP_DIR+"pfamA_"+PFAM_VERSION+"_matches_"+SP+".tsv.gz"
    ELM_HITS_FILE     = SP_DIR+"elm_hits_"+ELM_VERSION+"_"+SP+".tsv.gz"
    BIOGRID_FILE      = (SP_DIR+"BIOGRID-ORGANISM-"+SP+"-"+BIOGRID_VERSION
                        +".tab2.txt.gz")
    PPI_FILE          = SP_DIR+"ppi_db_"+SP+".tsv.gz"
    ASSOC_PROB_FILE   = SP_DIR+"prot_ele_association_prob_"+SP+".tsv.gz"
    PROT_DATA_FILE    = SP_DIR+"protein_data_"+SP+"_mongo.json.gz"
    IPRETS_FILE       = SP_DIR+"interprets_results_"+SP+".txt.gz"
    IPRETS_HSA_FILE   = SP_DIR+"human_aaa_biogrid_i2.txt.gz" # Hsa
    HIPPIE_FILE       = SP_DIR+"hippie_v2-2.tsv.gz" # Hsa
    HIPPIE_MAP_FILE   = SP_DIR+"hippie_v2-2_uniprot_mapping_table.tsv.gz" # Hsa
    PSP_FILE          = SP_DIR+"PSP_ptms_"+PSP_VERSION+"_"+SP+".tsv.gz" # Only Hsa and Mmu

    ### 2. Check existence/make required data files
    ## Check/Edit common files:
    # Pfam.dat
    if os.path.isfile(EDITED_PFAM_DAT_FILE):
        print_status(EDITED_PFAM_DAT_FILE, "exists")
    elif check_file_exists(PFAM_DAT_FILE):
        edit_pfam_dat(PFAM_DAT_FILE, EDITED_PFAM_DAT_FILE)
        print_status(EDITED_PFAM_DAT_FILE, "created")
    # 3did flat
    if os.path.isfile(EDITED_3DID_FILE):
        print_status(EDITED_3DID_FILE, "exists")
    elif check_file_exists(FLAT_3DID_FILE):
        edit_3did(FLAT_3DID_FILE, EDITED_3DID_FILE)
        print_status(EDITED_3DID_FILE, "created")
    # elm classes & domain-ints
    check_file_exists(ELM_CLASSES_FILE)
    check_file_exists(ELM_INTDOM_FILE)

    ## Check/Edit species files:
    for f in [UNI_TEXT_FILE, UNI_FASTA_FILE, PFAM_MATCHES_FILE, BIOGRID_FILE]:
        check_file_exists(f)
    # Check/Generate ELM annotation
    if os.path.isfile(ELM_HITS_FILE):
        print_status(ELM_HITS_FILE, "exists")
    else:
        find_all_elms.main(UNI_FASTA_FILE, ELM_CLASSES_FILE,
                           print_out=True, outfile=ELM_HITS_FILE)
        print_status(ELM_HITS_FILE, "created")

    ### 3. Get various data
    # Extract protein ID-dictionary, sequences and masks.
    (prot_dict, alt_ids, masks, ptms, uni_features,
        regions) = extract_protein_data_from_uniprot_text(UNI_TEXT_FILE)
    sys.exit()
    # Extract Pfam domains.
    pfams, pfam_names = extract_pfam_doms(PFAM_MATCHES_FILE, prot_dict)
    # Extract elm info.
    elm_names = get_elm_names(ELM_CLASSES_FILE)
    ints_elm = extract_elm_interactions(ELM_INTDOM_FILE, elm_names)
    for elm_ac in elm_names:
        if elm_ac not in ints_elm:
            print "\tNo interaction found for:", elm_ac, elm_names[elm_ac]
    elms = extract_linear_motifs(ELM_HITS_FILE, elm_names, masks, max_overlap=2)

    # Extract PTMs positions from PhosphoSite
    if SP in ["Hsa", "Mmu"]:
        check_file_exists(PSP_FILE)
        ptms = extract_ptms(PSP_FILE, prot_dict, ptms)

    # Get PP interaction from BioGRID file
    bio_ppi, bio_set, bio_id = extract_biogrid_interactions(BIOGRID_FILE,
                                                            prot_dict)
    if SP=="Hsa":
        check_file_exists(HIPPIE_MAP_FILE)
        check_file_exists(HIPPIE_FILE)
        hippie_map = get_HIPPIE_UniProt_map(HIPPIE_MAP_FILE, prot_dict)
        hippie_ppi = extract_HIPPIE_interactions(HIPPIE_FILE, hippie_map)

    ### Create PPI database (only BioGRID now)
    if os.path.isfile(PPI_FILE) and force_ppi==False:
        print_status(PPI_FILE, "exists")
    else:
        create_ppi_database(PPI_FILE, bio_ppi, prot_dict)
        print_status(PPI_FILE, "created")

    ### 4. Make protein data JSON file
    if os.path.isfile(PROT_DATA_FILE) and force_prot_data==False:
        print_status(PROT_DATA_FILE, "exists")
    else:
        print "Creating protein data JSON file"
        create_protein_data_json(prot_dict, alt_ids, pfams, elms,
                               ptms, uni_features, regions, bio_id, bio_set,
                               PROT_DATA_FILE, mode="mongo")
        print_status(PROT_DATA_FILE, "updated")


    ###  5. Calculate Dom/ELM-Dom/ELM Interactions probabilities
    if os.path.isfile(ASSOC_PROB_FILE) and force_probs==False:
        print_status(ASSOC_PROB_FILE, "exists")
    else:
        print "Calculating probabilities"
        ppi = bio_ppi
        if SP=="Hsa":
            ppi = hippie_ppi

        (all_proteins, total_interactions,
            ele_pair_count) = count_protein_and_domain_pairs(ppi, pfams, elms)
        element_count = count_individual_element_frequency(all_proteins, pfams,
                                                           elms)
        print SP, "Total proteins:", len(all_proteins)
        print SP, "Total interactions:", total_interactions

        ints_3did = extract_3did_interactions(EDITED_3DID_FILE)
        ele_names = pfam_names
        ele_names.update(elm_names)
        create_prob_file(ASSOC_PROB_FILE, len(all_proteins), total_interactions,
                         ele_pair_count, element_count, ele_names, ints_3did,
                         ints_elm)
        print_status(ASSOC_PROB_FILE, "created")

    ### 6. Run InterPreTS on PPI
    if os.path.isfile(IPRETS_FILE) and force_iprets==False:
        print_status(IPRETS_FILE, "exists")
    else:
        ppi = bio_ppi
        if SP=="Hsa":
            ppi = hippie_ppi
        create_interprets_file(IPRETS_FILE, ppi, prot_dict, SP_DIR+"iprets/")
        # os.system("gzip "+IPRETS_FILE)
        print_status(IPRETS_FILE, "created")

    print "Finished without problems :-)"
    sys.exit()



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
    # with open_file(IPRETS_FILE) as f: ## lines in this file are already sorted from highest-to-lowest Z-score!!
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
    pdata, ppi, probs, iprets = False, False, False, False
    if "-prot_data" in sys.argv:
        pdata = True
    if "-ppi" in sys.argv:
        ppi = True
    if "-probs" in sys.argv:
        probs = True
    if "-iprets" in sys.argv:
        iprets = True
    main(SP=sys.argv[1], DATA_DIR="static/data/",
         force_prot_data=pdata, force_ppi=ppi, force_probs=probs,
         force_iprets=iprets)
    sys.exit()
