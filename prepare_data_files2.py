#!/usr/bin/env python

import sys, re, os, gzip, math, pprint, datetime, random, string, itertools
import find_all_slims
import os.path
from collections import defaultdict
# import run_interprets
from Bio import SwissProt

def open_file(input_file, mode="r"):
    """ Open file Zipped or not
    """
    if re.search(".gz$", input_file):
        infile = gzip.open(input_file, mode)
    else:
        infile = open(input_file, mode)
    return infile


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
            print("WARNING: UniProt entries not sorted Reviewed->Unreviewed")

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
                # if feat[0]=="MUTAGEN":
                #     if start!=end:
                #         print uni_ac, main_gene, feat
            elif feat[0]=="REGION":
                region[uni_ac].append(
                    {"start": start,
                     "end": end,
                     "info": info
                    }
                )

            # else:
            #     print feat

        # for ref in record.cross_references:
        #     if "PDB" in ref:
        #         chains = set()
        #         for x in ref[-1].split(", "):
        #             for y in x.split("=")[0].split("/"):
        #                 chains.add(y)
        #         for c in chains:
        #             pdb2uni["pdb|"+ref[1]+"|"+c].add(uni_ac)
    return D, alt_ids, masks, modres, uni_features, region


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
    COM_DIR         = DATA_DIR+"common/"
    SP_DIR          = DATA_DIR+"species/"+SP+"/"
    UNI_VERSION     = "Oct2020"
    PFAM_VERSION    = "r32.0"   # release 32
    ELM_VERSION     = "Mar2020" #"Oct2019"
    BIOGRID_VERSION = "3.5.178" # release
    PSP_VERSION     = "Mar2019"
    DB3DID_VERSION  = "2019_01"

    ## Common files:
    ### PFAM
    PFAM_DAT_FILE        = COM_DIR+"Pfam-A.hmm_"+PFAM_VERSION+".dat.gz"
    EDITED_PFAM_DAT_FILE = COM_DIR+"Pfam-A.hmm_"+PFAM_VERSION+".tsv.gz"
    PFAM_INT_FILE        = COM_DIR+"PfamA_interactions_"+PFAM_VERSION+".txt.gz"

    ### 3did
    FLAT_3DID_FILE       = COM_DIR+"3did_flat-"+DB3DID_VERSION+".gz"
    DMI_3DID_FILE        = COM_DIR+"3did_dmi_flat-"+DB3DID_VERSION+".gz"
    EDITED_DMI_3DID_FILE = COM_DIR+"3did_dmi_flat_edited-"+DB3DID_VERSION+".tsv.gz"

    ### DDI database = Pfam + 3did DDI's
    DDI_FILE             = COM_DIR+"DDI_db.tsv.gz"

    ### ELM
    ELM_CLASSES_FILE     = COM_DIR+"elm_classes_"+ELM_VERSION+".tsv"
    ELM_INSTANCES_FILE   = COM_DIR+"elm_instances_"+ELM_VERSION+".tsv"
    ELM_INTDOM_FILE      = (COM_DIR+"elm_interaction_domains_"
                            +ELM_VERSION+"_reviewed.tsv")

    PDB2UNIPROT_FILE     = COM_DIR+"pdb_chain_uniprot.tsv.gz"

    ## Species files:
    UNI_TEXT_FILE     = SP_DIR+"uniprot_"+UNI_VERSION+"_proteome_"+SP+".txt.gz"
    UNI_FASTA_FILE    = (SP_DIR+"uniprot_"+UNI_VERSION+"_proteome_"+SP
                        +".fasta.gz")
    PFAM_MATCHES_FILE = SP_DIR+"pfamA_"+PFAM_VERSION+"_matches_"+SP+".tsv.gz"
    LM_3DID_HITS_FILE = SP_DIR+"3did_slims_"+DB3DID_VERSION+"_"+SP+".tsv.gz"
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
    PFAMSCAN_FILE     = SP_DIR+"pfamscan_results_"+SP+".txt.gz"

    COV2_PFAM_FILE       = COM_DIR+"Pfam-A.SARS-CoV-2.full.gz"
    COV2_PFAM_SCAN_FILE  = COM_DIR+"SARSCoV2_pfam_matches.scan.txt"

    # PfamScan
    PFAMSCAN          = "pfamscan.py"
    TMP_DIR           = "tmp/"
    TMP_FASTA         = TMP_DIR+"no_pfam_seqs"
    TMP_PFAMOUT       = TMP_DIR+"pfamout"
    evalue = "0.001"
    email = "juan-carlos.gonzalez@bioquant.uni-heidelberg.de"

    # ### 2. Check existence/make required data files
    # ## Check/Edit common files:

    # # Pfam.dat
    # if os.path.isfile(EDITED_PFAM_DAT_FILE):
    #     print_status(EDITED_PFAM_DAT_FILE, "exists")
    # elif check_file_exists(PFAM_DAT_FILE):
    #     edit_pfam_dat(PFAM_DAT_FILE, COV2_PFAM_FILE, EDITED_PFAM_DAT_FILE)
    #     print_status(EDITED_PFAM_DAT_FILE, "created")
    # pfam_names, pfam_accs, pfam_types = read_pfam_dat(EDITED_PFAM_DAT_FILE)
     
    # # 3did flat
    # if os.path.isfile(DDI_FILE):
    #     print_status(DDI_FILE, "exists")
    # elif check_file_exists(FLAT_3DID_FILE):
    #     create_ddi_database(pfam_names, PFAM_INT_FILE, FLAT_3DID_FILE, DDI_FILE)
    #     print_status(DDI_FILE, "created")

    # # elm classes & domain-ints
    # check_file_exists(ELM_CLASSES_FILE)
    # elm_map, elm_classes = parse_elm_classes(ELM_CLASSES_FILE)
    # check_file_exists(ELM_INTDOM_FILE)
    # check_file_exists(ELM_INSTANCES_FILE)

    # # 3did DMI
    # check_file_exists(PDB2UNIPROT_FILE)
    # check_file_exists(DMI_3DID_FILE)
    # edit_lines, lm_3did, all_pdbs = get_dmi_3did(DMI_3DID_FILE, pfam_accs)

    # if os.path.isfile(EDITED_DMI_3DID_FILE):
    #     print_status(EDITED_DMI_3DID_FILE, "exists")
    # else:
    #     with open_file(EDITED_DMI_3DID_FILE, "w") as out:
    #         out.write("\t".join(["MOTIF","REGEX","DOMAIN", "DOMAIN_NAME", "PDBS"])+"\n" )
    #         for line in edit_lines:
    #             out.write(line+"\n")
    #     print_status(EDITED_DMI_3DID_FILE, "created")

    # ## Check/Edit species files:
    # for f in [UNI_TEXT_FILE, UNI_FASTA_FILE, PFAM_MATCHES_FILE, BIOGRID_FILE]:
    #     check_file_exists(f)
    # # Check/Generate ELM annotation
    # if os.path.isfile(ELM_HITS_FILE):
    #     print_status(ELM_HITS_FILE, "exists")
    # else:
    #     find_all_slims.main(UNI_FASTA_FILE, elm_classes,
    #                         print_out=True, outfile=ELM_HITS_FILE, tmp_file="/tmp/lm_search_"+SP+".txt")
    #     print_status(ELM_HITS_FILE, "created")
    
    # if os.path.isfile(LM_3DID_HITS_FILE):
    #     print_status(LM_3DID_HITS_FILE, "exists")
    # else:
    #     find_all_slims.main(UNI_FASTA_FILE, lm_3did,
    #                         print_out=True, outfile=LM_3DID_HITS_FILE, tmp_file="/tmp/lm_search_"+SP+".txt")
    #     print_status(LM_3DID_HITS_FILE, "created")

    ### 3. Get various data
    # Extract protein ID-dictionary, sequences and masks.
    (prot_dict, alt_ids, masks, ptms, uni_features,
        regions) = extract_protein_data_from_uniprot_text(UNI_TEXT_FILE)
    sys.exit()

    # # Extract ELM hits.
    # ints_elm = extract_elm_interactions(ELM_INTDOM_FILE, elm_map)
    # for elm_ac in elm_classes:
    #     if elm_ac not in ints_elm:
    #         print("\tNo interaction found for:", elm_ac, elm_map[elm_ac])
    # elm_tp, elm_fp = extract_elm_instances(ELM_INSTANCES_FILE, elm_map, prot_dict["seq"])
    # elm_hits = extract_slims(ELM_HITS_FILE, elm_tp, elm_fp)
    
    # # Extract 3DID_LM hits
    # pdb2uni = get_pdb2uniprot_map(PDB2UNIPROT_FILE, set(prot_dict["seq"].keys()))
    # tp = get_3did_dmi_instances(lm_3did, pdb2uni, prot_dict["seq"])
    # lm3d_hits = extract_slims(LM_3DID_HITS_FILE, tp, {})
    # # Extract Pfam domains.
    # pfams, pfam_names = extract_pfam_doms(PFAM_MATCHES_FILE, prot_dict)
    # # Calculate missing ones
    # if os.path.isfile(PFAMSCAN_FILE):
    #     print_status(PFAMSCAN_FILE, "exists")
    # else:
    #     if not os.path.exists(TMP_DIR):
    #         os.mkdir(TMP_DIR)
    #     n, i = 0, 0
    #     for uni_ac in prot_dict["seq"]:
    #         if uni_ac not in pfams:
    #             i += 1
    #             if i == 100:
    #                 n += 1
    #                 i = 0
    #             fasta = TMP_FASTA+str(n)+".fasta"
    #             with open_file(fasta, "a") as out:
    #                 out.write(">"+uni_ac+"\n")
    #                 out.write(prot_dict["seq"][uni_ac]+"\n")

    #     for x in range(n+1):
    #         print("Running PfamScan #"+str(x))
    #         fasta = TMP_FASTA+str(x)+".fasta"
    #         pfamout = TMP_PFAMOUT+str(x)
    #         if not os.path.isfile(TMP_PFAMOUT+str(x)+".out.txt"):
    #             cmd = ("python {} --sequence {} --database pfam-a --evalue {}".format(
    #                     PFAMSCAN, fasta, evalue)+
    #                     " --format txt --outfile {} --email {} --quiet".format(
    #                     pfamout, email))
    #             os.system(cmd)

    #     miss = []
    #     with open_file(PFAMSCAN_FILE, "w") as f:
    #         for x in range(n+1):
    #             pfamout_file = TMP_PFAMOUT+str(x)+".out.txt"
    #             if os.path.isfile(pfamout_file):
    #                 with open_file(pfamout_file) as f2:
    #                     for line in f2:
    #                         if line[0]!="#":
    #                             f.write(line)
    #                 os.unlink(pfamout_file)
    #                 os.unlink(TMP_PFAMOUT+str(x)+".sequence.txt")
    #             else:
    #                 miss.append(x)
    #             os.unlink(TMP_FASTA+str(x)+".fasta")

    #     print("Missing Pfamscan for files", miss)
    #     print_status(PFAMSCAN_FILE+".out.txt", "created")

    # pfams, pfam_names = extract_pfamscan(PFAMSCAN_FILE, pfams, pfam_names)
  
    # # Extract PTMs positions from PhosphoSite
    # if SP in ["Hsa", "Mmu"]:
    #     check_file_exists(PSP_FILE)
    #     ptms = extract_ptms(PSP_FILE, prot_dict, ptms)

    # # Get PP interaction from BioGRID file
    # bio_ppi, bio_ppi_all, bio_set, bio_id = extract_biogrid_interactions(BIOGRID_FILE,
    #                                                         prot_dict)
    # if SP=="Hsa":
    #     check_file_exists(HIPPIE_MAP_FILE)
    #     check_file_exists(HIPPIE_FILE)
    #     hippie_map = get_HIPPIE_UniProt_map(HIPPIE_MAP_FILE, prot_dict)
    #     hippie_ppi = extract_HIPPIE_interactions(HIPPIE_FILE, hippie_map)

    # ### Create PPI database (only BioGRID now)
    # if os.path.isfile(PPI_FILE) and force_ppi==False:
    #     print_status(PPI_FILE, "exists")
    # else:
    #     create_ppi_database(PPI_FILE, bio_ppi_all, prot_dict)
    #     print_status(PPI_FILE, "created")

    # ### 4. Make protein data JSON file
    # if os.path.isfile(PROT_DATA_FILE) and force_prot_data==False:
    #     print_status(PROT_DATA_FILE, "exists")
    # else:
    #     print("Creating protein data JSON file")
    #     create_protein_data_json(prot_dict, alt_ids, pfams, elm_hits, lm3d_hits,
    #                            ptms, uni_features, regions, bio_id, bio_set,
    #                            SP, PROT_DATA_FILE, mode="mongo")
    #     print_status(PROT_DATA_FILE, "updated")
    

    # ###  5. Calculate Interactions Probabilities between Protein Elements
    # if os.path.isfile(ASSOC_PROB_FILE) and force_probs==False:
    #     print_status(ASSOC_PROB_FILE, "exists")
    # else:
    #     print("Calculating probabilities")
    #     ppi = bio_ppi
    #     if SP=="Hsa":
    #         ppi = hippie_ppi

    #     eles = merge_protein_elements(pfams, elm_hits, lm3d_hits)
    #     ele_names = merge_protein_elements_names(pfam_names, elm_map, lm_3did)

    #     (all_proteins, total_interactions,
    #         ele_pair_count) = count_protein_and_domain_pairs(ppi, eles)
    #     element_count = count_individual_element_frequency(all_proteins, eles)

    #     ints_ddi = extract_ddi_interactions(DDI_FILE)
    #     ints_dmi = extract_3did_dmi_interactions(EDITED_DMI_3DID_FILE)
    #     print(SP, "Total proteins:", len(all_proteins))
    #     print(SP, "Total interactions:", total_interactions)
    #     create_prob_file(ASSOC_PROB_FILE, len(all_proteins), total_interactions,
    #                      ele_pair_count, element_count, ele_names, ints_ddi,
    #                      ints_elm, ints_dmi)
    #     print_status(ASSOC_PROB_FILE, "created")

    # ### 6. Run InterPreTS on PPI
    # if os.path.isfile(IPRETS_FILE) and force_iprets==False:
    #     print_status(IPRETS_FILE, "exists")
    # else:
    #     ppi = bio_ppi
    #     if SP=="Hsa":
    #         ppi = hippie_ppi
    #     create_interprets_file(IPRETS_FILE, ppi, prot_dict, SP_DIR+"iprets/")
    #     # os.system("gzip "+IPRETS_FILE)
    #     print_status(IPRETS_FILE, "created")

    # print("Finished without problems :-)")
    # sys.exit()



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
