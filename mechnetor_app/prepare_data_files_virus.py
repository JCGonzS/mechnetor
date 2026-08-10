#!/usr/bin/python3
import sys, re, os, gzip, math, pprint, datetime, random, string, itertools
import site; site.getsitepackages()
import os.path
from collections import defaultdict
from Bio import SwissProt, SeqIO
import find_all_slims
import prepare_data_files as prep
import pandas as pd
# import run_interprets

def open_file(input_file, mode="rt"):
    """ Open file Zipped or not
    """
    if re.search(".gz$", input_file):
        infile = gzip.open(input_file, mode)
    else:
        infile = open(input_file, mode)
    return infile

def extract_pfam_scov2(pfam_file):
    pfams = defaultdict(lambda: defaultdict(set) )

    with open(pfam_file) as f:
        for line in f:
            if line[0] != "#":
                t = line.rstrip().split()
                uni_ac = t[3].split("|")[1]
                pfam_ac = t[1]
                evalue = float(t[6])
                start, end = int(t[19]), int(t[20])
                pfams[uni_ac][pfam_ac].add((start, end, evalue))
    return pfams

def merge_pfam_files(full_file, matches_file, pfamscan_file):
    with open_file(full_file, "wt") as out:
        cols = ["#uni_ac", "ali_start", "ali_end", "env_start", "env_end", "hmm_ac",
                "hmm_name", "type", "hmm_start", "hmm_end", "hmm_len", "bit_score",
                "e-value", "clan"]
        out.write("\t".join(cols)+"\n")

        with open_file(matches_file, "rt") as f:
            for line in f:
                if line.strip() and line[0]!="#":
                    t = line.rstrip().split()
                    out.write( "\t".join([t[3].split("|")[1],t[17],t[18],t[19],t[20],t[1],t[0],"Domain",t[15],t[16],t[2],t[7],t[6],"-"])+"\n" )

        if os.path.isfile(pfamscan_file):
            with open_file(pfamscan_file, "rt") as f:
                for line in f:
                    if line.strip() and line[0]!="#":
                        t = line.rstrip().split()
                        out.write( "\t".join(t[:5]+[t[5].split(".")[0]]+t[6:-2]+t[-1:])+"\n")
    return

def main( SP="SARS2",
          DATA_DIR="",
          max_overlap=0.2,
          no_overlap_between_doms_and_lms=True,
          force_prot_data=False,
          force_ppi=False,
          force_probs=False,
          force_iprets=False
        ):

    ### DEFINE required data files
    COM_DIR         = DATA_DIR+"common/"
    HSA_DIR         = DATA_DIR+"species/Hsa/"
    SP_DIR          = DATA_DIR+"species/"+SP+"/"
    UNI_VERSION     = "2020-05"
    PFAM_VERSION    = "r33-1"
    ELM_VERSION     = "Mar2020"

    DB3DID_VERSION  = "2020-01"
    BIOGRID_VERSION = "4.2.191" # release
    PSP_VERSION     = "Mar2019"


    UNI_TEXT_FILE     = SP_DIR+"uniprot_"+UNI_VERSION+"_proteome_"+SP+".txt.gz"
    UNI_FEAT_FILE     = SP_DIR+"uniprot_"+UNI_VERSION+"_features_"+SP+".txt.gz"
    ID_MAP_FILE       = SP_DIR+"id2uniprot_mapping_table_"+UNI_VERSION+"_"+SP+".tsv.gz"
    SEQ_FILE          = SP_DIR+"sequences_"+SP+".fasta"
    PTM_SQL_FILE      = SP_DIR+"ptms_SQL_"+SP+".tsv.gz"
    BIOGRID_FILE      = SP_DIR+"BIOGRID-PROJECT-covid19_coronavirus_project-INTERACTIONS-4.3.194.tab3.txt"
    PPI_FILE          = SP_DIR+"ppi_db_"+SP+".tsv.gz"
    PROT_DATA_FILE    = SP_DIR+"protein_data_"+SP+".tsv.gz"
    PFAM_MATCHES_FILE = SP_DIR+"SARSCoV2_pfam_matches.scan.txt"
    PFAMSCAN_FILE     = SP_DIR+"PfamScan_"+PFAM_VERSION+"_"+SP+".txt.gz"
    PFAM_FULL_FILE    = SP_DIR+"Pfam-A_"+PFAM_VERSION+"_matches_full_"+SP+".tsv.gz"
    ELM_HITS_FILE     = SP_DIR+"elm_hits_"+ELM_VERSION+"_"+SP+".tsv.gz"  
    LM_3DID_HITS_FILE = SP_DIR+"3did_slims_"+DB3DID_VERSION+"_"+SP+".tsv.gz"
    LM_HITS_SQL_FILE  = SP_DIR+"all_lm_hits_SQL_"+SP+".tsv.gz"
    IPRETS_REV_FILE   = SP_DIR+"interprets_reviewed_results_"+SP+".txt.gz"
    ASSOC_PROB_FILE   = SP_DIR+"prot_ele_association_prob_"+SP+".tsv.gz"
    NEW_ASSOC_PROB_FILE = SP_DIR+"prot_ele_association_prob_edited_"+SP+".tsv.gz"    

    ## HUMAN
    HUMAN_PROT_DATA = HSA_DIR+"protein_data_Hsa.tsv.gz"

    ## Common files:
    EDITED_PFAM_DAT_FILE = COM_DIR+"Pfam-A_"+PFAM_VERSION+".hmm.tsv.gz"
    ELM_CLASSES_FILE     = COM_DIR+"elm_"+ELM_VERSION+"_classes.tsv"
    ELM_INSTANCES_FILE   = COM_DIR+"elm_"+ELM_VERSION+"_instances.tsv"
    ELM_DOM_EDIT_FILE    = COM_DIR+"elm_"+ELM_VERSION+"_int_domains_final.tsv"
    DMI_3DID_FILE        = COM_DIR+"3did_dmi_flat_"+DB3DID_VERSION+".gz"
    UNI_DISEASE_FILE     = COM_DIR+"diseases-all.tab.gz"
    UNI2BIO_FILE         = COM_DIR+"UniProt_to_BioGRID_"+BIOGRID_VERSION+".tab.txt"
    PDB2UNIPROT_FILE     = COM_DIR+"pdb_chain_uniprot_2020-11.tsv.gz"
    DDI_FILE             = COM_DIR+"DDI_db.tsv.gz"

    pfam_names, pfam_accs, pfam_types = prep.read_pfam_dat(EDITED_PFAM_DAT_FILE)
    elm_map, elm_classes = prep.parse_elm_classes(ELM_CLASSES_FILE)
    ints_elm = prep.extract_elm_interactions(ELM_DOM_EDIT_FILE)    
    edit_lines, lm_3did, all_pdbs = prep.get_dmi_3did(DMI_3DID_FILE, pfam_accs)

    cov2pp = {
        "Host translation inhibitor nsp1": ["nsp1"],
        "Non-structural protein 2": ["nsp2"],
        "Non-structural protein 3": ["PL-PRO", "nsp3"],
        "Non-structural protein 4": ["nsp4"],
        "3C-like proteinase": ["3CL-PRO", "nsp5"],
        "Non-structural protein 6": ["nsp6"],
        "Non-structural protein 7": ["nsp7"],
        "Non-structural protein 8": ["nsp8"],
        "Non-structural protein 9": ["nsp9"],
        "Non-structural protein 10": ["nsp10", "GFL"],
        "Non-structural protein 11": ["nsp11"],
        "RNA-directed RNA polymerase": ["RdRp", "nsp12"],
        "Helicase": ["Hel", "nsp13"],
        "Proofreading exoribonuclease": ["ExoN", "nsp14"],
        "Uridylate-specific endoribonuclease": ["NendoU", "nsp15"],
        "2'-O-methyltransferase": ["nsp16"],
        "Spike protein S1": ["S1"],
        "Spike protein S2": ["S2"],
        "Spike protein S2'": ["S2'"]
    }


    (prot_data, uni_id_map, uni_seqs, 
        modres, uni_feats) = prep.extract_protein_data_from_uniprot_text(UNI_TEXT_FILE)

    for uni_ac in uni_feats:
        if "CHAIN" in uni_feats[uni_ac]:
            for i, chain in enumerate(uni_feats[uni_ac]["CHAIN"]):
                i+=1
                uni_id = uni_id_map[uni_ac][0]
                start, end, name = chain[0], chain[1], chain[3]
                seq = uni_seqs[uni_ac][start-1:end]
                if uni_id=="SPIKE_SARS2" and start==13 and end==1273:
                    continue
                new_uni_id = uni_id+":"+str(start)+"-"+str(end)
                new_uni_ac = uni_ac+":"+str(start)+"-"+str(end)
                prot_data[new_uni_id] = {
                    "gn": cov2pp[name][0],
                    "ac": new_uni_ac,
                    "des": name,
                    "dc": prot_data[uni_id]["dc"],
                    "seq": seq,
                    "length": len(seq),
                    "biogrid_id": set(),
                    "sorted_ints": []
                }
                uni_seqs[new_uni_ac] = seq
                for x in [new_uni_id, new_uni_ac]+cov2pp[name]:
                    uni_id_map[x.upper()].append(new_uni_id)

    org_map = {}
    for uni_id in prot_data:
        org_map[uni_id] = SP
        if uni_id=="R1A_SARS2":
            prot_data[uni_id]["gn"] = "1a"
        elif uni_id=="ORF3B_SARS2":
            prot_data[uni_id]["gn"] = "3b"
        elif uni_id=="ORF3C_SARS2":
            prot_data[uni_id]["gn"] = "3c"    
        elif uni_id=="ORF3D_SARS2":
            prot_data[uni_id]["gn"] = "3d"

    with prep.open_file(SEQ_FILE, "w") as out:
            for uni in uni_seqs:
                out.write(">"+uni+"\n")
                out.write(uni_seqs[uni]+"\n")

    prep.print_status(SEQ_FILE, "created")

    prep.create_id_map_table(ID_MAP_FILE, uni_id_map, SP)
    prep.print_status(ID_MAP_FILE, "created")

    uni_dis_names = prep.get_uniprot_diseases(UNI_DISEASE_FILE)
    prep.create_uni_features_table(UNI_FEAT_FILE, uni_feats, uni_dis_names)
    prep.print_status(UNI_FEAT_FILE, "created")

    prep.create_ptms_table(PTM_SQL_FILE, modres, uni_seqs)
    prep.print_status(PTM_SQL_FILE, "created")

    prot_data = prep.add_uniprot_biogrid_mapping(UNI2BIO_FILE, uni_id_map, prot_data)
    bio2uni = prep.get_bio_2_uni_mapping(prot_data)

    prot_data2 = prot_data.copy()

    with prep.open_file(HUMAN_PROT_DATA) as f:
        for i, line in enumerate(f):
            if i > 0:
                t = line.rstrip().split("\t")
                uni_id = t[0]
                prot_data2[uni_id] = {}
                prot_data2[uni_id]["ac"] = t[1]
                prot_data2[uni_id]["seq"] = t[7]
                prot_data2[uni_id]["length"] = len(t[7])
                org_map[uni_id] = "Hsa"
                try:
                    bio_ids = t[8]
                    for bio_id in bio_ids.split(", "):
                        bio2uni[bio_id].add(uni_id)
                except:
                    pass               

    # Get PP interaction from BioGRID file
    biogrid_ppi, biogrid_ppi_all, interaction_info = prep.extract_biogrid_interactions(BIOGRID_FILE, bio2uni)    

    prep.create_ppi_database(PPI_FILE, biogrid_ppi, interaction_info)
    prep.print_status(PPI_FILE, "created")

    prot_data = prep.add_sorted_interactors(prot_data, biogrid_ppi_all, bio2uni)
    prep.create_protein_data_table(PROT_DATA_FILE, prot_data, SP)
    prep.print_status(PROT_DATA_FILE, "created")
    
    # Check for PfamScan
    if not os.path.isfile(PFAMSCAN_FILE):
        pfam_matches = extract_pfam_scov2(PFAM_MATCHES_FILE)
        prep.run_pfamscan(PFAMSCAN_FILE, uni_seqs, pfam_matches)
        prep.print_status(PFAMSCAN_FILE, "created")

    # Check for Complete Pfam file (matches+pfamscan)
    merge_pfam_files(PFAM_FULL_FILE, PFAM_MATCHES_FILE, PFAMSCAN_FILE)
    prep.print_status(PFAM_FULL_FILE, "created")

    pfam_matches = prep.extract_pfam_doms(PFAM_FULL_FILE)

    # Check for ELM hits
    if not os.path.isfile(ELM_HITS_FILE):
        find_all_slims.main(SEQ_FILE, elm_classes,
                            print_out=True, outfile=ELM_HITS_FILE, tmp_file="/tmp/lm_search_"+SP+".txt")
        prep.print_status(ELM_HITS_FILE, "created")

    # Check for 3did-motifs hits
    if not os.path.isfile(LM_3DID_HITS_FILE):
        find_all_slims.main(SEQ_FILE, lm_3did,
                            print_out=True, outfile=LM_3DID_HITS_FILE, tmp_file="/tmp/lm_search_"+SP+".txt")
        prep.print_status(LM_3DID_HITS_FILE, "created")

    elm_tp, elm_fp = prep.extract_elm_instances(ELM_INSTANCES_FILE, elm_map, uni_seqs.keys())
    elm_hits = prep.extract_slim_hits(ELM_HITS_FILE, elm_tp, elm_fp)

    # Extract 3DID_LM hits
    pdb2uni = prep.get_pdb2uniprot_map(PDB2UNIPROT_FILE, uni_seqs.keys())
    tp = prep.get_3did_dmi_instances(lm_3did, pdb2uni, uni_seqs)
    lm3d_hits = prep.extract_slim_hits(LM_3DID_HITS_FILE, tp, {})

    if not os.path.isfile(LM_HITS_SQL_FILE):
        prep.reformat_slim_hit_file(LM_HITS_SQL_FILE , elm_hits, lm3d_hits)
        prep.print_status(LM_HITS_SQL_FILE , "created")
    
    sys.exit()

    ppi = defaultdict(set)
    for a in biogrid_ppi:
        for uni_a in bio2uni[a]:
            for b in biogrid_ppi[a]:
                for uni_b in bio2uni[b]:
                    if uni_a < uni_b:
                        ppi[uni_a].add(uni_b)
                    else:
                        ppi[uni_b].add(uni_a)


    # if prep.check_file_exists(ASSOC_PROB_FILE):
    #     prep.edit_assoc_file(ASSOC_PROB_FILE, NEW_ASSOC_PROB_FILE, org)
    # else:
    #     prep.check_file_exists(ASSOC_PROB_FILE)

    # sys.exit()
    
    ###  5. Calculate Interactions Probabilities between Protein Elements
    if os.path.isfile(ASSOC_PROB_FILE):
        prep.print_status(ASSOC_PROB_FILE, "exists")
    else:
        ele_names = prep.merge_protein_elements_names(pfam_names, elm_map, lm_3did)
        eles = prep.merge_protein_elements(pfam_matches, elm_hits, lm3d_hits)

        # (all_proteins, total_interactions,
        #     ele_pair_count) = count_protein_and_domain_pairs(ppi, eles)

    sys.exit()
    #     element_count = count_individual_element_frequency(all_proteins, eles)

    #     ints_ddi = extract_ddi_interactions(DDI_FILE)
    #     ints_dmi = extract_3did_dmi_interactions(EDITED_DMI_3DID_FILE)
    #     print SP, "Total proteins:", len(all_proteins)
    #     print SP, "Total interactions:", total_interactions
    #     create_prob_file(ASSOC_PROB_FILE, len(all_proteins), total_interactions,
    #                      ele_pair_count, element_count, ele_names, ints_ddi,
    #                      ints_elm, ints_dmi)
    #     print_status(ASSOC_PROB_FILE, "created")


    ### Run InterPreTS on PPI
    if os.path.isfile(IPRETS_REV_FILE):
        prep.print_status(IPRETS_REV_FILE, "exists")
    else:
        prep.create_interprets_file(IPRETS_REV_FILE, ppi, prot_data2, org_map,
                                    SP_DIR+"iprets/", DATA_DIR+"species/",
                                    "/net/home.isilon/ds-russell/blastdb/pdbseq")
        prep.print_status(IPRETS_REV_FILE, "created")
    sys.exit()
  


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
    try:
        sp = sys.argv[1]
    except:
        sp = None
    main(SP=sp, DATA_DIR="/net/home.isilon/ag-russell/bq_jgonzalez/int2mech/data/",
         force_prot_data=pdata, force_ppi=ppi, force_probs=probs,
         force_iprets=iprets)
    sys.exit()
