#!/usr/bin/env python

import os, sys, re
import pymongo
from pymongo import MongoClient

client = MongoClient("localhost", 27017)
data_dir = ("/net/home.isilon/ag-russell/bq_jgonzalez/int2mech/piv_app/"+
            "static/data/")

UNI_VERSION     = "Nov2019" # date of download
PFAM_VERSION    = "r32.0"   # release 32
ELM_VERSION     = "Mar2020"#"Oct2019"
BIOGRID_VERSION = "3.5.178" # release
PSP_VERSION     = "Mar2019"

key = sys.argv[1]

if key == "common":
    com_dir = data_dir+"common/"
    print "Creating common databases"

    # # Dom-Dom database
    # data_file = com_dir+"ddi_db.tsv.gz"
    # ddi_data = client["common"]["ddi_db"]
    # if (os.path.isfile(data_file)
    # and ("import" in sys.argv[1:] or ddi_data.count()==0)):
    #     os.system("zcat "+data_file+" | mongoimport -d common"+
    #                 " -c ddi_db --type tsv --headerline --drop")
    #     ddi_data.drop_indexes()
    #     ddi_data.create_index([("Pfam_Acc_A", pymongo.ASCENDING),
    #                             ("Pfam_Acc_B", pymongo.ASCENDING)])

    # ## 3DID DDI
    # db3did_data = client["common"]["db_3did"]
    # data_file = com_dir+"3did_flat_edited-2019_01.tsv.gz"
    # if (os.path.isfile(data_file)
    # and ("import" in sys.argv[1:] or db3did_data.count()==0)):
    #     os.system("zcat "+data_file+" | mongoimport -d common"+
    #               " -c db_3did --type tsv --headerline --drop")
    # db3did_data.drop_indexes()
    # db3did_data.create_index([("Pfam_Acc_A", pymongo.ASCENDING),
    #                           ("Pfam_Acc_B", pymongo.ASCENDING)])

    # ## 3DID DMI
    # dmi_data = client["common"]["dmi_3did"]
    # data_file = com_dir+"3did_dmi_flat_edited-2019_01.tsv.gz"
    # if (os.path.isfile(data_file)
    # and ("import" in sys.argv[1:] or dmi_data.count()==0)):
    #     os.system("zcat "+data_file+" | mongoimport -d common"+
    #               " -c dmi_3did --type tsv --headerline --drop")
    # sys.exit()
    ## ELM-Dom
    elm_int_data = client["common"]["elm_int_dom"]
    data_file = com_dir+"elm_interaction_domains_"+ELM_VERSION+"_reviewed.tsv"
    if (os.path.isfile(data_file)
    and ("import" in sys.argv[1:] or elm_int_data.count()==0)):
        os.system("cat "+data_file+" | mongoimport -d common"+
                  " -c elm_int_dom --type tsv --headerline --drop")
    elm_int_data.drop_indexes()

    # ## ELM classes
    # elm_classes = client["common"]["elm_classes"]
    # data_file = com_dir+"elm_classes_"+ELM_VERSION+".tsv"
    # if (os.path.isfile(data_file)
    # and ("import" in sys.argv[1:] or elm_classes.count()==0)):
    #     os.system("cat "+data_file+" | mongoimport -d common"+
    #               " -c elm_classes --type tsv --headerline --drop")
    # elm_classes.drop_indexes()

    # ## Pfam info
    # pfam_data = client["common"]["pfamA_data"]
    # data_file = com_dir+"Pfam-A.hmm_"+PFAM_VERSION+".tsv.gz"
    # if (os.path.isfile(data_file)
    # and ("import" in sys.argv[1:] or pfam_data.count()==0)):
    #     os.system("zcat "+data_file+" | mongoimport -d common"+
    #               " -c pfamA_data --type tsv --headerline --drop")
    # pfam_data.drop_indexes()

# for sps in ["Xla", "Ath", "Cel", "Dme", "Dre", "Sce"]:#"Hsa", Mmu"]:
elif key in ["Xla", "Ath", "Cel", "Dme", "Dre", "Sce", "Hsa", "Mmu"]:
    sps = key
    sps_dir = data_dir+"species/"+sps+"/"
    print "Creating databases for", sps

    ## Protein Data (ORGANISM-specific)
    # protein_data = client[sps]["protein_data"]
    # data_file = sps_dir+"protein_data_"+sps+"_mongo.json.gz"
    # if (os.path.isfile(data_file)
    # and ("import" in sys.argv[1:] or protein_data.count()==0)):
    #     os.system("zcat "+data_file+" | mongoimport -d "+sps+
    #               " -c protein_data --type json --drop")
    # protein_data.drop_indexes()
    # protein_data.create_index("uni_ac")
    # protein_data.create_index("uni_id")
    # protein_data.create_index("data_class")

    ## Protein Data (COMMON)
    # protein_data = client["common"]["protein_data"]
    # data_file = sps_dir+"protein_data_"+sps+"_mongo.json.gz"
    # if os.path.isfile(data_file):
    #     os.system("zcat "+data_file+" | mongoimport -d common -c protein_data --type json "+
    #                 "--mode upsert --upsertFields uni_id")

    ## PPI database (compendium)
    ppi_data = client[sps]["ppi_db"]
    data_file = sps_dir+"ppi_db_"+sps+".tsv.gz"
    if (os.path.isfile(data_file)
    and ("import" in sys.argv[1:] or ppi_data.count()==0)):
        os.system("zcat "+data_file+" | mongoimport -d "+sps+
                  " -c ppi_db --type tsv --headerline --drop")
    ppi_data.drop_indexes()
    ppi_data.create_index([("Acc_A", pymongo.ASCENDING),
                           ("Acc_B", pymongo.ASCENDING)])


    ## Prot Element Association
    prob_data = client[sps]["association_probabilities"]
    data_file = sps_dir+"prot_ele_association_prob_"+sps+".tsv.gz"
    if (os.path.isfile(data_file)
    and ("import" in sys.argv[1:] or prob_data.count()==0)):
        os.system("zcat "+data_file+" | mongoimport -d "+sps+
                " -c association_probabilities --type tsv --headerline --drop")
    prob_data.drop_indexes()
    prob_data.create_index([("Acc_A", pymongo.ASCENDING),
                            ("Acc_B", pymongo.ASCENDING)])

    # ## InterPreTS
    # iprets_data = client[sps]["interprets"]
    # data_file = sps_dir+"interprets_results_"+sps+".txt.gz"
    # if (os.path.isfile(data_file)
    # and ("import" in sys.argv[1:] or iprets_data.count()==0)):
    #     os.system("zcat "+data_file+" | mongoimport -d "+sps+
    #               " -c interprets --type tsv --headerline --drop")
    # iprets_data.drop_indexes()
    # iprets_data.create_index([("#Gene1", pymongo.ASCENDING),
    #                          ("Gene2", pymongo.ASCENDING)])

    # ## COSMIC
    # if sps == "Hsa":
    #     cosmic_data = client["cosmic_v87"]["genome_screens"]
    #     data_file = sps_dir+"Cosmicv87_GenomeScreens_parsed.tsv.gz"
    #     if (os.path.isfile(data_file)
    #     and ("import" in sys.argv[1:] or cosmic_data.count()==0)):
    #         os.system("zcat "+data_file+" | mongoimport -d cosmic_v87"+
    #                   " -c genome_screens --type tsv --headerline --drop")
    #     cosmic_data.drop_indexes()
    #     cosmic_data.create_index("uni_ac")


# ## COMMON protein data indexes
# protein_data = client["common"]["protein_data"]
# protein_data.drop_indexes()
# protein_data.create_index("uni_ac")
# protein_data.create_index("uni_id")
# protein_data.create_index("data_class")
# protein_data.create_index("organism")
# protein_data.create_index( [("organism", pymongo.ASCENDING),
#                             ("uni_ac", pymongo.ASCENDING)])