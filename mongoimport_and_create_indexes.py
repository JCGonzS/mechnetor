#!/usr/bin/env python

import os, sys, re
import pymongo
from pymongo import MongoClient

client = MongoClient("localhost", 27017)

## Import DBs
data_dir = "/net/home.isilon/ag-russell/bq_jgonzalez/int2mech/piv_app/static/data/"
com_dir = data_dir+"common/"

for sps in ["Hsa", "Dme", "Mmu"]:
    sps_dir = data_dir+"species/"+sps+"/"
    print "Creating databases for", sps
    ## Protein Data
    protein_data = client[sps]["protein_data"]
    data_file = sps_dir+"protein_data_"+sps+"_mongo.json.gz"
    if (os.path.isfile(data_file) and ("import" in sys.argv[1:] or protein_data.count()==0)):
        os.system("zcat "+data_file+" | mongoimport -d "+sps+" -c protein_data --type json --drop")
    protein_data.drop_indexes()
    protein_data.create_index("uni_ac")
    protein_data.create_index("uni_id")
    protein_data.create_index("data_class")
    sys.exit()
    ## BioGRID
    biogrid_data = client[sps]["biogrid"]
    data_file = sps_dir+"BIOGRID-ORGANISM-"+sps+"-3.5.172.tab2.txt.gz"
    if (os.path.isfile(data_file) and ("import" in sys.argv[1:] or biogrid_data.count()==0)):
        os.system("zcat "+data_file+" | mongoimport -d "+sps+" -c biogrid --type tsv --headerline --drop")
    biogrid_data.drop_indexes()
    biogrid_data.create_index([("Official Symbol Interactor A", pymongo.ASCENDING),("Official Symbol Interactor B", pymongo.ASCENDING)])

    ## Dom-Dom Association
    dd_ass_data = client[sps]["dom_dom_ass"]
    data_file = sps_dir+"dom_dom_association_"+sps+".tsv.gz"
    if (os.path.isfile(data_file) and ("import" in sys.argv[1:] or dd_ass_data.count()==0)):
        os.system("zcat "+data_file+" | mongoimport -d "+sps+" -c dom_dom_ass --type tsv --headerline --drop")
    dd_ass_data.drop_indexes()
    dd_ass_data.create_index([("dom_ac_a", pymongo.ASCENDING),("dom_ac_b", pymongo.ASCENDING)])

    ## InterPreTS
    iprets_data = client[sps]["interprets_biogrid"]
    data_file = sps_dir+"i2_biogrid_"+sps+".txt.gz"
    if (os.path.isfile(data_file) and ("import" in sys.argv[1:] or iprets_data.count()==0)):
        os.system("zcat "+data_file+" | mongoimport -d "+sps+" -c interprets_biogrid --type tsv --headerline --drop")
        iprets_data.drop_indexes()
        iprets_data.create_index([("#Gene1", pymongo.ASCENDING),("Gene2", pymongo.ASCENDING)])

    if sps == "Hsa":
        ## COSMIC
        cosmic_data = client["cosmic_v87"]["genome_screens"]
        data_file = sps_dir+"Cosmicv87_GenomeScreens_parsed.tsv.gz"
        if (os.path.isfile(data_file) and ("import" in sys.argv[1:] or cosmic_data.count()==0)):
            os.system("zcat "+data_file+" | mongoimport -d cosmic_v87 -c genome_screens --type tsv --headerline --drop")
        cosmic_data.drop_indexes()
        cosmic_data.create_index("uni_ac")

print "Creating common databases"
## 3did
db3did_data = client["common"]["3did"]
data_file = com_dir+"3did_flat_edited-2018_04.tsv.gz"
if (os.path.isfile(data_file) and ("import" in sys.argv[1:] or db3did_data.count()==0)):
    os.system("zcat "+data_file+" | mongoimport -d common -c 3did --type tsv --headerline --drop")
db3did_data.drop_indexes()
db3did_data.create_index([("Pfam_Acc_A", pymongo.ASCENDING),("Pfam_Acc_B", pymongo.ASCENDING)])

## ELM-Dom
elm_int_data = client["common"]["elm_int_dom"]
data_file = com_dir+"elm_interaction_domains_edited_Jan18.tsv"
if (os.path.isfile(data_file) and ("import" in sys.argv[1:] or elm_int_data.count()==0)):
    os.system("cat "+data_file+" | mongoimport -d common -c elm_int_dom --type tsv --headerline --drop")
elm_int_data.drop_indexes()

## ELM classes
elm_classes = client["common"]["elm_classes"]
data_file = com_dir+"elm_classes_May2019.tsv.gz"
if (os.path.isfile(data_file) and ("import" in sys.argv[1:] or elm_classes.count()==0)):
    os.system("zcat "+data_file+" | mongoimport -d common -c elm_classes --type tsv --headerline --drop")
elm_classes.drop_indexes()

## Pfam info
pfam_data = client["common"]["pfamA_data"]
data_file = com_dir+"Pfam-A.hmm_r32.0.tsv.gz"
if (os.path.isfile(data_file) and ("import" in sys.argv[1:] or pfam_data.count()==0)):
    os.system("zcat "+data_file+" | mongoimport -d common -c pfamA_data --type tsv --headerline --drop")
pfam_data.drop_indexes()
# no indexes needed
