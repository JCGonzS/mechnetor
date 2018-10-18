#!/bin/bash
# From static/data
# Define files paths
sp="Hsa"
biogrid = "BIOGRID-ORGANISM-3.5.165.tab2.txt.gz"
#
#
gunzip species/Hsa/BIOGRID-ORGANISM-3.5.165.tab2.txt.gz
mongoimport --db interactions_Hsa --collection biogrid_Hsa --file species/Hsa/BIOGRID-ORGANISM-3.5.165.tab2.txt --type tsv --headerline
gzip species/Hsa/BIOGRID-ORGANISM-3.5.165.tab2.txt
#
#
gunzip species/Hsa/human_aaa_biogrid_i2.txt.gz
mongoimport --db interactions_Hsa --collection iprets_Hsa --file species/Hsa/human_aaa_biogrid_i2.txt --type tsv --headerline
gzip species/Hsa/human_aaa_biogrid_i2.txt
#
#
gunzip species/Hsa/protein_data_new_Hsa_mongo.json.gz
mongoimport --db protein_data --collection Hsa --file species/Hsa/protein_data_new_Hsa_mongo.json --type json
gzip species/Hsa/protein_data_new_Hsa_mongo.json
#
#
gunzip common/3did_flat_edited-2018_04.tsv.gz
mongoimport --db interactions_common --collection db3did --file common/3did_flat_edited-2018_04.tsv --type tsv --headerline
gzip common/3did_flat_edited-2018_04.tsv
#
#
