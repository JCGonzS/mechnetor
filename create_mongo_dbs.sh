#!/bin/bash
# From static/data
# Define files paths
sp="Hsa"
biogrid = "BIOGRID-ORGANISM-3.5.165.tab2.txt.gz"
gunzip species/Hsa/BIOGRID-ORGANISM-3.5.165.tab2.txt.gz
mongoimport --db interactions_Hsa --collections biogrid_Hsa --file species/Hsa/BIOGRID-ORGANISM-3.5.165.tab2.txt --type tsv --headerline
gunzip species/Hsa/BIOGRID-ORGANISM-3.5.165.tab2.txt
