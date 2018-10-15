#!/usr/bin/env python

import sys, re
import gzip, json
import pandas as pd


## Data Files
data_dir = "static/data/"
gen_data_dir = data_dir+"general/"
species = "Hsa"
sp_data_dir = data_dir+"species/"+species+"/"
biogrid_file = sp_data_dir+"BIOGRID-ORGANISM-species-3.4.156.tab.txt.gz"
biogrid_js = sp_data_dir+"BIOGRID-ORGANISM-species-3.4.15.json.gz"

def open_file(input_file, mode="r"):
    """ Open file Zipped or not
    """
    if re.search(".gz$", input_file):
        infile = gzip.open(input_file, mode)
    else:
        infile = open(input_file, mode)
    return infile


### BioGRID database
# 1. Read data file into a dataframe
df = pd.read_csv(biogrid_file, compression="gzip", sep="\t")

# 2. Remove columns that will not be needed
remove = ["INTERACTOR_A", "INTERACTOR_B", "ORGANISM_A_ID", "ORGANISM_B_ID"]
df = df.drop(remove, axis=1)

# 3. Save as JSON
with open_file(biogrid_js, "w") as out:
    json.dump(df.to_json(orient="index"), out)
