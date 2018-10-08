#!/usr/bin/env python
### test PyMongo to connect to mongoDB ###

import sys,json

json_file = 'data/species/Hsa/protein_data_new_Hsa_normal.json'

protein_data = json.load(open(json_file))
print protein_data["P04637"]

sys.exit()
print(json.dumps(file_data, sort_keys=True,indent=4))
