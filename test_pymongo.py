#!/usr/bin/env python

### test PyMongo to connect to mongoDB ###

## Data uploaded by:
# mongoimport --db test_3 --collection protein_data --file static/data/protein_data_new_Hsa.json
# from working directory:
# /net/home.isilon/ag-russell/bq_kweise/kyle_piv/kyle_app

## Related links ##
# https://docs.mongodb.com/manual/tutorial/query-documents/
# https://docs.mongodb.com/manual/tutorial/query-embedded-documents/
from pymongo import MongoClient
import datetime, sys
import pprint
import json
import bson
import base64
from bson.binary import Binary

# connection to MongoDB database #
client = MongoClient('localhost', 27017)

# Get database #
db = client['interactions_Hsa']
# Get collection #
data = db['domain_propensities_Hsa']

pfam_a = "XkdN"
pfam_b = "SH3_9"
obs_min=1
lo_min=2.0
ndom_min=1


cursor = data.find({"$or": [{"#DOM1": pfam_a, "DOM2": pfam_b},
					{"#DOM1": pfam_b, "DOM2": pfam_a}],
					"OBS": {"$gte": obs_min}, "LO": {"$gte": lo_min},
					"N_DOM1": {"$gte": ndom_min}, "N_DOM2": {"$gte": ndom_min}},
					{ "_id": 0})

for c in cursor:
	pprint.pprint(c)
	# print c["LO"]

print "done"
sys.exit()
# sys.exit()
# pipeline = [
# 	{"$unwind": "$pfams"},
# 	{"$match": {"uniprot_acc": a }},
# 	# {"$match": {"$and":[{"pfams":{"$ne":[]}} , { "pfams.e-val":{"$gt":0}}]}},
# 	# {"$sort":{"pfams.e-val":1}}
# 	# {"$limit":1}
# 	{"$project": {"_id": 0, "pfams.name": 1 }}
# 	]
# pprint.pprint(list(data.aggregate(pipeline)))
# sys.exit()



examples = 0
if examples == 1:
	# 1. Get one document from collection
	data.find_one()

	# 2. Get all documents from collection
	# for gene in data.find():
	pprint.pprint(gene)

	# 3. Query document for uniprot ID ~#
	#    Query protein_data_collection for documents with uniprot_id = ZCRB1_HUMAN
	cursor = data.find({"uniprot_id": "ZCRB1_HUMAN"})
	for doc in cursor:
		pprint.pprint(doc)

	# 4. Query for a nested value in document ~#
	#    Query protein_data_collection for the e-val field embedded in the domains field, less than 2.2e-50 #
	cursor = protein_data.find({"domains.e-val":{"$lt": 2.2e-75}})
	for doc in cursor:
		pprint.pprint(doc)

	# 5.Find document with the minimum domain.e-val value
	#   Approach 1 - sort and get top 1
	#   needs to be indexed
	cursor = protein_data.find().sort(["domains.e-val", 1]).limit(1)

	# 6. MongoDB aggregation framework
	pipeline = [
		{"$unwind": "$pfams"},
		{"$match": {"$and":[{"pfams":{"$ne":[]}} , { "pfams.e-val":{"$gt":0}}]}},
		{"$sort":{"pfams.e-val":1}},
		{"$limit":1}
		]
	pprint.pprint(list(protein_data.aggregate(pipeline)))

	pprint.pprint(list(protein_data.find({"domains":{"$ne":[]}}).sort({"domains.e-val":1}).limit(1)))


	db.command('aggregate', 'protein_data', pipeline=pipeline, explain=True)
	cursor = prot_data_collection.aggregate(pipeline)

	for doc in cursor:
		pprint.pprint(doc)
