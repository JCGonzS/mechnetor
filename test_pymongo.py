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
db = client['test_3']
# Get collection #
protein_data = db['protein_data']

#~ get one document from collection 
pprint.pprint(protein_data.find_one("P04637"))

#~ get all documents from collection 
#for gene in protein_data.find():
#	pprint.pprint(gene)	

#~ query document for uniprot ID ~#
# Query protein_data_collection for documents with uniprot_id = ZCRB1_HUMAN

cursor = protein_data.find({"uniprot_id": "ZCRB1_HUMAN"})
for doc in cursor:
	pprint.pprint(doc)
sys.exit()
#~ query for a nested value in document ~#
# Query protein_data_collection for the e-val field embedded in the domains field, less than 2.2e-50 #

#~ cursor = protein_data.find({"domains.e-val":{"$lt": 2.2e-75}})
#~ for doc in cursor:
	#~ pprint.pprint(doc)



#~ Find document with the minimum domain.e-val value

#~ Approach 1 - sort and get top 1
#~ needs to be indexed 
#~ cursor = protein_data.find().sort(["domains.e-val", 1]).limit(1)


#~ MongoDB aggregation framework
pipeline = [
	{"$unwind": "$domains"},
	{"$match": {"$and":[{"domains":{"$ne":[]}} , { "domains.e-val":{"$gt":0}}]}},
	{"$sort":{"domains.e-val":1}},
	{"$limit":1}
	]

pprint.pprint(list(protein_data.aggregate(pipeline)))

#~ pprint.pprint(list(protein_data.find({"domains":{"$ne":[]}}).sort({"domains.e-val":1}).limit(1)))


#~ db.command('aggregate', 'protein_data', pipeline=pipeline, explain=True)
#~ cursor = prot_data_collection.aggregate(pipeline)

#~ for doc in cursor:
	#~ pprint.pprint(doc)