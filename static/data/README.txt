-----------------------------
DATA FILES FOR INT2MECH & PIV
-----------------------------
Oct 2018
-----------------------------



Short instructions summary
--------------------------
1. Download all required sources files. See list below

2. Run "prepare_data_files.py" to generate some additional files.

3. Run "create_mongo_dbs.sh" to create Mongo databases and collections

4. make protein data?


Data Files
----------
Data files are divided into:
- species-specific (eg. Homo sapiens or 'Hsa'): --> static/data/species/Hsa/
  - common (same for any species) --> static/data/common/

Original data files must be obtained from their original sources and put into
one these two directories.
Required source files are the following:
> Species:
  - BioGRID:
  - InterPreTs results (only BioGRID interacting pairs):

> Common
  - 3did:



The following instructions explains in detail how this files need to be processed
so they can be used by PIV.


---------------------------
Importing files to MongoDB
---------------------------

MongoDB accepts files in json, csv or tsv formats.

********************
A. SPECIES-SPECIFIC
********************

I. BioGRID -  Biological General Repository for Interaction Datasets
--------------------------------------------------------------------
* thebiogrid.org/
* version 3.5.165 (Oct 2018)
* Source file: "species/sp/BIOGRID-ORGANISM-'SpeciesName'-3.5.165.tab2.txt.gz"
  (Columns definitions at: wiki.thebiogrid.org/doku.php/biogrid_tab_version_2.0)

1. Use 'create_mongo_dbs.sh' to create collection in the Mongo interactions database
2. Create indexes for this database in mongo shell:
    > use interactions_Hsa
    > db.biogrid_Hsa.createIndex( {"Official Symbol Interactor A": 1,
                                   "Official Symbol Interactor B": 1} )


II. InterPreTS - Interaction Prediction through Tertiary Structure
------------------------------------------------------------------
* russelllab.org/cgi-bin/tools/interprets.pl/
* Results of the predictions for all interaction pairs according to BioGRID
* Source file: "species/sp/human_aaa_biogrid_i2.txt.gz"

1. Use 'create_mongo_dbs.sh' to create collection in the Mongo interactions database
2. Create indexes for this database in mongo shell:
    > use interactions_Hsa
    > db.iprets_Hsa.createIndex( { "#Gene1": 1, "Gene2": 1 } )


*********
B. COMMON
*********

II. 3did - Database of three-dimensional interacting domains
------------------------------------------------------------
* 3did.irbbarcelona.org/
* version 2018_04
* Source file: "common/3did_flat-2018-04.gz"

1. File needs to be parsed into a tsv format with "prepare_data_files.py"
   Generates "3did_flat_edited-2018_04.tsv.gz"
2. Use 'create_mongo_dbs.sh' to create collection in the Mongo interactions database
3. Create indexes for this database in mongo shell:
   > use interactions_gen
   > db.db3did.createIndex( { "Pfam_Name_A": 1, "Pfam_Name_B": 1} )




"elm_classes.tsv"
	Source: elm.eu.org
	Contains all the classes of LM in the database with their regular expressions (RE)
	and probabilities.

	This file is used by "find_all_elms.py", which uses the REs to find all ELMs
	present in a provided fasta file. The results is the file "elm_hits.tsv", which
	is required by int2mech.py

"elm_instances.tsv"
	Source: elm.eu.org
	Contains all annotated ELM instances.

	Not used.

####
