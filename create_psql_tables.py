import sys, os, re, gzip
from collections import defaultdict
import psycopg2

def db_exists(cursor, table_name):
    cursor.execute("select exists(select * from information_schema.tables where table_name=\'"+table_name+"\')") 
    return cursor.fetchone()[0]

def create_psql_database(conn, cursor, table_name, columns, drop=False, index=None):
    exists = db_exists(cursor, table_name)

    if exists and drop:
        cursor.execute("DROP TABLE "+table_name+";")
        print("Database \'"+table_name+"\' dropped")
        exists = None
   
    if not exists:
        ## Create Database
        cursor.execute( "CREATE TABLE "+table_name+"("+columns+");" )
        print("Database \'"+table_name+"\' created")
        if index:
            for idx in index:
                cursor.execute(idx)
                print("Index created for \'"+table_name+"\'")
        conn.commit()
    return

def populate_psql_database(conn, cursor, table_name, import_file):
    if ".gz" in import_file:
        tmp_file = "/tmp/"+table_name
        os.system("zcat "+import_file+" > "+tmp_file)
        cursor.execute("COPY "+table_name+" FROM \'"+tmp_file+"\' DELIMITER \'\t\' CSV HEADER;")
        os.unlink(tmp_file)
    else:
        cursor.execute("COPY "+table_name+" FROM \'"+import_file+"\' DELIMITER \'\t\' CSV HEADER;")

    print("File \'"+import_file+"\' imported to db \'"+table_name+"\'")

    # Alternative (but should skip first line [header])
    # os.system("zcat ../../int2mech/piv_app/static/data/common/DDI_db.tsv.gz | "+
    #           "psql -U bq_jgonzalez -d testdb -c 'COPY DDI_db FROM stdin'")

    conn.commit()
    return

def open_file(input_file, mode="r"):
    """ Open file Zipped or not
    """

    if re.search(".gz$", input_file):
        infile = gzip.open(input_file, mode)
    else:
        infile = open(input_file, mode)
    return infile

def print_table_max_lengths(table_file, sep="\t"):
    d = defaultdict(list)
    with open_file(table_file) as f:
        for i, line in enumerate(f):
            values = line.rstrip().split("\t")
            if i == 0:
                columns = values
            else:
                for col, val in zip(columns, values):
                    d[col].append(val)
    for col in columns:
        max_len = 0
        max_val = ""
        for val in d[col]:
            if len(val) > max_len:
                max_len = len(val)
                max_val = val
        print(col,"\t", max_len, "\t", max_val)

## Connect to PSQL (needs DB)
db = sys.argv[1]
conn = psycopg2.connect(database=db, 
                        user="bq_jgonzalez")
cursor = conn.cursor()

com_dir = "/net/home.isilon/ag-russell/bq_jgonzalez/int2mech/data/common/"
sp_dir = "/net/home.isilon/ag-russell/bq_jgonzalez/int2mech/data/species/"
if "all" in db:
    sps = ["Ath","Cel","Dme","Dre","Mmu","Sce","Xtr"]
elif "human" in db:
    sps = ["Hsa", "SARS2"]
drop = True

###########################################
## MAKE SURE ALL TSV FILES HAVE A HEADER ##
###########################################

# ### Pfam Data
# import_file = com_dir+"Pfam-A_r33-1.hmm.tsv.gz"
# table_name = "pfam_a_data"
# columns = '''
# accession   VARCHAR(10) PRIMARY KEY,
# type        TEXT        NOT NULL,
# identifier  TEXT        NOT NULL,
# description TEXT        NOT NULL
# '''
# create_psql_database(conn, cursor, table_name, columns, drop=drop)
# populate_psql_database(conn, cursor, table_name, import_file)


# ### DDI DATABASE
# import_file = com_dir+"DDI_db.tsv.gz"
# table_name = "domain_domain_ints"
# columns = '''
# pfam_acc_a  VARCHAR(10) NOT NULL,
# pfam_ide_a  VARCHAR(20) NOT NULL,
# pfam_acc_b  VARCHAR(10) NOT NULL,
# pfam_ide_b  VARCHAR(20) NOT NULL,
# source      VARCHAR(20) NOT NULL,
# pdbs        TEXT,
# CONSTRAINT  domain_domain_ints_pfam_acc_pair_pkey PRIMARY KEY (pfam_acc_a, pfam_acc_b)'''
# create_psql_database(conn, cursor, table_name, columns, drop=drop)
# populate_psql_database(conn, cursor, table_name, import_file)


# ### ELM CLASSES
# import_file = com_dir+"elm_Mar2020_classes.tsv"
# table_name = "elm_classes"
# columns = '''
# Accession           CHAR(10)        NOT NULL,
# Identifier          VARCHAR(30)     PRIMARY KEY,
# Name                VARCHAR(100)    NOT NULL,
# Description         TEXT            NOT NULL,
# Regex               TEXT            NOT NULL,
# Probability         REAL            NOT NULL,
# Instances           INT             NOT NULL,
# Instances_in_PDB    INT             NOT NULL'''
# create_psql_database(conn, cursor, table_name, columns, drop=drop)
# populate_psql_database(conn, cursor, table_name, import_file)


# ### ELM Interaction Domains
# import_file = com_dir+"elm_Mar2020_int_domains_final.tsv"
# table_name = "elm_domain_ints"
# columns = '''
# elm_ac                      TEXT    NOT NULL,
# elm_id                      TEXT    NOT NULL,
# domain_ids                  TEXT    NOT NULL,
# domain_name                 TEXT    NOT NULL,
# domain_description          TEXT    NOT NULL,
# present_in_taxon            TEXT,
# not_present_in_taxon        TEXT,
# elm_containing_genes        TEXT,
# dom_containing_genes_hsa    TEXT,
# dom_containing_genes        TEXT,
# phosphosites                TEXT,
# observations                TEXT,
# other_elm_required          TEXT
# '''
# create_psql_database(conn, cursor, table_name, columns, drop=drop)
# populate_psql_database(conn, cursor, table_name, import_file)


# ### 3did DMI
# import_file = com_dir+"3did_dmi_flat_edited_2020-01.tsv.gz"
# table_name = "domain_3did_motif_ints"
# columns = '''
# motif       TEXT    PRIMARY KEY NOT NULL,
# regex       TEXT                NOT NULL,
# domain_acc  VARCHAR(10)         NOT NULL,
# domain_name TEXT                NOT NULL,
# pdb_number  INT                 NOT NULL
# '''
# create_psql_database(conn, cursor, table_name, columns, drop=drop)
# populate_psql_database(conn, cursor, table_name, import_file)


# ### COSMIC mutations
# import_file = com_dir+"Cosmicv87_GenomeScreens_parsed.tsv.gz"
# table_name = "cosmic_genome_screens"
# columns = '''
# uniprot_acc     VARCHAR(15) NOT NULL,
# cosmic_id       VARCHAR(10) NOT NULL,
# enst            VARCHAR(20) NOT NULL,
# cds_mut         VARCHAR(20) NOT NULL,
# aa_mut          VARCHAR(10) NOT NULL,
# sample_num      INT         NOT NULL,
# cancer_types    TEXT        NOT NULL 
# '''
# indices = [
#     "CREATE INDEX cosmic_genomes_screens_uniprot_acc_idx ON "+table_name+" (uniprot_acc);"]
# create_psql_database(conn, cursor, table_name, columns, index=indices, drop=drop)
# populate_psql_database(conn, cursor, table_name, import_file)


# ######## ORGANISM-specific ##############

### ID MAPPING DATA
table_name = "protein_id_mapping"
columns = '''
protein_id   TEXT    NOT NULL,
uniprot_id   TEXT    NOT NULL,
uniprot_ids  TEXT    NOT NULL,
organism     TEXT    NOT NULL,
CONSTRAINT  protein_id_mapping_protein_id_organism_pkey PRIMARY KEY (protein_id, organism)
'''
indices = [
    "CREATE INDEX protein_id_mapping_organism_idx ON "+table_name+" (organism);",
    "CREATE INDEX protein_id_mapping_protein_id_idx ON "+table_name+" (protein_id);"]
create_psql_database(conn, cursor, table_name, columns, index=indices, drop=drop)
base_file = "id2uniprot_mapping_table_2020-05_SPS.tsv.gz"
for sp in sps:
    import_file = sp_dir+sp+"/"+base_file.replace("SPS", sp)
    populate_psql_database(conn, cursor, table_name, import_file)


### Protein Data
table_name = "protein_data"
columns = '''
uniprot_id  TEXT        PRIMARY KEY,
uniprot_acc TEXT        NOT NULL,
gene        TEXT        NOT NULL,
description TEXT        NOT NULL,
data_class  VARCHAR(15) NOT NULL,
organism    TEXT        NOT NULL,
length      INT         NOT NULL,
sequence    TEXT        NOT NULL,
biogrid_id  TEXT,
sorted_ints TEXT
'''
create_psql_database(conn, cursor, table_name, columns, drop=drop)
base_file = "protein_data_SPS.tsv.gz"
for sp in sps:
    import_file = sp_dir+sp+"/"+base_file.replace("SPS", sp)
    populate_psql_database(conn, cursor, table_name, import_file)


## UniProt Features
table_name = "uniprot_features"
columns = '''
uniprot_acc TEXT        NOT NULL,
type        VARCHAR(15) NOT NULL,
id          TEXT,
start_pos   INT         NOT NULL,
end_pos     INT         NOT NULL,
change      TEXT,
note        TEXT,
evidence    TEXT
'''
indices = ["CREATE INDEX uniprot_features_uniprot_acc_type_idx ON "+table_name+" (uniprot_acc, type);"]
create_psql_database(conn, cursor, table_name, columns, drop=drop, index=indices)
base_file = "uniprot_2020-05_features_SPS.tsv.gz"
for sp in sps:
    import_file = sp_dir+sp+"/"+base_file.replace("SPS", sp)
    populate_psql_database(conn, cursor, table_name, import_file)


### PFAM matches
table_name = "pfam_a_matches"
columns = '''
uniprot_acc TEXT        NOT NULL,
ali_start   INT         NOT NULL,
ali_end     INT         NOT NULL,
env_start   INT         NOT NULL,
env_end     INT         NOT NULL,
hmm_ac      VARCHAR(10) NOT NULL,
hmm_name    TEXT        NOT NULL,
type        VARCHAR(15) NOT NULL,
hmm_start   INT         NOT NULL,
hmm_end     INT         NOT NULL,
hmm_len     INT         NOT NULL,
bit_score   REAL        NOT NULL,
e_value     VARCHAR(15) NOT NULL,
clan        VARCHAR(15) NOT NULL,
CONSTRAINT pfam_a_matches_match_idx UNIQUE(uniprot_acc, env_start, env_end, hmm_ac)
'''
indices = [
    "CREATE INDEX pfam_a_matches_uniprot_acc_idx ON "+table_name+" (uniprot_acc);"]
create_psql_database(conn, cursor, table_name, columns, drop=drop, index=indices)
base_file = "Pfam-A_r33-1_matches_full_SPS.tsv.gz"
for sp in sps:
    import_file = sp_dir+sp+"/"+base_file.replace("SPS", sp)
    populate_psql_database(conn, cursor, table_name, import_file)


### LINEAR MOTIFS
table_name = "linear_motifs_matches"
columns = '''
uniprot_acc TEXT  NOT NULL,
motif       TEXT  NOT NULL,
source      TEXT  NOT NULL,
motif_start INT   NOT NULL,
motif_end   INT   NOT NULL,
sequence    TEXT  NOT NULL,
status      TEXT
'''
indices = ["CREATE INDEX linear_motifs_matches_uniprot_acc_idx ON "+table_name+" (uniprot_acc);",
           "CREATE INDEX linear_motifs_matches_uniprot_acc_source_idx ON "+table_name+" (uniprot_acc, source);"]
create_psql_database(conn, cursor, table_name, columns, index=indices, drop=drop)
base_file = "all_lm_hits_SQL_SPS.tsv.gz"
for sp in sps:
    import_file = sp_dir+sp+"/"+base_file.replace("SPS", sp)
    populate_psql_database(conn, cursor, table_name, import_file)


# ## PTMS
# table_name = "ptms"
# columns = '''
# uniprot_acc VARCHAR(15) NOT NULL,
# position    INT         NOT NULL,
# residue     VARCHAR(15) NOT NULL,
# ptm         VARCHAR(30) NOT NULL,
# sources     VARCHAR(30) NOT NULL
# '''
# indices = ["CREATE INDEX ptms_uniprot_acc_idx ON "+table_name+" (uniprot_acc);"]
# create_psql_database(conn, cursor, table_name, columns, index=indices, drop=drop)
# base_file = "ptms_SQL_SPS.tsv.gz"
# for sp in sps:
#     import_file = sp_dir+sp+"/"+base_file.replace("SPS", sp)
#     populate_psql_database(conn, cursor, table_name, import_file)


# ### PPI db
# table_name = "protein_protein_ints"
# columns = '''
# interaction_id  VARCHAR(25) PRIMARY KEY,
# protein_a       VARCHAR(25) NOT NULL,
# protein_b       VARCHAR(25) NOT NULL,
# throughput      VARCHAR(10) NOT NULL,
# pubmed          TEXT        NOT NULL
# '''
# indices = [
#     "CREATE INDEX protein_protein_ints_protein_a_protein_b_idx ON "+table_name+" (protein_a, protein_b);"]
# create_psql_database(conn, cursor, table_name, columns, index=indices, drop=drop)
# base_file = "ppi_db_SPS.tsv.gz"
# for sp in sps:
#     import_file = sp_dir+sp+"/"+base_file.replace("SPS", sp)
#     populate_psql_database(conn, cursor, table_name, import_file)


# ### Association Scores
# table_name = "association_scores"
# columns = '''
# type    TEXT        NOT NULL,
# acc_a   TEXT        NOT NULL,
# name_a  TEXT        NOT NULL,
# n_a     INT         NOT NULL,
# prob_a  REAL        NOT NULL,
# acc_b   TEXT        NOT NULL,
# name_b  TEXT        NOT NULL,
# n_b     INT         NOT NULL,
# prob_b  REAL        NOT NULL,    
# prob_ab REAL        NOT NULL,
# exp     REAL        NOT NULL,
# obs     INT         NOT NULL,
# p_val   REAL        NOT NULL,
# ratio   REAL        NOT NULL,
# lo      REAL        NOT NULL,
# organism VARCHAR(10) NOT NULL
# '''
# indices = [
#     "CREATE INDEX association_scores_acc_a_acc_b_organism_idx ON "+table_name+" (acc_a, acc_b, organism);"]
# create_psql_database(conn, cursor, table_name, columns, index=indices, drop=drop)
# base_file = "prot_ele_association_prob_edited_SPS.tsv.gz"
# for sp in sps:
#     import_file = sp_dir+sp+"/"+base_file.replace("SPS", sp)
#     populate_psql_database(conn, cursor, table_name, import_file)


### InterPreTS
# table_name = "interprets_pdb2019"
# columns = '''
# uniprot_id_a    TEXT        NOT NULL,
# uniprot_ac_a    VARCHAR(15) NOT NULL,
# pdb_a           VARCHAR(15),
# blast_eval_a    VARCHAR(15),
# blast_pcid_a    REAL,
# uni_start_a     INT,
# uni_end_a       INT,
# pdb_start_a     INT,
# pdb_end_a       INT,
# uniprot_id_b    TEXT        NOT NULL,
# uniprot_ac_b    VARCHAR(15) NOT NULL,
# pdb_b           VARCHAR(15),
# blast_eval_b    VARCHAR(15),
# blast_pcid_b    REAL,
# uni_start_b     INT,
# uni_end_b       INT,
# pdb_start_b     INT,
# pdb_end_b       INT,
# i2_raw          REAL,
# rand            VARCHAR(10),
# rand_mean       REAL,
# rand_sd         REAL,
# z_score         REAL,
# p_value         REAL,
# not_sure1       REAL,
# not_sure2       REAL
# '''
# indices = [
#     "CREATE INDEX interprets_pdb2019_uniprot_id_a_uniprot_id_b_idx ON "+table_name+" (uniprot_id_a, uniprot_id_b);"]
# create_psql_database(conn, cursor, table_name, columns, index=indices, drop=drop)
# base_file = "interprets_reviewed_results_SPS.txt.gz"
# for sp in sps:
#     import_file = sp_dir+sp+"/"+base_file.replace("SPS", sp)
#     populate_psql_database(conn, cursor, table_name, import_file)



conn.close()
sys.exit()






