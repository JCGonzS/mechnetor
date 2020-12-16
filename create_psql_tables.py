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
        print("Database", table_name, "dropped")
        exists = None
   
    if not exists:
        ## Create Database
        cursor.execute( "CREATE TABLE "+table_name+"("+columns+");" )
        print("Database",table_name,"created")
        if index:
            cursor.execute(index)
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

def get_elm_info(cursor, table_name):
    data = {}
    cursor.execute("SELECT Accession, Identifier, Name, Description, Regex, Probability "+
                   "FROM "+table_name+";")
    for row in cursor.fetchall():
        data[row[0]] = {
            "ide":   row[1],
            "name":  row[2],
            "des":   row[3],
            "regex": row[4],
            "prob":  row[5]
            }
    return data

def get_ddi(cursor, table_name):
    data = {}
    cursor.execute("SELECT pfam_acc_a, pfam_acc_b, source, pdbs "+
                    "FROM "+table_name+";") 
    for row in cursor.fetchall():
        pfam_a, pfam_b, source = row[0], row[1], row[2]
        pdbs   = [x for x in str(row[3]).split(";") if x!=""]
        # if pfam_a in all_pfams and pfam_b in all_pfams:
        data[(pfam_a, pfam_b)] = {"dbs":source, "pdbs": pdbs}
        data[(pfam_b, pfam_a)] = {"dbs":source, "pdbs": pdbs}
    return data

def get_ddi_2(cursor, table_name, a, b):
    a = "\'"+a+"\'"
    b = "\'"+b+"\'"
    cursor.execute("SELECT pfam_acc_a, pfam_acc_b, source, pdbs "+
                    "FROM "+table_name+
                    " WHERE pfam_acc_pair = ("+a+", "+b+")"
                    " OR pfam_acc_pair = ("+b+", "+a+");")
    row = cursor.fetchone()
    pdbs = [x for x in str(row[3]).split(";") if x!=""]
    data = {"dbs": row[2], "pdbs": pdbs}
    return data

def get_3did_dmi(cursor, table_name):
    data = {}
    cursor.execute("SELECT motif, domain_acc, regex, pdb_number "+
                    "FROM "+table_name+";")
    for row in cursor.fetchall():
        data[(row[0].upper(), row[1])] = (row[2], row[3])
    return data

def get_pfam_info(cursor, table_name):
    data = {}
    cursor.execute("SELECT accession, type, identifier, description "+
                   "FROM "+table_name+";")
    for row in cursor.fetchall():
        data[row[0]] = {
            "type": row[1],
            "ide":  row[2],
            "des":  row[3]
        }
    return data

def get_elm_dom(cursor, table_name):
    data = defaultdict(list)
    cursor.execute(
        "SELECT elm_id, domain_ids, present_in_taxon, not_present_in_taxon, elm_containing_genes, "+
        "dom_containing_genes_hsa, dom_containing_genes, phosphosites, observations, other_elm_required "+
        "FROM "+table_name+";")

    for row in cursor.fetchall():
        elm = row[0].upper()
        doms = string2list_fix(row[1].upper())
        for dom in doms:
            data[(elm, dom)].append( {
                "in_taxon"        : string2list_fix(row[2]),
                "not_in_taxon"    : string2list_fix(row[3]),
                "elm_genes"       : string2list_fix(row[4]),
                "human_dom_genes" : string2list_fix(row[5]),
                "dom_genes"       : string2list_fix(row[6]),
                "phos"            : string2list_fix(row[7]),
                "obs"             : format_none(row[8]),
                "req_elms"        : string2list_fix(row[9])
            })

    return data

def format_none(string):
    if string==None:
        return ""
    else:
        return string

def string2list_fix(string):
    return [str(x) for x in str(string).upper().replace(" ","").split(",") if x!="" and x!="NONE"]

def get_PPI(cursor, table_name, acc_a, acc_b):
    #TRY WITH CONSTRAINT KEY TO COMPARE SPEED
    acc_a = "\'"+acc_a+"\'"
    acc_b = "\'"+acc_b+"\'"
    cursor.execute("SELECT interaction_id FROM "+table_name+
                    " WHERE (interactor_a = "+acc_a+" AND interactor_b = "+acc_b+")"+
                    " OR (interactor_a = "+acc_b+" AND interactor_b = "+acc_a+");")
    ints = []
    for row in cursor.fetchall():
        ints.append(row[0])

    return ints


## Connect to PSQL (needs DB)
conn = psycopg2.connect(database="piv", 
                        user="bq_jgonzalez")
cursor = conn.cursor()

com_dir = "/net/home.isilon/ag-russell/bq_jgonzalez/int2mech/data/common/"
sp_dir = "/net/home.isilon/ag-russell/bq_jgonzalez/int2mech/data/species/"
sps = ["Sce"]#,"Hsa"]

###########################################
## MAKE SURE ALL TSV FILES HAVE A HEADER ##
###########################################

# prots_a = ["118093", "113010", "112055", "107985"]
# prots_b = ["117030", "110358", "111384", "107028"]
# for _ in range(0, 10):
#     for a in prots_a:
#         for b in prots_b:
#             ints = get_PPI(cursor, "ppi_db", a, b)

## Pfam Data
import_file = com_dir+"Pfam-A_r33-1.hmm.tsv.gz"
table_name = "pfam_a_data"
columns = '''
accession   VARCHAR(10) PRIMARY KEY,
type        TEXT        NOT NULL,
identifier  TEXT        NOT NULL,
description TEXT        NOT NULL
'''
create_psql_database(conn, cursor, table_name, columns)
#populate_psql_database(conn, cursor, table_name, import_file)

## ELM CLASSES
import_file = com_dir+"elm_Mar2020_classes.tsv"
table_name = "elm_classes"
columns = '''
Accession           CHAR(10)                NOT NULL,
Identifier          VARCHAR(30) PRIMARY KEY NOT NULL,
Name                VARCHAR(100)            NOT NULL,
Description         TEXT                    NOT NULL,
Regex               TEXT                    NOT NULL,
Probability         REAL                    NOT NULL,
Instances           INT                     NOT NULL,
Instances_in_PDB    INT                     NOT NULL'''
create_psql_database(conn, cursor, table_name, columns)
# populate_psql_database(conn, cursor, table_name, import_file)

## ELM Interaction Domains
import_file = com_dir+"elm_Mar2020_int_domains_final.tsv"
table_name = "elm_int_dom"
columns = '''
elm_ac                      TEXT    NOT NULL,
elm_id                      TEXT    NOT NULL,
domain_ids                  TEXT    NOT NULL,
domain_name                 TEXT    NOT NULL,
domain_description          TEXT    NOT NULL,
present_in_taxon            TEXT,
not_present_in_taxon        TEXT,
elm_containing_genes        TEXT,
dom_containing_genes_hsa    TEXT,
dom_containing_genes        TEXT,
phosphosites                TEXT,
observations                TEXT,
other_elm_required          TEXT
'''
create_psql_database(conn, cursor, table_name, columns)
# populate_psql_database(conn, cursor, table_name, import_file)

# ## DDI DATABASE
# import_file = com_dir+"DDI_db.tsv.gz"
# table_name = "ddi_db"
# columns = '''
# pfam_acc_a  VARCHAR(10) NOT NULL,
# pfam_ide_a  VARCHAR(20) NOT NULL,
# pfam_acc_b  VARCHAR(10) NOT NULL,
# pfam_ide_b  VARCHAR(20) NOT NULL,
# source      VARCHAR(20) NOT NULL,
# pdbs        TEXT,
# CONSTRAINT  pfam_acc_pair PRIMARY KEY (pfam_acc_a,pfam_acc_b)'''
# create_psql_database(conn, cursor, table_name, columns)
# # populate_psql_database(conn, cursor, table_name, import_file)

# ## 3did DMI
# import_file = com_dir+"3did_dmi_flat_edited_2020-01.tsv.gz"
# table_name = "dmi_3did"
# columns = '''
# motif       TEXT    PRIMARY KEY NOT NULL,
# regex       TEXT                NOT NULL,
# domain_acc  VARCHAR(10)         NOT NULL,
# domain_name TEXT                NOT NULL,
# pdb_number  INT                 NOT NULL
# '''
# create_psql_database(conn, cursor, table_name, columns)
# # populate_psql_database(conn, cursor, table_name, import_file)

# ## COSMIC mutations
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
# create_psql_database(conn, cursor, table_name, columns)
# # populate_psql_database(conn, cursor, table_name, import_file)


# ######## ORGANISM-specific ##############

# ## ID MAPPING DATA
# table_name = "id_mapping"
# columns = '''
# id          TEXT    NOT NULL,
# uniprot_id  TEXT    NOT NULL,
# uniprot_ids TEXT    NOT NULL,
# organism    TEXT    NOT NULL,
# CONSTRAINT  id_and_org PRIMARY KEY (id, organism)
# '''
# create_psql_database(conn, cursor, table_name, columns)

# base_file = "id2uniprot_mapping_table_2020-05_SPS.tsv.gz"
# for sp in sps:
#     import_file = sp_dir+sp+"/"+base_file.replace("SPS", sp)
#     # populate_psql_database(conn, cursor, table_name, import_file)

# ## Protein Data
# table_name = "protein_data"
# columns = '''
# uniprot_id  TEXT        PRIMARY KEY,
# uniprot_acc VARCHAR(15) NOT NULL,
# gene        TEXT        NOT NULL,
# description TEXT        NOT NULL,
# data_class  VARCHAR(15) NOT NULL,
# organism    TEXT        NOT NULL,
# length      INT         NOT NULL,
# sequence    TEXT        NOT NULL,
# biogrid_id  TEXT,
# sorted_ints TEXT
# '''
# create_psql_database(conn, cursor, table_name, columns)
# base_file = "protein_data_SPS.tsv.gz"
# for sp in sps:
#     import_file = sp_dir+sp+"/"+base_file.replace("SPS", sp)
#     # populate_psql_database(conn, cursor, table_name, import_file)

# table_name = "uniprot_features"
# columns = '''
# uniprot_acc VARCHAR(15) NOT NULL,
# type        VARCHAR(15) NOT NULL,
# id          TEXT,
# start_pos       INT         NOT NULL,
# end_pos         INT         NOT NULL,
# note        TEXT,
# evidence    TEXT
# '''
# index = "CREATE INDEX idx_main ON uniprot_features (uniprot_acc, type);"
# create_psql_database(conn, cursor, table_name, columns, drop=True, index=index)
# base_file = "uniprot_2020-05_features_SPS.tsv.gz"
# for sp in sps:
#     import_file = sp_dir+sp+"/"+base_file.replace("SPS", sp)
#     populate_psql_database(conn, cursor, table_name, import_file)


# ## PFAM matches
# table_name = "pfam_a_matches"
# columns = '''
# uniprot_acc VARCHAR(15) NOT NULL,
# ali_start   INT         NOT NULL,
# ali_end     INT         NOT NULL,
# env_start   INT         NOT NULL,
# env_end     INT         NOT NULL,
# hmm_ac      VARCHAR(10) NOT NULL,
# hmm_name    TEXT        NOT NULL,
# type        VARCHAR(15) NOT NULL,
# hmm_start   INT         NOT NULL,
# hmm_end     INT         NOT NULL,
# hmm_len     INT         NOT NULL,
# bit_score   REAL        NOT NULL,
# e_value     VARCHAR(15) NOT NULL,
# clan        VARCHAR(15) NOT NULL,
# CONSTRAINT Pfam_Match UNIQUE(uniprot_acc, env_start, env_end, hmm_ac)
# '''
# create_psql_database(conn, cursor, table_name, columns, drop=True)
# base_file = "Pfam-A_r33-1_matches_full_SPS.tsv.gz"
# for sp in sps:
#     import_file = sp_dir+sp+"/"+base_file.replace("SPS", sp)
#     populate_psql_database(conn, cursor, table_name, import_file)

# ### PTMS
# # uniprot_acc", "residue", "ptm", "sources
# table_name = "ptms"
# columns = '''
# uniprot_acc VARCHAR(15) NOT NULL,
# position    INT         NOT NULL,
# residue     VARCHAR(15) NOT NULL,
# ptm         VARCHAR(30) NOT NULL,
# sources     VARCHAR(30) NOT NULL
# '''
# create_psql_database(conn, cursor, table_name, columns)
# base_file = "ptms_SQL_SPS.tsv.gz"
# for sp in sps:
#     import_file = sp_dir+sp+"/"+base_file.replace("SPS", sp)
#     # populate_psql_database(conn, cursor, table_name, import_file)


# ### LINEAR MOTIFS
# table_name = "lm_hits"
# columns = '''
# uniprot_acc TEXT  NOT NULL,
# motif       TEXT  NOT NULL,
# source      TEXT  NOT NULL,
# motif_start INT   NOT NULL,
# motif_end   INT   NOT NULL,
# sequence    TEXT  NOT NULL,
# status      TEXT
# '''
# create_psql_database(conn, cursor, table_name, columns)

# base_file = "all_lm_hits_SQL_SPS.tsv.gz"
# for sp in sps:
#     import_file = sp_dir+sp+"/"+base_file.replace("SPS", sp)
#     # populate_psql_database(conn, cursor, table_name, import_file)


# ## PPI db (all SPs)
# table_name = "ppi_db"
# columns = '''
# interaction_id  VARCHAR(15) PRIMARY KEY,
# interactor_a    VARCHAR(15) NOT NULL,
# interactor_b    VARCHAR(15) NOT NULL,
# throughput      VARCHAR(10) NOT NULL,
# pubmed          TEXT        NOT NULL
# '''
# create_psql_database(conn, cursor, table_name, columns, drop=True)
# # cursor.execute("CREATE INDEX idx_main ON ppi_db (interactor_a, interactor_b);")
# base_file = "ppi_db_SPS.tsv.gz"
# for sp in sps:
#     import_file = sp_dir+sp+"/"+base_file.replace("SPS", sp)
#     populate_psql_database(conn, cursor, table_name, import_file)


### InterPreTS
table_name = "interprets_pdb2019"
columns = '''
uniprot_id_a    TEXT        NOT NULL,
uniprot_ac_a    VARCHAR(15) NOT NULL,
pdb_a           VARCHAR(15),
blast_eval_a    VARCHAR(15),
blast_pcid_a    REAL,
uni_start_a     INT,
uni_end_a       INT,
pdb_start_a     INT,
pdb_end_a       INT,
uniprot_id_b    TEXT        NOT NULL,
uniprot_ac_b    VARCHAR(15) NOT NULL,
pdb_b           VARCHAR(15),
blast_eval_b    VARCHAR(15),
blast_pcid_b    REAL,
uni_start_b     INT,
uni_end_b       INT,
pdb_start_b     INT,
pdb_end_b       INT,
i2_raw          REAL,
rand            VARCHAR(10),
rand_mean       REAL,
rand_sd         REAL,
z_score         REAL,
p_value         REAL,
not_sure1       REAL,
not_sure2       REAL
'''
index = "CREATE INDEX interprets_2019_main_idx ON "+table_name+" (uniprot_ac_a, uniprot_ac_b);"
drop = True
create_psql_database(conn, cursor, table_name, columns, index=index, drop=drop)
base_file = "interprets_reviewed_results_SPS.txt.gz"
for sp in sps:
    import_file = sp_dir+sp+"/"+base_file.replace("SPS", sp)
    populate_psql_database(conn, cursor, table_name, import_file)

conn.close()
sys.exit()




# ###
# pfam_info = get_pfam_info(cursor, "pfam_a_data")
# elm_info = get_elm_info(cursor, "elm_classes")
# elm_dom = get_elm_dom(cursor, "elm_int_dom")
# ddi = get_ddi(cursor, "ddi_db")
# dmi = get_3did_dmi(cursor, "dmi_3did")
# ###








