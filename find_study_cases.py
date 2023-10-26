import sys, os, re, gzip, psycopg2
from collections import defaultdict

def open_file(input_file, mode="rt"):
    """ Open file Zipped or not
    """
    if re.search(".gz$", input_file):
        infile = gzip.open(input_file, mode)
    else:
        infile = open(input_file, mode)
    return infile

def identify_protein(cursor, id_map, protein_id, org):
    if org:
        cursor.execute("SELECT uniprot_id FROM "+id_map+
                       " WHERE (id = \'"+protein_id+"\' AND organism = \'"+org+"\');")
        results = cursor.fetchone()
        if results:
            return results[0], org
        else:
            return None, None
    else:
        cursor.execute("SELECT uniprot_id, organism FROM "+id_map+
                       " WHERE (id = \'"+protein_id+"\');")
        results = cursor.fetchall()
        if results and len(results) == 1:
            return results[0], results[1]
        else:
            return None, None

def get_protein_data_sql(cursor, sql_table, uni_id):
    cursor.execute("SELECT uniprot_acc, gene, description, length, biogrid_id, sorted_ints"+
                    " FROM "+sql_table+
                    " WHERE (uniprot_id = \'"+uni_id+"\');")
    row = cursor.fetchone()
    bioids = ""
    if row[4]:
        bioids = row[4]
    ints = []
    if row[5]:
        row[5]
        ints = row[5].split(", ")
    data = {
        "uni_ac":       row[0],
        "genes":        row[1].split("; "),
        "description":  row[2],
        "length":       row[3],
        "biogrid_ids":  bioids,
        "sorted_ints":  ints
    }
    return data

def get_pfam_matches_sql(cursor, sql_table, uni_ac):
    matches = []
    cursor.execute("SELECT env_start, env_end, hmm_ac, e_value"+
                    " FROM "+sql_table+
                    " WHERE (uniprot_acc = \'"+uni_ac+"\');")
    for row in cursor.fetchall():
        matches.append({
            "start": int(row[0]),
            "end": int(row[1]),
            "acc": row[2],
            "e-value": float(row[3])
        })
    return matches

def get_lm_hits_sql(cursor, sql_table, uni_ac):
    data = defaultdict(list)
    cursor.execute("SELECT motif, source, motif_start, motif_end, sequence, status"+
                " FROM "+sql_table+
                " WHERE (uniprot_acc = \'"+uni_ac+"\');")
    for row in cursor.fetchall():
        data[row[0]].append({
            "source": row[1],
            "start":  row[2],
            "end":    row[3],
            "seq":    row[4],
            "status": row[5]
        })
    return data

def get_uni_vars_sql(cursor, sql_table, uni_ac):
    data = defaultdict(list)
    uni_dis = []
    for role in ["VARIANT"]:
        cursor.execute("SELECT type, start_pos, end_pos, note"+
                    " FROM "+sql_table+
                    " WHERE (uniprot_acc=\'"+uni_ac+"\' AND type=\'"+role+"\');");
        for row in cursor.fetchall():
            try:
                var, info = re.search("(.+)\((.+)\)", row[3]).group(1, 2)
            except:
                continue
    
            if len([x for x in ["cancer","melanoma","carcinoma"] if x in info])==0:
                dis = re.findall("in ([A-Z][^\s]+)", info)
                if len(dis)>0:
                    var = var.replace(" ","")
                    wt, mt = var.split("->")
                    start = int(row[1])
                    end = int(row[2])
                    data[row[0]].append({
                        "var": var,
                        "info": info,
                        "start": start,
                        "end": end,
                        "dis": dis
                    })
                    uni_dis += dis
                    print(uni_ac+"/"+wt+str(end)+mt)

                        # print(uni_ac, start,end,var,dis)
    return data, set(uni_dis)


conn = psycopg2.connect(database="piv", 
                            user="bq_jgonzalez")
cursor = conn.cursor()

cosmic_file = "/net/home.isilon/ag-russell/bq_jgonzalez/Projects/background_mut_models/MutantExport_Parsed.tsv.gz"

target = sys.argv[1]
a, b = get_uni_vars_sql(cursor, "uniprot_features", target)
sys.exit()

census = []
with open_file("Census_all2020.tsv.txt") as f:
    for line in f:
        t = line.rstrip().split("\t")
        gene = t[0]
        uni_id, org = identify_protein(cursor, "id_mapping", gene, "HUMAN")
        if uni_id:
            census.append(uni_id)

protein_data, pfam_matches, lms = {}, {}, {}
uni_var = {}
uni_acs = []
for uni_id in census:
    protein_data[uni_id] = get_protein_data_sql(cursor, "protein_data", uni_id)
    uni_ac = protein_data[uni_id]["uni_ac"]
    uni_acs.append(uni_ac)
    # pfam_matches[uni_id] = get_pfam_matches_sql(cursor, "pfam_a_matches", uni_ac)        
    # lms[uni_id] = get_lm_hits_sql(cursor, "lm_hits", uni_ac)

    uni_var[uni_id], uni_dis = get_uni_vars_sql(cursor, "uniprot_features", uni_ac)
    # if len(uni_dis)>1:
    #     print(uni_id, uni_dis)


# with gzip.open(cosmic_file):


