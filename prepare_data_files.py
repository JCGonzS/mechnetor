#!/usr/bin/python3
import sys, re, os, gzip, math, pprint, datetime, random, string, itertools
import site; site.getsitepackages()
import os.path
from collections import defaultdict
from Bio import SwissProt, SeqIO
import find_all_slims
# import pandas as pd
# import run_interprets

def open_file(input_file, mode="rt"):
    """ Open file Zipped or not
    """
    if re.search(".gz$", input_file):
        infile = gzip.open(input_file, mode)
    else:
        infile = open(input_file, mode)
    return infile

def check_file_exists(somefile):
    if os.path.isfile(somefile):
        print_status(somefile, "exists")
        return True
    else:
        sys.exit("FileError: '"+somefile+"' does not exist.\n"+
            "Make sure the file is located and named correctly\n")

def print_status(thefile, status):
    thefile = thefile.split("/")[-1]
    l = 80 - len(thefile) - len(status)
    print ("'"+thefile+"'"+"."*l+status)

def edit_pfam_dat(pfam_dat_file, cov2_pfam_file, out_file):
    """ Writes a simple TSV file with a few selected columns from Pfam data.
    """
    done = []
    with open_file(out_file, "wt") as out:
        cols = ["Accession", "Type", "Identifier", "Description"]
        out.write("\t".join(cols)+"\n")

        for pfam_file in [cov2_pfam_file, pfam_dat_file]: 
            with open_file(pfam_file) as f:
                for line in f:
                    if line.startswith("#=GF ID"):
                        ide = re.search("#=GF\s+ID\s+(.+)", line.rstrip()).group(1)
                    elif line.startswith("#=GF AC"):
                        acc = re.search("#=GF\s+AC\s+(.+)", line.rstrip()).group(1)
                        acc = acc.split(".")[0]
                    elif line.startswith("#=GF DE"):
                        des = re.search("#=GF\s+DE\s+(.+)", line.rstrip()).group(1)
                    elif line.startswith("#=GF TP"):
                        tp = re.search("#=GF\s+TP\s+(.+)", line.rstrip()).group(1)
                    elif line.startswith("//"):
                        if acc not in done:
                            out.write("\t".join([acc, tp, ide, des])+"\n")
                            done.append(acc)

def read_pfam_dat(pfam_dat_file):
    pfam_names, pfam_accs, pfam_types = {}, {}, {}
    with open_file(pfam_dat_file) as f:
        for line in f:
            t = line.rstrip().split("\t")
            pfam_names[t[0]] = t[2]
            pfam_accs[t[2]] = t[0]
            pfam_types[t[0]] = t[1]
    return pfam_names, pfam_accs, pfam_types

def extract_protein_data_from_uniprot_text(uniprot_file):
    """Extracts information from a UniProt text file.
    Requires Bio::SwissProt module
    """

    data = {}
    primary_id_map = defaultdict(list)
    secondary_id_map = defaultdict(list)
    uni_ac_seqs, masks = {}, {}
    modres = defaultdict(lambda: defaultdict(dict))
    uni_features = defaultdict(lambda: defaultdict(list))
    region = defaultdict(list)
    pdb2uni = defaultdict(set)
    reviewed = set()
    previous_data_class = ""
    s=set()
    for record in SwissProt.parse(open_file(uniprot_file)):

        if previous_data_class == "Unreviewed" and record.data_class == "Reviewed":
            print("WARNING: UniProt entries not sorted Reviewed->Unreviewed")

        uni_id = record.entry_name.upper()
        uni_acs = record.accessions

        main_genes, genes, syns = [], [], []
        main_gene = "N/A"
        if record.gene_name.strip():
           
            ## Main gene names
            main_genes = [match.split()[0] for match in re.findall("Name=([^;]+)",
                                                            record.gene_name)]
        
            ## ORF & Locus names
            orfnames = []
            for match in (re.findall("ORFNames=([^;]+)", record.gene_name)+
                    re.findall("OrderedLocusNames=([^;]+)", record.gene_name)):
                for name in match.split(", "):
                    name = re.search("([^{]+)", name).group(1).replace(" ","")
                    if not name.startswith("ECO:"):
                        for name2 in name.split("/"):
                            orfnames.append(name2)

            ## Other Synonyms
            for match in re.findall("Synonyms=([^;]+)", record.gene_name):
                for syn in match.split(", "):
                    syn = re.search("([^{]+)", syn).group(1).replace(" ","")
                    if not syn.startswith("ECO:"):
                        syns.append(syn)

            genes = main_genes+orfnames+syns
            if len(main_genes)>0:
                main_gene = "; ".join(main_genes)                                                    
            else:
                main_gene = genes[0]

        main_uni_ac = uni_acs[0].upper()
        data[uni_id] = {
            "gn": main_gene,
            "ac": main_uni_ac,
            "des": re.search("Name: Full=([^;|{]+)", record.description).group(1),
            "dc": record.data_class,
            "seq": record.sequence,
            "biogrid_id": set(),
            "sorted_ints": []
        }
        uni_ac_seqs[main_uni_ac] = record.sequence

        for key in [uni_id, main_uni_ac]+main_genes:
            if uni_id not in primary_id_map[key]:
                primary_id_map[key].append(uni_id)

        for key in uni_acs[1:]+genes[1:]:
            if uni_id not in secondary_id_map[key]:
                secondary_id_map[key].append(uni_id)
        
        masks[uni_id] = ["0"] * len(record.sequence)
        
        # for ref in record.cross_references:
        #     if ref[0].upper()=="BIOGRID":
        #         data[uni_id]["biogrid_id"].append(ref[1])
        #         data[uni_id]["biogrid_ints"] += int(ref[2])

    # Feature types:    
    # 'REGION', 'CROSSLNK', 'STRAND', 'BINDING', 'TOPO_DOM', 'MOTIF', 'NP_BIND', 'SIGNAL', 'VARIANT', 
    # 'PROPEP', 'INIT_MET', 'TURN', 'DOMAIN', 'TRANSIT', 'CONFLICT', 'NON_CONS', 'CHAIN', 'HELIX', 
    # 'MOD_RES', 'CARBOHYD', 'DNA_BIND', 'VAR_SEQ', 'REPEAT', 'LIPID', 'NON_TER', 'ACT_SITE', 'MUTAGEN', 
    # 'CA_BIND', 'COILED', 'SITE', 'COMPBIAS', 'INTRAMEM', 'PEPTIDE', 'NON_STD', 'ZN_FING', 'METAL', 'DISULFID',
    #  'TRANSMEM'
        for feat in record.features:
            start = str(feat.location.start)
            end = str(feat.location.end)
            ac = str(feat.location).split("[")[0]
            s.add(feat.type)

            if feat.type=="MOD_RES" and not ac and end!="UnkownPosition()":
                end = int(end)
                info = feat.qualifiers["note"].split(";")[0]
                if "Phospho" in info:
                    modres[main_uni_ac][end]["Phosphorylation"] = ["UniProt"]
                elif "acetyl" in info:
                    modres[main_uni_ac][end]["Acetylation"] = ["UniProt"]
                # else:
                #     modres[main_uni_ac]["other"][str(end)+"-"+res] = info

            elif feat.type in ["VARIANT", "MUTAGEN", "METAL", "BINDING"]:
                evidence = feat.qualifiers.get("evidence", "")
                note = feat.qualifiers.get("note", "")
                ide = feat.id
                if not feat.id:
                    ide = "-"
                uni_features[main_uni_ac][feat.type].append( (start, end, ide, evidence, note) )

        #     elif feat[0]=="REGION":
        #         region[uni_ac].append(
        #             {"start": start,
        #              "end": end,
        #              "info": info
        #             }
        #         )
    for key in secondary_id_map:
        for uni_id in secondary_id_map[key]:
            if uni_id not in primary_id_map[key]:
                primary_id_map[key].append(uni_id) 

    return data, primary_id_map, uni_ac_seqs, modres, uni_features #, region

def create_uni_features_table(outfile, uni_feats):
    with open_file(outfile, "wt") as out:
        vals = ["#Uni_Ac", "Type", "Id", "Start", "End", "Note", "Evidence"]
        out.write("\t".join(vals)+"\n")
        for uni_ac in uni_feats:
            for feat_type in uni_feats[uni_ac]:
                for feat in uni_feats[uni_ac][feat_type]:
                    (start, end, feat_id, evidence, note) = feat
                    vals = [uni_ac, feat_type, feat_id, start, end, note, evidence]
                    out.write("\t".join(vals)+"\n")
    return

def add_uniprot_biogrid_mapping(infile, uni_map, prot_data):
    with open_file(infile) as f:
        for line in f:
            if line[0]!="#":
                t = line.rstrip().split("\t")
                uni_ac, bio_id, ints = t[0], t[1], int(t[2])
                if uni_ac in uni_map:
                    uni_id = uni_map[uni_ac][0]
                    prot_data[uni_id]["biogrid_id"].add(bio_id)
    return prot_data

def get_bio_2_uni_mapping(prot_data):
    bio2uni = defaultdict(set)
    for uni_id in prot_data:
        bio_ids = sorted(list(prot_data[uni_id]["biogrid_id"]))
        if len(bio_ids) > 0:
            for bio_id in bio_ids:
                bio2uni[bio_id].add(uni_id)
    return bio2uni

def extract_biogrid_interactions(biogrid_file, bio2uni):
    biogrid_ppi_all = defaultdict(lambda: defaultdict(set))
    biogrid_ppi = defaultdict(lambda: defaultdict(set))
    interaction_info = {}

    with open_file(biogrid_file) as f:
        for line in f:
            t = line.rstrip().split("\t")
            if line[0] == "#":
                cols = t
            else:
                # for i,(c, x) in enumerate(zip(cols,t)):
                #     print(i, c, x)
            
                interaction_id = t[0]
                bio_id_a, bio_id_b = t[3], t[4]
                pubmed, throughput = t[14], t[17].replace(" Throughput","")
                # org_a, org_b = t[-2], t[1]
                # uni_acs_a, uni_acs_b = t[23].split("|"), t[26].split("|")
    
                if bio_id_a in bio2uni and bio_id_b in bio2uni:

                    biogrid_ppi_all[bio_id_a][bio_id_b].add(interaction_id)
                    biogrid_ppi_all[bio_id_b][bio_id_a].add(interaction_id)

                    if bio_id_a < bio_id_b:
                        biogrid_ppi[bio_id_a][bio_id_b].add(interaction_id)
                    else:
                        biogrid_ppi[bio_id_b][bio_id_a].add(interaction_id)
                        
                    interaction_info[interaction_id] = {
                        "pubmed": pubmed,
                        "exp": throughput
                    }

    return biogrid_ppi, biogrid_ppi_all, interaction_info

def create_ppi_database(outfile, bio_ppi, interaction_info):
    with open_file(outfile, "wt") as out:
        out.write("\t".join(["Int_ID", "Bio_ID_A", "Bio_ID_B", "Throughput", "PubMed"])+"\n")

        for bio_a in bio_ppi:
            for bio_b in bio_ppi[bio_a]:
                # they are sorted bio_a < bio_b
                for int_id in bio_ppi[bio_a][bio_b]:
                    pm = interaction_info[int_id]["pubmed"]
                    th = interaction_info[int_id]["exp"]
                    cols = [int_id, bio_a, bio_b, th, pm]
                    out.write("\t".join(cols)+"\n")
    return

def add_sorted_interactors(prot_data, biogrid_ppi_all, bio2uni):
    for uni_id in prot_data:
        bio_ids = sorted(list(prot_data[uni_id]["biogrid_id"]))
        if len(bio_ids) > 0:
            pp = defaultdict(set)
            for bio_id_a in bio_ids:
                for bio_id_b in biogrid_ppi_all[bio_id_a]:
                    pp[bio_id_b].update(biogrid_ppi_all[bio_id_a][bio_id_b])
            
            for bio_id in sorted(pp, key=lambda k: len(pp[k]), reverse=True):
                bio_uni_id = sorted(list(bio2uni[bio_id]))[0]
                # ints = len(pp[bio_id])
                prot_data[uni_id]["sorted_ints"].append(bio_uni_id)
    return prot_data

def create_protein_data_table(outfile_name, protein_data, org):

    with open_file(outfile_name, "wt") as out:
        cols = ["Uni_Id", "Uni_Ac", "Gene", "Description", "Data_Class", 
                "Organism", "Seq_Length", "Sequence", "Biogrid IDs"]
        out.write("\t".join(cols)+"\n")

        for uni_id in protein_data:
            data = protein_data[uni_id]
            vals = [uni_id, data["ac"], data["gn"], data["des"], data["dc"], 
                    org, str(len(data["seq"])), data["seq"], 
                    ", ".join(data["biogrid_id"]), ", ".join(data["sorted_ints"])]
            out.write("\t".join(vals)+"\n")
    return

def create_id_map_table(outfile, id_map, org):

    with open_file(outfile, "wt") as out:
        cols = ["Id", "Uni_Id", "Uni_Ids", "Organism"]
        out.write("\t".join(cols)+"\n")

        for key in id_map:
            ids = id_map[key]
            vals = [key, ids[0], "; ".join(ids), org]
            out.write("\t".join(vals)+"\n")

def create_ptms_table(outfile, ptms, uni_seqs):
    with open_file(outfile, "wt") as out:
        out.write("\t".join(["uniprot_acc", "position", "residue", "ptm", "sources"])+"\n")
        for uni_ac in ptms:
            for pos in sorted(ptms[uni_ac]):
                res = uni_seqs[uni_ac][pos-1]
                for ptm_type in ptms[uni_ac][pos]:
                    sources = ",".join(ptms[uni_ac][pos][ptm_type])
                    out.write("\t".join([uni_ac, str(pos), res, ptm_type, sources])+"\n")
    return

def extract_pfam_doms(pfam_file, max_eval=999):
    """Pfam-A matches in species proteome. File downloaded from PFAM.
    """
    pfams = defaultdict(lambda: defaultdict(set) )
    pfams_temp = defaultdict(lambda: defaultdict(set) )

    with open_file(pfam_file) as f:
        for line in f:
            if line[0] != "#":
                t = line.rstrip().split("\t")

                uni_ac = t[0].upper()
                start, end = int(t[3]), int(t[4])
                pfam_ac, domain_e_val = t[5].split(".")[0], float(t[12])

                pfams[uni_ac][pfam_ac].add((start, end, domain_e_val))

    return pfams

def run_pfamscan(PFAMSCAN_FILE, uni_ac_seqs, pfam_matches):
    pfamscan_script  = "pfamscan.py"
    tmp_dir          = "/tmp/" # This is the systems TEMP directory
    tmp_fasta        = tmp_dir+"no_pfam_seqs"
    tmp_pfamout      = tmp_dir+"pfamout"
    evalue           = "0.001"
    email            = "juan-carlos.gonzalez@bioquant.uni-heidelberg.de"

    n, i = 0, 0
    for uni_ac, seq in uni_ac_seqs.items():
        if uni_ac not in pfam_matches:
            i += 1
            if i > 100:
                n += 1
                i = 1
            fasta_file = tmp_fasta+"_"+str(n).zfill(3)+".fasta"
            with open_file(fasta_file, "at") as out:
                out.write(">"+uni_ac+"\n")
                out.write(seq+"\n")

    for x in range(n+1):
        x = str(x).zfill(3)
        print("Running PfamScan #"+x)
        fasta = tmp_fasta+"_"+x+".fasta"
        pfamout = tmp_pfamout+"_"+x
        if not os.path.isfile(tmp_pfamout+"_"+x+".out.txt"):
            cmd = ("python3 {} --sequence {} --database pfam-a --evalue {}".format(
                    pfamscan_script, fasta, evalue)+
                    " --format txt --outfile {} --email {} --quiet".format(
                    pfamout, email))
            os.system(cmd)
            os.unlink(fasta)

    miss = []
    with open_file(PFAMSCAN_FILE, "at") as f:
        for x in range(n+1):
            x = str(x).zfill(3)
            pfamout_file = tmp_pfamout+"_"+x+".out.txt"
            if os.path.isfile(pfamout_file):
                with open_file(pfamout_file) as f2:
                    for line in f2:
                        if line[0]!="#" and line.strip():
                            f.write(line)
                os.unlink(pfamout_file)
                os.unlink(tmp_pfamout+"_"+x+".sequence.txt")
            else:
                miss.append(x)

    if len(miss)>0:
        print("Missing Pfamscan for files", miss)

def merge_pfam_files(full_file, matches_file, pfamscan_file):
    with open_file(full_file, "wt") as out:
        cols = ["#uni_ac", "ali_start", "ali_end", "env_start", "env_end", "hmm_ac",
                "hmm_name", "type", "hmm_start", "hmm_end", "hmm_len", "bit_score",
                "e-value", "clan"]
        out.write("\t".join(cols)+"\n")

        with open_file(matches_file, "rt") as f:
            for line in f:
                if line.strip() and line[0]!="#":
                    out.write(line)

        with open_file(pfamscan_file, "rt") as f:
            for line in f:
                if line.strip() and line[0]!="#":
                    t = line.rstrip().split()
                    out.write( "\t".join(t[:5]+[t[5].split(".")[0]]+t[6:-2]+t[-1:])+"\n")
    return

def create_ddi_database(pfam_names, pfam_int_file, db3did_file, out_file):

    dom_ints = defaultdict(set)
    with open_file(pfam_int_file, "rt") as f:
        for line in f:
            a, b = line.rstrip().split()
            pair = (a, b)
            if b < a:
                pair = (b, a)
            dom_ints[pair].add("Pfam")

    pdbs_3did = {}
    with open_file(db3did_file, "rt") as f:
        for line in f:
            #=ID    1-cysPrx_C      1-cysPrx_C       (PF10417.4@Pfam       PF10417.4@Pfam)
            if line.startswith("#=ID"):
                pfam_acc_a = line.rstrip().split()[3].split(".")[0].split("(")[1]
                pfam_acc_b = line.rstrip().split()[4].split(".")[0]
                pdbs = set()
   
            #=3D    1n8j    E:153-185       O:153-185       0.99    1.35657 0:0
            elif line.startswith("#=3D"):
                pdb = line.rstrip().split()[1]
                pdbs.add(pdb)

            elif line.startswith("//"):
                pair = (pfam_acc_a, pfam_acc_b)
                if pfam_acc_b < pfam_acc_a:
                    pair = (pfam_acc_b, pfam_acc_a)
                dom_ints[pair].add("3did")
                pdbs_3did[pair] = ";".join(sorted(list(pdbs)))

    with open_file(out_file, "wt") as out:
        cols = ["Pfam_Acc_A", "Pfam_Ide_A", "Pfam_Acc_B", "Pfam_Ide_B", "Source", "PDBs"]
        out.write("\t".join(cols)+"\n")
        for pair in dom_ints:
            a, b = pair
            if a in pfam_names and b in pfam_names:
                name_a, name_b = pfam_names[a], pfam_names[b]
                sources = ", ".join(sorted(list(dom_ints[pair])))
                pdbs = pdbs_3did.get(pair, "")
                out.write("\t".join([a, name_a, b, name_b, sources, pdbs])+"\n")
    return

def get_dmi_3did(dmi_file, pfam_acc):
    outlines = []
    lm_3did = {}
    all_pdbs = set()
    with open_file(dmi_file) as f:
        for line in f:
            t = line.rstrip().split("\t")
            if line.startswith("#=ID"):
                pfam, lm = t[1], t[3]
                pdb_set = set()
            elif line.startswith("#=PT"):
                regex = t[1].split()[0]
            elif line.startswith("#=3D"):
                pdb = t[1].upper()+":"+t[3].split(":")[0]
                seq = t[4]
                pdb_set.add(pdb+":"+t[3].split(":")[1]+":"+t[4])
                all_pdbs.add(pdb)
            elif line.startswith("//"):
                if pfam == "Pkinase_Tyr":
                    pfam = "PK_Tyr_Ser-Thr"
                if pfam == "CTD_bind":
                    pfam = "CID"
                outlines.append( "\t".join([lm, regex, pfam_acc[pfam], pfam, str(len(pdb_set)) ]) )
                lm_3did[lm] = {"ide": lm, "regex": regex, "prob": "-", "pdbs":list(pdb_set)}

    return outlines, lm_3did, all_pdbs

def parse_elm_classes(elm_classes_file):
    elm_classes = {}
    elm_map = {}
    with open_file(elm_classes_file) as f:
        for line in f:
            if line.startswith("ELM"):
                t = line.rstrip().split("\t")
                elm_map[t[0]] = t[1].upper()
                elm_map[t[1].upper()] = t[0]
                elm_classes[t[0]] = {"ide": t[1], "regex": t[4], "prob": t[5]}
    return elm_map, elm_classes

def edit_elm_interactions(old_file, new_file, elm_map):
    elm_acs = []
    with open_file(new_file, "wt") as out:
        tsv = pd.read_csv(old_file, sep='\t')
        for elm_id in tsv["ELM identifier"]:
            elm_acs.append(elm_map.get(elm_id.upper(), "N/A"))
        tsv.insert(loc=0, column="ELM accesion", value=elm_acs)
        tsv.to_csv(new_file, sep="\t", index=False)

def extract_elm_interactions(elm_intdom_file):
    elm_int = defaultdict(set)
    with open_file(elm_intdom_file, "rt") as f:
        for i, line in enumerate(f):
            t = line.rstrip().split("\t")
            elm_ac, doms = t[0], t[2]
            if elm_ac != "N/A":
                for dom in doms.split(", "):
                    elm_int[elm_ac].add(dom)
                    elm_int[dom].add(elm_ac)
    return elm_int

def extract_elm_instances(instances_file, elm_map, target_acs):
    tp = defaultdict(lambda: defaultdict(list))
    fp = defaultdict(lambda: defaultdict(list))
    with open_file(instances_file, "rt") as f:
        for line in f:
            if line[0]!="#" and not line.startswith("\"Acc"):
                t = line.rstrip().replace("\"", "").split("\t")
                elm_id = t[2].upper()
                uni_acs = t[5].replace("-1", "").split()
                start, end = int(t[6]), int(t[7])
                status = t[10]
                elm_ac = elm_map[elm_id]
                for uni_ac in uni_acs:
                    if uni_ac in target_acs:
                        if status == "true positive":
                            tp[uni_ac][elm_ac].append((start, end))
                            fp[uni_ac][elm_ac].append((start, end))
    return tp, fp

def extract_slim_hits(hits_file, tp, fp):
    elms = defaultdict(lambda: defaultdict(list) )
    with open_file(hits_file, "rt") as f:
        for line in f:
            if line[0] != "#" and line.strip():
                t = line.rstrip().split("\t")
                lm_name, lm_ac = t[0], t[1]
                uni_ac, uni_id = t[2].split("|")[1:]
                start, end = int(t[4]), int(t[5])
                prob_score = t[6] # elm probability
                seq = t[-1]
                status = ""
                if (start, end) in tp.get(uni_ac, {}).get(lm_ac, []):
                    status = "TP"
                elif (start, end) in fp.get(uni_ac, {}).get(lm_ac, []):
                    status = "FP"
                elms[uni_ac][lm_ac].append({
                    "start": start, "end": end, 
                    "seq": seq,
                    "status": status
                    })
    return elms

def get_pdb2uniprot_map(mapfile, allowed_prts):
    d = defaultdict(lambda: defaultdict(list))
    with open_file(mapfile, "rt") as f:
        for line in f:
            if line[0]=="#" or line.startswith("PDB"):
                continue
            t = line.rstrip().split()
            pdb = t[0].upper()+":"+t[1]
            seq_start, seq_end = t[3], t[4]
            pdb_start, pdb_end = t[5], t[6]
            uni = t[2]
            uni_start, uni_end = int(t[7]), int(t[8])
            if uni in allowed_prts:
                d[pdb][uni].append( ";".join(t[3:9]) )
    return d

def get_3did_dmi_instances(lms, pdb2uni, uni_seqs):
    tp = defaultdict(lambda: defaultdict(list))
    for lm in lms:
        pdbs = lms[lm]["pdbs"]
        for pdb in pdbs:
            t = pdb.split(":")
            pdb_id, pdb_coor, seq = t[0]+":"+t[1], t[2], t[3]

            p = re.compile(seq)
            if pdb_id in pdb2uni:
                for uni_ac in pdb2uni[pdb_id]:
                    for coord in pdb2uni[pdb_id][uni_ac]:
                        uni_start, uni_end = coord.split(";")[4:]
                        sub_seq = uni_seqs[uni_ac][int(uni_start)-1:int(uni_end)]
              
                        for m in p.finditer(sub_seq):
                            m_start = int(uni_start)+int(m.start())
                            m_end = m_start + len(m.group()) - 1
                            tp[uni_ac][lm].append((m_start, m_end))
    return tp

def reformat_slim_hit_file(edited_file, elms, lm3d):
    with open_file(edited_file, "wt") as out:
        cols = ["#Uni_Ac", "LM_Ac", "Source", "Start", "End", "Seq", "Status"]
        out.write("\t".join(cols)+"\n")

        for lms, source in zip([elms, lm3d], ["ELM", "3did"]):
            for uni_ac in lms:
                for lm_ac in lms[uni_ac]:
                    for hit in lms[uni_ac][lm_ac]:
                        cols = [uni_ac, lm_ac, source, str(hit["start"]), str(hit["end"]), 
                                hit["seq"], hit["status"]]
                        out.write("\t".join(cols)+"\n")
    return

def edit_interprets_raw_file(raw_file, rev_file, uni_map, uni_seqs, prot_data):
    hits = defaultdict(lambda: defaultdict(list))
    s = set()
    with open_file(raw_file, "rt") as f:
        for i, line in enumerate(f):
            t = line.rstrip().split("\t")
            if i==0:
                cols = t
            elif len(t) == 24:
                ### They are sorted: uniac_a < uniac_b
                ac_a, ac_b = t[0], t[8]
                id_a, id_b = uni_map[ac_a][0], uni_map[ac_b][0]
                ac_a, ac_b = prot_data[id_a]["ac"], prot_data[id_b]["ac"]
                start_a, end_a = int(t[4]), int(t[5])
                start_b, end_b = int(t[12]), int(t[13])
                len_a = len(uni_seqs[ac_a])
                len_b = len(uni_seqs[ac_b])
                if start_a > len_a:
                    continue
                elif end_a > len_a:
                    end_a = len_a
                if start_b > len_b:
                    continue
                elif end_b > len_b:
                    end_b = len_b
                info_a = [ac_a]+t[1:5]+[str(end_a)]+t[6:8] 
                info_b = [ac_b]+t[9:13]+[str(end_b)]+t[14:16] 
                scores = t[16:]
                z = float(t[-4])
                if id_a <= id_b:
                    hits[(id_a, id_b)][z].append({
                        "info_a": info_a, "info_b": info_b, "scores": scores                
                        })
                else:
                    hits[(id_b, id_a)][z].append({
                        "info_a": info_b, "info_b": info_a, "scores": scores                
                        })

    mask, mask2 = {}, {}
    final = defaultdict(list)
    for pair in hits:
        (id_a, id_b) = pair
        for z in sorted(hits[pair], reverse=True):
            for hit in hits[pair][z]:
                # Add if first one, compare with previous ones if it's not
                if pair not in final:
                    final[pair].append(hit)
                else:
                    ac_a, ac_b = prot_data[id_a]["ac"], prot_data[id_b]["ac"]
                    start_a, end_a = int(hit["info_a"][4]), int(hit["info_a"][5])
                    start_b, end_b = int(hit["info_b"][4]), int(hit["info_b"][5])
                    mask[ac_a] = ["0"]*len(uni_seqs[ac_a])
                    mask[ac_b] = ["0"]*len(uni_seqs[ac_b])
                    mask2[ac_a] = fill_mask(start_a, end_a, ["0"]*len(uni_seqs[ac_a]))
                    mask2[ac_b] = fill_mask(start_b, end_b, ["0"]*len(uni_seqs[ac_b]))
                    
                    for fhit in final[pair]:
                        fstart_a, fend_a = int(fhit["info_a"][4]), int(fhit["info_a"][5])
                        fstart_b, fend_b = int(fhit["info_b"][4]), int(fhit["info_b"][5])
                        mask[ac_a] = fill_mask(fstart_a, fend_a, mask[ac_a])
                        mask[ac_b] = fill_mask(fstart_b, fend_b,mask[ac_b])
                        ov2_a = calculate_overlap(fstart_a, fend_a, mask2[ac_a])
                        ov2_b = calculate_overlap(fstart_b, fend_b, mask2[ac_b])
                        if ov2_a>=0.75 and ov2_b>=0.75:
                            break
                    else:      
                        ov_a = calculate_overlap(start_a, end_a, mask[ac_a])
                        ov_b = calculate_overlap(start_b, end_b, mask[ac_b])
                        if ov_a<0.75 and ov_b<0.75:
                            final[pair].append(hit)

    with open_file(rev_file, "wt") as out:
        v = ["Uni-ID1"]+cols[1:8]+["Uni-ID2"]+cols[9:]
        for (id_a, id_b) in final:
            for hit in final[(id_a, id_b)]:
                v = [id_a]+hit["info_a"]+[id_b]+hit["info_b"]+hit["scores"]
                out.write("\t".join(v)+"\n")

def fill_mask(start, end, mask):
    for i in range(start-1, end):
        mask[i] = "1"
    return mask

def calculate_overlap(start, end, mask):
    sub_mask = mask[start-1:end]
    overlap = sub_mask.count("1") / float(len(sub_mask))
    return overlap
#################################################

def merge_gene_names(genes):
    """For a list of similar gene names ("GENE1A", "GENE2B", "GENE3C"),
    returns a string following the format: "GENE(1A,2B,3C)"
    """
    pattern = ""
    flag = 0
    for i in range(1, len(genes[0])+1):
        reg = genes[0][:i]
        for g in genes[1:]:
            if reg != g[:i]:
                flag = 1
        if flag == 0:
            pattern = reg
        else:
            break
    var = []
    for g in genes:
        v = g.split(pattern)[1]
        if v.strip():
            var.append(v)

    return pattern+"("+",".join(var)+")"





def compare_ints(t, gene1, gene2, ints):
    s1, e1 = int(t[4]), int(t[5])
    s2, e2 = int(t[12]), int(t[13])
    for info in ints[gene1][gene2]:
        i = info.split("\t")
        sA, eA = int(i[4]), int(i[5])
        sB, eB = int(i[12]), int(i[13])

        mask1 = [0 for _ in range(0,eA)]
        mask2 = [0 for _ in range(0,eB)]
        if e1 > eA:
            mask1 = [0 for _ in range(0,e1)]
        if e2 > eB:
            mask2 = [0 for _ in range(0,e2)]

        mask1 = mask1[:sA-1]+[n+1 for n in mask1[sA-1:eA]]+mask1[eA:]
        mask1 = mask1[:s1-1]+[n+1 for n in mask1[s1-1:e1]]+mask1[e1:]
        overlap1 = mask1.count(2)

        mask2 = mask2[:sB-1]+[n+1 for n in mask2[sB-1:eB]]+mask2[eB:]
        mask2 = mask2[:sB-1]+[n+1 for n in mask2[sB-1:eB]]+mask2[eB:]
        overlap2 = mask2.count(2)


        print("OLD",gene1,":",sA,"-",eA,"\t",gene2,":",sB,"-",eB)
        print("NEW",gene1,":",s1,"-",e1,"\t",gene2,":",s2,"-",e2,"\n")
 
def create_protein_data_json(prot_dict, alt_ids, pfams, elms, elms2, ptms,
                             uni_features, regions, bio_id, bio_set,
                             sp, outfile_name, mode="mongo"):

    pp = pprint.PrettyPrinter(indent=4)
    json_data = {}
    with open_file(outfile_name, "w") as out:
        for uni_ac in prot_dict["seq"]:
            seq = prot_dict["seq"][uni_ac]
            protein_data = {
                    "uni_ac": uni_ac,
                    "uni_id": list(prot_dict["ID"][uni_ac])[0], # always, only 1
                    "genes": list(prot_dict["GN"][uni_ac]),
                    "alt_ids": alt_ids[uni_ac],
                    "description": prot_dict["des"][uni_ac],
                    "data_class": prot_dict["dc"][uni_ac],
                    "length": len(seq),
                    "sequence": seq,
                    "organism": sp,
                    "pfams": [],
                    "elms": [],
                    "3dlms": [],
                    "phosphorylation": [],
                    "acetylation": [],
                    "uni_features": [],
                    "regions": [],
                    "biogrid_id": "NA",
                    "biogrid_interactors": []
            }

            if uni_ac in bio_id:
                protein_data["biogrid_id"] = bio_id[uni_ac]
                sorted_ints = sorted(bio_set[uni_ac], key=bio_set[uni_ac].get,
                                     reverse=True)
                protein_data["biogrid_interactors"] = sorted_ints

            for pfam_ac in sorted(pfams.get(uni_ac, [])):
                for (start, end, e_val) in sorted(pfams[uni_ac][pfam_ac],
                                                  key=lambda x: int(x[0])):
                    protein_data["pfams"].append(
                        {
                            "acc":   pfam_ac,
                            "start": int(start),
                            "end":   int(end),
                            "e-val": float(e_val)
                        })

            for elm_ac in sorted(elms.get(uni_ac, [])):
                for (start, end, status) in sorted(elms[uni_ac][elm_ac],
                                                key=lambda x: int(x[0])):
                    protein_data["elms"].append(
                        {
                            "acc":    elm_ac,
                            "start":  int(start),
                            "end":    int(end),
                            "seq":    seq[int(start)-1:int(end)],
                            "status": status
                        })
            
            for elm_ac in sorted(elms2.get(uni_ac, [])):
                for (start, end, status) in sorted(elms2[uni_ac][elm_ac],
                                                key=lambda x: int(x[0])):
                    protein_data["3dlms"].append(
                        {
                            "acc":    elm_ac,
                            "start":  int(start),
                            "end":    int(end),
                            "seq":    seq[int(start)-1:int(end)],
                            "status": status
                        })

            for pos_res in ptms.get(uni_ac, []):
                pos, res = pos_res.split("-")

                for mod, source in ptms[uni_ac][pos_res].iteritems():
                    if mod == "p":
                        ptm_type = "phosphorylation"
                    elif mod == "ac":
                        ptm_type = "acetylation"
                    protein_data[ptm_type].append(
                        {
                            "pos": int(pos),
                            "res": res,
                            "source": ";".join(source)
                        })
            # else:
            #     for i in range(9):
            #         iso_uni_ac = uni_ac+"-"+str(i)
            #         if iso_uni_ac in ptms:
            #             for (res, pos, mod) in sorted(ptms[iso_uni_ac],
            #                                          key=lambda x: int(x[1])):
            #                 if mod == "p":
            #                     ptm_type = "phosphorylation"
            #                 elif mod == "ac":
            #                     ptm_type = "acetylation"
            #                 protein_data[ptm_type].append(
            #                     {
            #                         "pos" : int(pos),
            #                         "res" : res,
            #                         "source": "PSP,-"+str(i)
            #                     })
            #             break

            for feat in uni_features:
                for pos in uni_features[feat].get(uni_ac, []):
                    protein_data["uni_features"].append(
                        {
                            "feat": feat,
                            "pos": pos,
                            "info": ";".join(uni_features[feat][uni_ac][pos])
                        }
                    )

            protein_data["regions"] = regions.get(uni_ac, [])

            if mode == "mongo":
                out.write(str(protein_data)+"\n")
            elif mode == "json":
                json_data[uni_ac] = protein_data

        if mode == "json":
            json.dump(json_data, out)

def print_uniprot_species_map(prot_dict, alt_ids, sp, outfile):
    with open_file(outfile, "w") as f:
        for uni_ac in prot_dict["seq"]:
            uni_id = list(prot_dict["ID"][uni_ac])[0]
            alts = alt_ids[uni_ac]
            for x in [uni_ac]+[uni_id]+alts:
                print(x+"\t"+sp)

def hasNumbers(inputString):
	return any(char.isdigit() for char in inputString)

def assign_Number(inputString):
	if ";" in inputString:
		return inputString.split(";")[0]
	else:
		return inputString

def extract_ptms(psp_file, prot_dict, modres):
    psp = defaultdict(set)
    psp_temp = defaultdict(set)
    with open_file(psp_file) as f:
        for line in f:
            t = line.rstrip().split("\t")
            pos, mod = t[4].split("-")
            res, pos = re.search("(\w)([\d]+)", pos).group(1,2)
            lt_lit = "."
            ms_lit = "."
            ms_cst = "."
            cst_cat = "."
            uni_ac = t[2]
            if len(t)>10 and hasNumbers(t[10]):
                lt_lit = assign_Number(t[10])
            if len(t)>11 and hasNumbers(t[11]):
                ms_lit = assign_Number(t[11])
            if len(t)>12 and hasNumbers(t[12]):
                ms_cst = assign_Number(t[12])
            if len(t)>13 and hasNumbers(t[13]):
                cst_cat = assign_Number(t[13])

            if uni_ac in prot_dict["seq"]:
                seq = prot_dict["seq"][uni_ac]
                if len(seq)>=int(pos) and seq[int(pos)-1]==res:
                    psp[uni_ac].add( (res, pos, mod) )
                    if (uni_ac in modres and pos+"-"+res in modres[uni_ac]
                    and mod in modres[uni_ac][pos+"-"+res]):
                        modres[uni_ac][pos+"-"+res][mod].append("PSP")
                    else:
                        modres[uni_ac][pos+"-"+res][mod] = ["PSP"]
            elif uni_ac in prot_dict["AC"]:
                psp_temp[uni_ac].add( (res, pos, mod) )

    for alt_ac in psp_temp:
        for uni_ac in prot_dict["AC"][alt_ac]:
            if uni_ac not in psp:
                seq = prot_dict["seq"][uni_ac]
                for (res,pos,mod) in psp_temp[alt_ac]:
                    if len(seq)>=int(pos) and seq[int(pos)-1]==res:
                        psp[uni_ac].add( (res, pos, mod) )
                        if (uni_ac in modres and pos+"-"+res in modres[uni_ac]
                        and mod in modres[uni_ac][pos+"-"+res]):
                            modres[uni_ac][pos+"-"+res][mod].append("PSP")
                        else:
                            modres[uni_ac][pos+"-"+res][mod] = ["PSP"]

    return modres

def get_HIPPIE_UniProt_map(map_file, prot_dict):
    """
    # Source of PPI: HIPPIE
    # not all ids from the hippie db matched a uniprot_ac, thus Uniprot mapping
    # tool was used and the results need to be read before:
    """
    hippie_map = defaultdict(set)
    with open_file(map_file) as f:
        for line in f:
            if (line.startswith("yourlist") or "Deleted." in line
                or "Homo sapiens" not in line):
                continue
            t = line.rstrip().split("\t")
            uni_ac = t[1]
            uni_ac = list(prot_dict["AC"][uni_ac])[0] #always 1 (checked)
            for protein in t[0].split(","):
                hippie_map[protein].add(uni_ac)

    return hippie_map

def extract_HIPPIE_interactions(hippie_file, hippie_map, score_cutoff=0.63):
    """
        (from HIPPIE) Medium confidence = 0.63 / High confidence = 0.73
    """
    ppi = defaultdict(set)
    with open_file(hippie_file) as f:
        for line in f:
            t = line.rstrip().split("\t")
            protsA = t[0].split(",")
            protsB = t[2].split(",")
            score = float(t[4])
            if score >= score_cutoff:
                for protA in protsA:
                    if protA in hippie_map:
                        a = list(hippie_map[protA])[0]
                        for protB in protsB:
                            if protB in hippie_map:
                                b = list(hippie_map[protB])[0]
                                if b < a:
                                    a = list(hippie_map[protB])[0]
                                    b = list(hippie_map[protA])[0]
                                ppi[a].add(b)
    return ppi

def extract_ddi_interactions(ddi_file):
    ints = defaultdict(dict)
    with open_file(ddi_file) as f:
        for i, line in enumerate(f):
            if i > 0:
                t = line.rstrip().split("\t")
                ac_a, ac_b, source = t[0], t[2], t[4]
                ints[ac_a][ac_b] = source
                ints[ac_b][ac_a] = source

    return ints

def count_protein_and_domain_pairs(ppi, eles):
    # From non-redundant PPI set, get:
    #  - total number of proteins
    #  - total number of interactions (or protein pairs)

    ## WARNING: THIS TAKES VERY LONG
    all_proteins = set()
    total_interactions = 0
    ele_pair_count = defaultdict(lambda: defaultdict(int))
    for a in ppi:
        all_proteins.add(a)
        eles_A = list(eles[a])
        for b in ppi[a]:
            all_proteins.add(b)
            total_interactions += 1
            eles_B = list(eles[b])
            ele_pairs = set()
            for (ele_a, ele_b) in itertools.product(eles_A, eles_B):
                eleA, eleB = ele_a, ele_b
                if ele_b <= ele_a:
                    eleA, eleB = ele_b, ele_a
    
                if (eleA, eleB) not in ele_pairs:
                    ele_pair_count[eleA][eleB] += 1
                    ele_pairs.add((eleA, eleB))

    return all_proteins, total_interactions, ele_pair_count

def count_individual_element_frequency(all_proteins, eles):
    element_count = defaultdict(int)
    for protein in all_proteins:
        for ele in eles[protein]:
            element_count[ele] += 1

    return element_count

def create_prob_file(out_file, total_proteins, total_interactions,
                     ele_pair_count, ele_count, ele_names,
                     ints_ddi, ints_elm, ints_dmi, min_npair=1):
    """ Calculate element-element* probabilities and log-odds
        *elements are Pfam domains and ELMs

        Only for those appearing in some PPI. To the rest we can assign the
        lowestLO score)

        Parameters:
            - Minimum number of protein pairs with the domain-domain signature
    """

    lo = set()
    # with gzip.open(out_file, "wb") as out:
    with open_file(out_file, "w") as out:
        line = ("{:<8} {:<10} {:<24} {:<6} {:<10}"+
                " {:<10} {:<24} {:<6} {:<10}"+
                " {:<9} {:<12} {:<8} {:<9} {:<12} {:<8}").format(
                "#Type", "Acc_A", "Name_A", "n_A", "Prob_A",
                          "Acc_B", "Name_B", "n_B", "Prob_B",
                "Prob_A&B", "Exp", "Obs", "P-val", "Ratio", "LO")
        # print line
        cols = ["#Type", "Acc_A", "Name_A", "n_A", "Prob_A",
                         "Acc_B", "Name_B", "n_B", "Prob_B",
                "Prob_A&B", "Exp", "Obs", "P-val", "Ratio", "LO"]
        out.write("\t".join(cols)+"\n")


        for ele_a in sorted(ele_pair_count):
            for ele_b in sorted(ele_pair_count[ele_a]):
                if ele_b in ints_ddi.get(ele_a, []):
                    cat = ints_ddi[ele_a][ele_b]
                elif ele_b in ints_elm.get(ele_a, []):
                    cat = "ELM"
                elif ele_b in ints_dmi.get(ele_a, []):
                    cat = "3did:DMI"
                else:
                    cat = "Pred"

                prob_a = float(ele_count[ele_a]) / total_proteins
                prob_b = float(ele_count[ele_b]) / total_proteins
                prob_ab = prob_a * prob_b
                exp = prob_ab * total_interactions
                obs = ele_pair_count[ele_a][ele_b]
                p_val = float(obs)/total_interactions
                ratio = float(obs) / exp
                log2 = math.log(ratio, 2)

                # if obs >= min_npair:
                line = ("{:<8} {:<10} {:<24} {:<6} {:<.3E} "+
                     " {:<10} {:<24} {:<6} {:<.3E}"+
                     " {:<.3E} {:<12.3f} {:<8} {:<.3E} {:<12.3f} {:<8.4f}").format(
                    cat, ele_a, ele_names[ele_a], ele_count[ele_a], prob_a,
                    ele_b, ele_names[ele_b], ele_count[ele_b], prob_b,
                    prob_ab, exp, obs, p_val, ratio, log2)
                # print line
                cols = [cat,
                        ele_a, ele_names[ele_a], str(ele_count[ele_a]), str(prob_a),
                        ele_b, ele_names[ele_b], str(ele_count[ele_b]), str(prob_b),
                        str(prob_ab), str(exp), str(obs), str(p_val),
                        str(ratio), str(log2)]
                out.write("\t".join(cols)+"\n")

def create_interprets_file(outfile, ppi, prot_dict, iprets_dir):

    with open_file(outfile, "w") as out:
        cols = ["#Gene1","PDB1","Blast-E1","Blast-PCID1","qstart1","qend1",
                "pdbstart1","pdbend1","Gene2","PDB2","Blast-E2","Blast-PCID2",
                "qstart2","qend2","pdbstart2","pdbend2","i2-raw","rand",
                "rand-mean","rand-sd","Z","p-value","not-sure1","not-sure2"]
        out.write("\t".join(cols)+"\n")

    for a in sorted(ppi):
        input_seqs = {}
        input_seqs[a] = prot_dict["seq"][a]
        pairs = []
        for b in sorted(ppi[a]):
            pairs.append((a,b))
            input_seqs[b] = prot_dict["seq"][b]

        ide = ''.join(random.choice(string.ascii_uppercase +
                        string.ascii_lowercase + string.digits) for _ in range(8))
        print( "Running InterPrets for {} and {} partners: {}".format(a,
                                        str(len(ppi[a])), ";".join(ppi[a])) )
        try:
            run_interprets.main(input_seqs, iprets_dir, outfile, ide,
                                these_pairs=pairs, fasta_as_input=False,
                                verbose=False)
        except:
            print("InterPreTS failed")
        print( "[IntePreTs finished for {} protein pairs]".format(len(pairs)) )

def merge_protein_elements(*argv):
    eles = defaultdict(set)
    for arg in argv:
        for prot in arg:
            for k in arg[prot]:
                eles[prot].add(k)

    return eles

def extract_3did_dmi_interactions(dmi_file):
    dmi = defaultdict(set)
    with open_file(dmi_file) as f:
        for i, line in enumerate(f):
            t = line.rstrip().split("\t")
            if i > 0:
                motif, dom = t[0], t[2]
                dmi[motif].add(dom)
                dmi[dom].add(motif)
    return dmi

def merge_protein_elements_names(pfam_names, elm_names, lm_3did):
    names = pfam_names
    names.update(elm_names)
    for lm in lm_3did.keys():
        names[lm] = lm
    return names

def main( SP="Hsa",
          DATA_DIR="",
          max_overlap=0.2,
          no_overlap_between_doms_and_lms=True,
          force_prot_data=False,
          force_ppi=False,
          force_probs=False,
          force_iprets=False
        ):

    organisms = {
    "Ath": ["ARATH", "3702", "Arabidopsis thaliana"],
    "Cel": ["CAEEL", "6239"],
    "Dme": ["DROME", "7227"],
    "Dre": ["DANRE", "7955"],
    "Hsa": ["HUMAN", "9606", "Homo sapiens"],
    "HUMAN": ["HUMAN", "9606", "Homo sapiens"],
    "Mmu": ["MOUSE", "10090"],
    "Sce": ["YEAST", "559292", "Saccharomyces cerevisiae"],
    "Xla": ["XENLA", "8355"],
    "Xtr": ["XENTR", "8364"]
    }
    org = organisms[SP][0]

    ### DEFINE required data files
    COM_DIR         = DATA_DIR+"common/"
    SP_DIR          = DATA_DIR+"species/"+SP+"/"
    UNI_VERSION     = "2020-05" # Release
    PFAM_VERSION    = "r33-1"   # Release 33.1
    ELM_VERSION     = "Mar2020" #"Oct2019"
    DB3DID_VERSION  = "2020-01"
    BIOGRID_VERSION = "4.2.191" # release
    #
    PSP_VERSION     = "Mar2019"

    ## Common files:
    ### PFAM
    PFAM_DAT_FILE        = COM_DIR+"Pfam-A_"+PFAM_VERSION+".hmm.dat.gz"
    EDITED_PFAM_DAT_FILE = COM_DIR+"Pfam-A_"+PFAM_VERSION+".hmm.tsv.gz"
    PFAM_INT_FILE        = COM_DIR+"Pfam-A_"+PFAM_VERSION+".interactions.txt.gz"
    COV2_PFAM_SCAN_FILE  = COM_DIR+"SARSCoV2_pfam_matches.scan.txt"
    COV2_PFAM_FILE       = COM_DIR+"Pfam-A.SARS-CoV-2.full.gz"
    ### ELM
    ELM_CLASSES_FILE     = COM_DIR+"elm_"+ELM_VERSION+"_classes.tsv"
    ELM_INSTANCES_FILE   = COM_DIR+"elm_"+ELM_VERSION+"_instances.tsv"
    ELM_DOM_FILE         = COM_DIR+"elm_"+ELM_VERSION+"_int_domains_reviewed.tsv"
    ELM_DOM_EDIT_FILE    = COM_DIR+"elm_"+ELM_VERSION+"_int_domains_final.tsv"

    ### 3did
    FLAT_3DID_FILE       = COM_DIR+"3did_flat_"+DB3DID_VERSION+".gz"
    DMI_3DID_FILE        = COM_DIR+"3did_dmi_flat_"+DB3DID_VERSION+".gz"
    EDITED_DMI_3DID_FILE = COM_DIR+"3did_dmi_flat_edited_"+DB3DID_VERSION+".tsv.gz"
    ### DDI database = Pfam + 3did DDI's
    DDI_FILE             = COM_DIR+"DDI_db.tsv.gz"
    PDB2UNIPROT_FILE     = COM_DIR+"pdb_chain_uniprot_2020-11.tsv.gz"
    BIOGRID_FILE         = COM_DIR+"BIOGRID-ALL-"+BIOGRID_VERSION+".tab3.gz"
    UNI2BIO_FILE         = COM_DIR+"UniProt_to_BioGRID_"+BIOGRID_VERSION+".tab.txt"

    ## Species files:
    UNI_TEXT_FILE     = SP_DIR+"uniprot_"+UNI_VERSION+"_proteome_"+SP+".txt.gz"
    UNI_FASTA_FILE    = SP_DIR+"uniprot_"+UNI_VERSION+"_proteome_"+SP+".fasta.gz"
    PROT_DATA_FILE    = SP_DIR+"protein_data_"+SP+".tsv.gz"
    PFAM_MATCHES_FILE = SP_DIR+"Pfam-A_"+PFAM_VERSION+"_matches_"+SP+".tsv.gz"
    PFAMSCAN_FILE     = SP_DIR+"PfamScan_"+PFAM_VERSION+"_"+SP+".txt.gz"
    PFAM_FULL_FILE    = SP_DIR+"Pfam-A_"+PFAM_VERSION+"_matches_full_"+SP+".tsv.gz"
    ID_MAP_FILE       = SP_DIR+"id2uniprot_mapping_table_"+UNI_VERSION+"_"+SP+".tsv.gz"
    ELM_HITS_FILE     = SP_DIR+"elm_hits_"+ELM_VERSION+"_"+SP+".tsv.gz"
    LM_3DID_HITS_FILE = SP_DIR+"3did_slims_"+DB3DID_VERSION+"_"+SP+".tsv.gz"
    LM_HITS_SQL_FILE  = SP_DIR+"all_lm_hits_SQL_"+SP+".tsv.gz"
    PTM_SQL_FILE      = SP_DIR+"ptms_SQL_"+SP+".tsv.gz"
    PPI_FILE          = SP_DIR+"ppi_db_"+SP+".tsv.gz"
    UNI_FEAT_FILE     = SP_DIR+"uniprot_"+UNI_VERSION+"_features_"+SP+".tsv.gz"
    IPRETS_RAW_FILE   = SP_DIR+"interprets_raw_results_"+SP+".txt.gz"
    IPRETS_REV_FILE   = SP_DIR+"interprets_reviewed_results_"+SP+".txt.gz"
    ####---------updated til here


    ASSOC_PROB_FILE   = SP_DIR+"prot_ele_association_prob_"+SP+".tsv.gz"    
    HIPPIE_FILE       = SP_DIR+"hippie_v2-2.tsv.gz" # Hsa
    HIPPIE_MAP_FILE   = SP_DIR+"hippie_v2-2_uniprot_mapping_table.tsv.gz" # Hsa
    PSP_FILE          = SP_DIR+"PSP_ptms_"+PSP_VERSION+"_"+SP+".tsv.gz" # Only Hsa and Mmu


    ### CHECK existence/make required common data files
    # Pfam.dat
    if os.path.isfile(EDITED_PFAM_DAT_FILE):
        print_status(EDITED_PFAM_DAT_FILE, "exists")
    elif check_file_exists(PFAM_DAT_FILE):
        edit_pfam_dat(PFAM_DAT_FILE, COV2_PFAM_FILE, EDITED_PFAM_DAT_FILE)
        print_status(EDITED_PFAM_DAT_FILE, "created")
    pfam_names, pfam_accs, pfam_types = read_pfam_dat(EDITED_PFAM_DAT_FILE)
    
    # elm classes & domain-ints
    check_file_exists(ELM_CLASSES_FILE)
    elm_map, elm_classes = parse_elm_classes(ELM_CLASSES_FILE)
   
    check_file_exists(ELM_INSTANCES_FILE)

    if os.path.isfile(ELM_DOM_EDIT_FILE):
        print_status(ELM_DOM_EDIT_FILE, "exists")
    elif check_file_exists(ELM_DOM_FILE):
        edit_elm_interactions(ELM_DOM_FILE, ELM_DOM_EDIT_FILE, elm_map)
        print_status(ELM_DOM_EDIT_FILE, "created")

    ints_elm = extract_elm_interactions(ELM_DOM_EDIT_FILE)

    # 3did flat
    if os.path.isfile(DDI_FILE):
        print_status(DDI_FILE, "exists")
    elif check_file_exists(FLAT_3DID_FILE):
        create_ddi_database(pfam_names, PFAM_INT_FILE, FLAT_3DID_FILE, DDI_FILE)
        print_status(DDI_FILE, "created")
   
    # 3did DMI
    check_file_exists(PDB2UNIPROT_FILE)
    check_file_exists(DMI_3DID_FILE)
    edit_lines, lm_3did, all_pdbs = get_dmi_3did(DMI_3DID_FILE, pfam_accs)
    if os.path.isfile(EDITED_DMI_3DID_FILE):
        print_status(EDITED_DMI_3DID_FILE, "exists")
    else:
        with open_file(EDITED_DMI_3DID_FILE, "wt") as out:
            out.write("\t".join(["MOTIF","REGEX","DOMAIN", "DOMAIN_NAME", "PDBS"])+"\n" )
            for line in edit_lines:
                out.write(line+"\n")
        print_status(EDITED_DMI_3DID_FILE, "created")

    check_file_exists(UNI2BIO_FILE)
   
    ## ORGANISM-specific:
    for f in [UNI_TEXT_FILE, UNI_FASTA_FILE, PFAM_MATCHES_FILE]:#, BIOGRID_FILE]:
        check_file_exists(f)
    
    prot_data, uni_id_map, uni_seqs, modres, uni_feats = extract_protein_data_from_uniprot_text(UNI_TEXT_FILE)
    prot_data = add_uniprot_biogrid_mapping(UNI2BIO_FILE, uni_id_map, prot_data)
    bio2uni = get_bio_2_uni_mapping(prot_data)

    # Get PP interaction from BioGRID file
    biogrid_ppi, biogrid_ppi_all, interaction_info = extract_biogrid_interactions(BIOGRID_FILE, bio2uni)    
    
    ### Create PPI database (only BioGRID now)
    if os.path.isfile(PPI_FILE):
        print_status(PPI_FILE, "exists")
    else:
        create_ppi_database(PPI_FILE, biogrid_ppi, interaction_info)
        print_status(PPI_FILE, "created")

    if os.path.isfile(ID_MAP_FILE):
        print_status(ID_MAP_FILE, "exists")
    else:
        create_id_map_table(ID_MAP_FILE, uni_id_map, org)
        print_status(ID_MAP_FILE, "created")

    if os.path.isfile(PROT_DATA_FILE):
        print_status(PROT_DATA_FILE, "exists")
    else:
        prot_data = add_sorted_interactors(prot_data, biogrid_ppi_all, bio2uni)
        create_protein_data_table(PROT_DATA_FILE, prot_data, org)
        print_status(PROT_DATA_FILE, "created")
    
    if os.path.isfile(UNI_FEAT_FILE):
        print_status(UNI_FEAT_FILE, "exists")
    else:
        create_uni_features_table(UNI_FEAT_FILE, uni_feats)
        print_status(UNI_FEAT_FILE, "created")

    if os.path.isfile(PTM_SQL_FILE):
        print_status(PTM_SQL_FILE, "exists")
    else:
        # FALTA PTMS from PhosphoSitePlus
        create_ptms_table(PTM_SQL_FILE, modres, uni_seqs)
        print_status(PTM_SQL_FILE, "created")
 
    # Check for PfamScan
    if os.path.isfile(PFAMSCAN_FILE):
        print_status(PFAMSCAN_FILE, "exists")
    else:
        pfam_matches = extract_pfam_doms(PFAM_MATCHES_FILE)
        run_pfamscan(PFAMSCAN_FILE, uni_seqs, pfam_matches)
        print_status(PFAMSCAN_FILE, "created")

    # Check for Complete Pfam file (matches+pfamscan)
    if os.path.isfile(PFAM_FULL_FILE):
        print_status(PFAM_FULL_FILE, "exists")
    else:
        merge_pfam_files(PFAM_FULL_FILE, PFAM_MATCHES_FILE, PFAMSCAN_FILE)
        print_status(PFAM_FULL_FILE, "created")

    # Check for ELM hits
    if os.path.isfile(ELM_HITS_FILE):
        print_status(ELM_HITS_FILE, "exists")
    else:
        print("Looking for ELMs v."+ELM_VERSION, "in", UNI_FASTA_FILE)
        find_all_slims.main(UNI_FASTA_FILE, elm_classes,
                            print_out=True, outfile=ELM_HITS_FILE, tmp_file="/tmp/lm_search_"+SP+".txt")
        print_status(ELM_HITS_FILE, "created")

    # Check for 3did-motifs hits
    if os.path.isfile(LM_3DID_HITS_FILE):
        print_status(LM_3DID_HITS_FILE, "exists")
    else:
        print("Looking for 3did motifs v."+DB3DID_VERSION, "in", UNI_FASTA_FILE)
        find_all_slims.main(UNI_FASTA_FILE, lm_3did,
                            print_out=True, outfile=LM_3DID_HITS_FILE, tmp_file="/tmp/lm_search_"+SP+".txt")
        print_status(LM_3DID_HITS_FILE, "created")

    # pfam_matches = extract_pfam_doms(PFAM_FULL_FILE)
    # elm_tp, elm_fp = extract_elm_instances(ELM_INSTANCES_FILE, elm_map, uni_seqs.keys())
    # elm_hits = extract_slim_hits(ELM_HITS_FILE, elm_tp, elm_fp)

    # Extract 3DID_LM hits
    # pdb2uni = get_pdb2uniprot_map(PDB2UNIPROT_FILE, uni_seqs.keys())
    # tp = get_3did_dmi_instances(lm_3did, pdb2uni, uni_seqs)
    # lm3d_hits = extract_slim_hits(LM_3DID_HITS_FILE, tp, {})

    if os.path.isfile(LM_HITS_SQL_FILE):
        print_status(LM_HITS_SQL_FILE, "exists")
    else:
        reformat_slim_hit_file(LM_HITS_SQL_FILE , elm_hits, lm3d_hits)
        print_status(LM_HITS_SQL_FILE , "created")
    
    ### Run InterPreTS on PPI
    if os.path.isfile(IPRETS_RAW_FILE):
        print_status(IPRETS_RAW_FILE, "exists")
    # else:
    #     ppi = bio_ppi
    #     if SP=="Hsa":
    #         ppi = hippie_ppi
    #     create_interprets_file(IPRETS_FILE, ppi, prot_dict, SP_DIR+"iprets/")
    #     # os.system("gzip "+IPRETS_FILE)
    #     print_status(IPRETS_FILE, "created")
  
    if os.path.isfile(IPRETS_REV_FILE):
        print_status(IPRETS_REV_FILE, "exists")
    else:
        edit_interprets_raw_file(IPRETS_RAW_FILE, IPRETS_REV_FILE, uni_id_map, uni_seqs, prot_data)  
        print_status(IPRETS_REV_FILE, "created")
   
    sys.exit()
    



    






    # # Extract PTMs positions from PhosphoSite
    # if SP in ["Hsa", "Mmu"]:
    #     check_file_exists(PSP_FILE)
    #     ptms = extract_ptms(PSP_FILE, prot_dict, ptms)

   
    # if SP=="Hsa":
    #     check_file_exists(HIPPIE_MAP_FILE)
    #     check_file_exists(HIPPIE_FILE)
    #     hippie_map = get_HIPPIE_UniProt_map(HIPPIE_MAP_FILE, prot_dict)
    #     hippie_ppi = extract_HIPPIE_interactions(HIPPIE_FILE, hippie_map)


    # ### 4. Make protein data JSON file
    # if os.path.isfile(PROT_DATA_FILE) and force_prot_data==False:
    #     print_status(PROT_DATA_FILE, "exists")
    # else:
    #     print "Creating protein data JSON file"
    #     create_protein_data_json(prot_dict, alt_ids, pfams, elm_hits, lm3d_hits,
    #                            ptms, uni_features, regions, bio_id, bio_set,
    #                            SP, PROT_DATA_FILE, mode="mongo")
    #     print_status(PROT_DATA_FILE, "updated")
    

    # ###  5. Calculate Interactions Probabilities between Protein Elements
    # if os.path.isfile(ASSOC_PROB_FILE) and force_probs==False:
    #     print_status(ASSOC_PROB_FILE, "exists")
    # else:
    #     print "Calculating probabilities"
    #     ppi = bio_ppi
    #     if SP=="Hsa":
    #         ppi = hippie_ppi

    #     eles = merge_protein_elements(pfams, elm_hits, lm3d_hits)
    #     ele_names = merge_protein_elements_names(pfam_names, elm_map, lm_3did)

    #     (all_proteins, total_interactions,
    #         ele_pair_count) = count_protein_and_domain_pairs(ppi, eles)
    #     element_count = count_individual_element_frequency(all_proteins, eles)

    #     ints_ddi = extract_ddi_interactions(DDI_FILE)
    #     ints_dmi = extract_3did_dmi_interactions(EDITED_DMI_3DID_FILE)
    #     print SP, "Total proteins:", len(all_proteins)
    #     print SP, "Total interactions:", total_interactions
    #     create_prob_file(ASSOC_PROB_FILE, len(all_proteins), total_interactions,
    #                      ele_pair_count, element_count, ele_names, ints_ddi,
    #                      ints_elm, ints_dmi)
    #     print_status(ASSOC_PROB_FILE, "created")



    # print "Finished without problems :-)"
    # sys.exit()



    ## Preprocessing InterPreTS info
    # chain2uni = {}
    ## WRONG
    ## I cant use fixed PDB-UniAcc conversion because InterPreTS uses homology
    ## to find templates!
    ## Don't do anything for the moment
    # s = set()
    # with open_file(pdbchain2uniprot) as f:
    #     for line in f:
    #         t = line.upper().rstrip().split()
    #         if len(t) == 3:
    #             pdb = t[0]
    #             chain = t[1]
    #             uni_ac = t[2]
    #             if uni_ac in prot_dict["AC"]:
    #                 for ac in prot_dict["AC"][uni_ac]:
    #                     chain2uni["pdb|"+pdb+"|"+chain].add(ac)

    # ints = defaultdict(lambda: defaultdict(list))
    # with open_file(IPRETS_FILE) as f: ## lines in this file are already sorted from highest-to-lowest Z-score!!
    #     for line in f:
    #         t = line.rstrip().split("\t")
    #         if line[0] == "#":
    #             columns = t
    #         else:
                # pdb1 = t[1]
                # pdb2 = t[9]
                # if pdb1 not in chain2uni:
                #     s.add(pdb1)
                # if pdb2 not in chain2uni:
                #     s.add(pdb2)
                # if pdb1 in chain2uni and pdb2 in chain2uni:
                #     for uni_ac1 in chain2uni[pdb1]:
                #         for uni_ac2 in chain2uni[pdb2]:
                #             print t[0]+"\t"+uni_ac1+"\t"+"\t".join(t[1:9])+"\t"+uni_ac2+"\t"+"\t".join(t[9:])

                # gene1, gene2 = t[0], t[8]
                # if len(ints[gene1][gene2])>0:
                #     compare_ints(t, gene1, gene2, ints)
                # else:
                #     ints[gene1][gene2].append("\t".join(t))


    # print  "# not found:", len(s)
    # for pdb in s:
    #     print pdb



if __name__ == "__main__":
    pdata, ppi, probs, iprets = False, False, False, False
    if "-prot_data" in sys.argv:
        pdata = True
    if "-ppi" in sys.argv:
        ppi = True
    if "-probs" in sys.argv:
        probs = True
    if "-iprets" in sys.argv:
        iprets = True
    main(SP=sys.argv[1], DATA_DIR="/net/home.isilon/ag-russell/bq_jgonzalez/int2mech/data/",
         force_prot_data=pdata, force_ppi=ppi, force_probs=probs,
         force_iprets=iprets)
    sys.exit()
