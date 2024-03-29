#!/usr/bin/env python

import os, re, sys, gzip, random, string, itertools, datetime
from collections import defaultdict
from Bio import SearchIO, SeqIO


def open_file(input_file, mode="r"):
    """ Open file Zipped or not
    """
    if re.search(".gz$", input_file):
        infile = gzip.open(input_file, mode)
    else:
        infile = open(input_file, mode)
    return infile

def make_dir(directory, verbose):
    if not os.path.exists(directory):
        try:
            os.mkdir(directory)
            if verbose:
                print( directory, "created")
        except OSError:
            print("Creation of the directory {} failed".format(directory))
    return

def run_blast(mode, psiblast, blastpgp, blastmat, db,
              seq_file, blast_out_file, verbose=False):
    
    if mode=="blastpgp": ## Not working with apache
        #    -m  alignment view options
        #    -d database
        #    -b  Number of database sequence to show alignments for (B) [Integer]
        #         default = 250
        #    -v  Number of database sequences to show one-line descriptions for (V) [Integer]
        #         default = 500
        options = "-m 0 -b 1000 -v 1000"
        if ".gz" in seq_file:
            com = "zcat {} | {} -d {} {}".format(seq_file, blastpgp, db, options)
        else:
            com = "{} -i {} -d {} {}".format(blastpgp, seq_file, db, options)
        os.environ["BLASTMAT"] = blastmat
    elif mode=="psiblast":
        options = "-outfmt 5 -num_iterations 1 -max_target_seqs 1000"
        if ".gz" in seq_file:
            com = "zcat {} | {} -db {} {}".format(seq_file, psiblast, db, options)
        else:
            com = "{} -query {} -db {} {}".format(psiblast, seq_file, db, options)
    else:
        print("ERROR: could not recognize Blast mode: \'"+mode+"\'")
        return

    if ".gz" in blast_out_file:
        com += " | gzip > "+blast_out_file
    else:
        com += " > "+blast_out_file
    # print "[{}] Running BLAST: {}".format(datetime.datetime.now(), com)
    if verbose:
        print(com)
    os.system(com)
    
def parse_blast(blast_pdb_file, max_E, min_pcid, max_pcid, hits, mode="blast-xml"):
    with open_file(blast_pdb_file, "rt") as f:
        for qresult in SearchIO.parse(f, mode):
            query = qresult.id#.split("|")[1]
            for hit in qresult:
                s = hit.id + hit.description
                hsp = hit[0] # Only the 1st one
                evalue = hsp.evalue
                pcid = float(hsp.ident_num)/hsp.aln_span*100
                if (evalue<=max_E
                and pcid>=min_pcid and pcid<=max_pcid):
                    # print "\t>HIT:",hit.id, set(re.findall("pdb\|\w\w\w\w\|\w", s))
                    # print "\t", hsp.evalue, "{:2.1f}".format(pcid)
                    # print hsp.query_start, hsp.query_end
                    # print hsp.hit_start+1, hsp.hit_end
                    for match in re.findall("pdb\|\w\w\w\w\|\w", s):
                        pdb, chain = match.split("|")[1:]
                        hits[query][pdb][chain]={
                                      "ide": "{:2.1f}".format(pcid),
                                      "e-val": evalue,
                                      "q-start": str(hsp.query_start+1),
                                      "q-end": str(hsp.query_end),
                                      "s-start": str(hsp.hit_start+1),
                                      "s-end": str(hsp.hit_end)
                                      }
                else:
                    break
    return hits

def main(input_seqs, output_dir, i2sum_file, ide, org_map,
    fasta_as_input=False, print_output=True, verbose=False, return_hits=False,
    force_new=False,
    these_pairs=[], hits=defaultdict(lambda: defaultdict(dict)),
    max_templates=10,
    mode="blastpgp", # Like Rob's original script
    psiblast="psiblast",
    blastpgp="/net/home.isilon/ag-russell/install/CentOS-7.3.1611-x86_64/bin/blastpgp",
    blastmat="/net/home.isilon/ag-russell/install/CentOS-7.3.1611-x86_64/blast-2.2.23/data",
    blastdb="/net/home.isilon/ds-russell/blastdb/pdbaa_2019",
    muscle="/net/home.isilon/ag-russell/install/CentOS-7.3.1611-x86_64/bin/muscle -quiet",
    i2="/net/home.isilon/ag-russell/install/CentOS-7.3.1611-x86_64/bin/interprets",
    i2_opts=" -rand 500 -show_muts -mode 4 -q",
    data_dir="static/data/"
    ):
    os.environ["I2DIR"] = "/net/home.isilon/ag-russell/code/interprets/data"
    mode = mode+"_v_"+blastdb.split("/")[-1]
    
    results = {}

    ### 1: Get input sequences (either dictionary or file)
    if fasta_as_input:
        fasta_file = input_seqs
        seq = {}
        with open_file(fasta_file) as f:
            for record in SeqIO.parse(f, "fasta"):
                if "|" in record.id:
                    seq[str(record.id).split("|")[1]]=str(record.seq)
                else:
                    seq[str(record.id)]=str(record.seq)
    else:
        seq = input_seqs
    if verbose:
        print("Received sequences for {} proteins".format(len(seq)))


    ### 2: Create neccessary directories
    seqs_dir = output_dir+"seqs/"
    temp_dir = output_dir+"temp/"
    blast_dir = output_dir+mode+"/"
    for directory in [output_dir, seqs_dir, blast_dir, temp_dir]:
        make_dir(directory, verbose)
    for org in set(org_map.values()):
        data_seqs_dir = data_dir+org+"/seqs/"
        data_blast_dir = data_dir+org+"/"+mode+"/"
        make_dir(data_seqs_dir, verbose)
        make_dir(data_blast_dir, verbose)

    ### 3: Print individual fasta
    seq_file = {}
    for prot_id in seq:
        file_name = prot_id.replace("|","_")+".fa"
        data_seqs_dir = data_dir+org_map[prot_id]+"/seqs/"
        old_file = data_seqs_dir+file_name
        new_file = seqs_dir+file_name
        if os.path.isfile(old_file) and not force_new:
            seq_file[prot_id] = old_file
        elif not os.path.isfile(new_file):
            with open_file(new_file, "w") as out:
                out.write(">"+prot_id+"\n")
                out.write(seq[prot_id]+"\n")
            os.system("cp "+new_file+" "+data_seqs_dir)
            seq_file[prot_id] = new_file
        
    ### 4: Run individual blastpgp
    blast_files = []
    for prot_id in seq:
        data_blast_dir = data_dir+org_map[prot_id]+"/"+mode+"/"
        if mode.startswith("psiblast"):
            file_name = prot_id.replace("|","_")+"_psiblast.xml.gz"
            blast = "psiblast"
        elif mode.startswith("blastpgp"):
            file_name = prot_id.replace("|","_")+"_blastpgp.txt.gz"
            blast = "blastpgp"
        old_file = data_blast_dir+file_name
        new_file = blast_dir+file_name
        if os.path.isfile(old_file) and not force_new:
            blast_files.append(old_file)
        elif not os.path.isfile(new_file):
            run_blast(blast, psiblast,  blastpgp, blastmat, blastdb,
                      seq_file[prot_id], new_file, verbose=verbose)
            os.system("cp "+new_file+" "+data_blast_dir)
            blast_files.append(new_file)
    if verbose:
        print("[{}] BLAST done".format(datetime.datetime.now()))

    ### 5: Read BLAST
    max_E = 0.01
    min_pcid, max_pcid = -1, 1000
    # hits = defaultdict(lambda: defaultdict(dict))
    for blast_out_file in blast_files:
        if mode.startswith("psiblast"):
            hits = parse_blast(blast_out_file, max_E, min_pcid, max_pcid,
                    hits, mode="blast-xml")
        elif mode.startswith("blastpgp"):
            hits = parse_blast(blast_out_file, max_E, min_pcid, max_pcid, 
                    hits, mode="blast-text")
    if verbose:
        print("[{}] BLAST parsed for:".format(datetime.datetime.now()), str(len(hits.keys())))

    ### 6: Find best hits for each pair
    pairs = defaultdict(lambda: defaultdict(list))
    if not these_pairs:
        these_pairs = sorted(itertools.combinations(hits.keys(), 2))


    for pair in these_pairs:
        q1, q2 = sorted(pair)
        results[(q1, q2)] = defaultdict(list)

        if q1 not in hits or q2 not in hits:
            continue
  
        for pdb in sorted(list(set(hits[q1]).intersection(set(hits[q2])))):
            for c1 in sorted(hits[q1][pdb]):
                hit1 = hits[q1][pdb][c1]
                for c2 in sorted(hits[q2][pdb]):
                    hit2 = hits[q2][pdb][c2]
                    if c1 != c2:
                        eval1 = hit1["e-val"]
                        eval2 = hit2["e-val"]
                        avg_eval = float(eval1)+float(eval2)/2
                        ele1 = "\t".join([q1,"pdb|"+pdb+"|"+c1,str(eval1),hit1["ide"]])
                        ele1 += "\t"+"\t".join([hit1["q-start"], hit1["q-end"]])
                        ele1 += "\t"+"\t".join([hit1["s-start"], hit1["s-end"]])
                        ele2 = "\t".join([q2,"pdb|"+pdb+"|"+c2,str(eval2),hit2["ide"]])
                        ele2 += "\t"+"\t".join([hit2["q-start"], hit2["q-end"]])
                        ele2 += "\t"+"\t".join([hit2["s-start"], hit2["s-end"]])
                        pairs[(q1, q2)][avg_eval].append( (pdb+":"+c1+":"+c2, ele1+"\t"+ele2) )

    ### 7: Run InterPreTS
    already = set()
    i = 0
    for pair in pairs:
        (q1, q2) = pair
        n_template = 0
        for avg in sorted(pairs[pair]):
            for match in sorted(pairs[pair][avg]):
                if n_template >= max_templates:
                    break
                i+=1
                pdb, c1, c2 = match[0].split(":")

                ## Files names
                dfile = temp_dir+ide+"_"+str(i)+"_"+pdb+".dom"
                cfile = temp_dir+ide+"_"+str(i)+"_"+pdb+".contact"
                # cfile = temp_dir+"see"
                f1 = temp_dir+ide+"_"+str(i)+"_"+pdb+"_0.fa"
                f2 = temp_dir+ide+"_"+str(i)+"_"+pdb+"_1.fa"
                m1 = temp_dir+ide+"_"+str(i)+"_"+pdb+"_0_muscle.txt"
                m2 = temp_dir+ide+"_"+str(i)+"_"+pdb+"_1_muscle.txt"
                ras = temp_dir+ide+"_"+str(i)+"_"+pdb+"_i2_ras.txt"
                fi2 = temp_dir+ide+"_"+str(i)+"_"+pdb+"_i2.txt"

                ## Make domain file
                with open_file(dfile, "w") as out:
                    out.write("UNK "+pdb+c1+" { CHAIN "+c1+" }\n")
                    out.write("UNK "+pdb+c2+" { CHAIN "+c2+" }\n")

                ## Test contacts
                com = i2+" -d "+dfile+" -count -q > "+cfile
                if verbose:
                    print("[{}] Testing contacts: {}".format(
                                                datetime.datetime.now(), com))
                os.system(com)
                n_contact = 0
                with open(cfile) as f:
                    for line in f:
                        t = line.rstrip().split()
                        r1 = t[3]+t[4]
                        r2 = t[12]+t[13]
                        if r1+"-"+r2 in already:
                            continue
                        if ((t[3]==c1 and t[12]==c2)
                        or  (t[3]==c2 and t[12]==c1)):
                            n_contact += 1

                        already.add(r1+"-"+r2)
                        already.add(r2+"-"+r1)

                if n_contact>5:
                    # 1. Get protein sequences
                    ide = ''.join(random.choice(string.ascii_uppercase +
                                    string.ascii_lowercase + string.digits) for _ in range(8))
                    tmp_fasta = "/tmp/tmp_"+ide+".fasta"
                    com = i2+" -fasta -q -d "+dfile+" > "+tmp_fasta
                    if verbose:
                        print("[{}] {}".format(datetime.datetime.now(), com))
                    os.system(com)
                    with open_file(tmp_fasta) as f:
                        for record in SeqIO.parse(f, "fasta"):
                            seq[str(record.id)] = str(record.seq)

                    # 2. Make pairwise files for alignment
                    with open_file(f1, "w") as out:
                        for x in [q1, pdb+c1]:
                            out.write(">"+x+"\n")
                            out.write(seq[x]+"\n")
                    with open_file(f2, "w") as out:
                        for x in [q2, pdb+c2]:
                            out.write(">"+x+"\n")
                            out.write(seq[x]+"\n")

                    # 3. Run muscle
                    com = muscle+" -in "+f1+" -clwstrict > "+m1
                    if verbose:
                        print("[{}] Align: {}".format(datetime.datetime.now(), com))
                    os.system(com)
                    com = muscle+" -in "+f2+" -clwstrict > "+m2
                    if verbose:
                        print("[{}] Align: {}".format(datetime.datetime.now(), com))
                    os.system(com)

                    # 4. Run interprets
                    com = i2+" -d "+dfile+" -a "+m1+" "+m2+" "+i2_opts+" -op \'"+q1+"\' \'"+q2+"\' > "+fi2
                    if verbose:
                        print("[{}] Final InterPreTS: {}".format(datetime.datetime.now(), com))
                    os.system(com)

                    # 5. Read i2 file
                    i2_sum = ""
                    with open_file(fi2) as f:
                        for line in f:
                            if line.startswith("SUM"):
                                t = line.rstrip().split()
                                i2_sum = "\t".join(t[5:])
                                z = float(t[-4])

                    # if n_template == 0: # only the first (and best) hit per pair
                    results[(q1, q2)][z].append(
                        { "info_a": match[1].split("\t")[1:8],
                          "info_b": match[1].split("\t")[9:],
                          "scores": i2_sum.split("\t")
                    })

                    if print_output==True:
                        with open_file(i2sum_file, "at") as out2:
                            l = [match[1], i2_sum]
                            out2.write("\t".join(l)+"\n")

                    n_template += 1

                    for fl in [dfile, cfile, f1, f2, m1, m2, fi2]:
                        try:
                            os.unlink(fl)
                        except:
                            print(f1,"wasn't created")
                else:
                    for fl in [dfile, cfile]:
                        os.unlink(fl)
    if return_hits:
        return results, hits
    else:
        return results

if __name__ == "__main__":
    fasta_file = sys.argv[1]
    output_dir = sys.argv[2]
    i2sum_file = "final_output.tsv.gz"
    main( fasta_file, output_dir, i2sum_file,
          fasta_as_input=True)
