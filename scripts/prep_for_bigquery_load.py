#!/usr/bin/env python3

import os
import csv
import gzip

outdir = 'bigquery'
def merops(indir="results/function/merops",force=False):
    # load MEROPS data
    outfile = os.path.join(outdir,os.path.basename(indir) + ".csv")
    if os.path.exists(outfile) and not force:
        return
    with open(outfile, "w", newline='') as of:
        writer = csv.writer(of)
        writer.writerow(['species_prefix','protein_id','merops_id','percent_identity','aln_length','mismatches','gap_openings','q_start','q_end',
                        's_start','s_end', 'evalue', 'bitscore'])
        for file in os.listdir(indir):
            if file.endswith(".blasttab.gz"):
                with gzip.open(os.path.join(indir,file), "rt") as infh:
                    reader = csv.reader(infh, delimiter='\t')
                    for row in reader:
                        prefix = row[0].split('_')[0]
                        newrow = [prefix]
                        newrow.extend(row)
                        writer.writerow(newrow)

def cazy_overview(indir="results/function/cazy",force=False):
    # load CAZY data
    outfile = os.path.join(outdir,os.path.basename(indir) + ".overview.csv")
    if os.path.exists(outfile) and not force:
        return
    with open(outfile, "w", newline='') as of:
        writer = csv.writer(of)
        # Gene_ID	EC	cazyme_fam	sub_fam	diamond_fam	Substrate	#ofTools
        writer.writerow(['species_prefix','protein_id','EC','cazyme_fam','sub_fam','diamond_fam','substrate','toolcount'])
        for spdir in os.listdir(indir):
            infile = os.path.join(indir,spdir,'overview.tsv.gz')
            if os.path.exists(infile):
                with gzip.open(infile, "rt") as infh:
                    reader = csv.reader(infh, delimiter='\t')
                    next(reader)
                    for row in reader:
                        prefix = row[0].split('_')[0]
                        newrow = [prefix]
                        newrow.extend(row)
                        writer.writerow(newrow)

def cazy_hmm(indir="results/function/cazy",force=False):
    # load CAZY data
    outfile = os.path.join(outdir,os.path.basename(indir) + ".cazymes_hmm.csv")
    if os.path.exists(outfile) and not force:
        return
    with open(outfile, "w", newline='') as of:
        writer = csv.writer(of)
        # HMM_Profile	Profile_Length	Gene_ID	Gene_Length	Evalue	Profile_Start	Profile_End	Gene_Start	Gene_End	Coverage
        writer.writerow(['species_prefix','HMM_id','profile_length','protein_id','protein_length','evalue',
                        'q_start','q_end','s_start','s_end', 'coverage'])
        for spdir in os.listdir(indir):
            infile = os.path.join(indir,spdir,'cazymes.tsv.gz')
            if os.path.exists(infile):
                with gzip.open(infile, "rt") as infh:
                    reader = csv.reader(infh, delimiter='\t')
                    next(reader)
                    for row in reader:
                        prefix = row[2].split('_')[0]
                        newrow = [prefix]
                        newrow.extend(row)
                        writer.writerow(newrow)

def signalp(indir="results/function/signalp",force=False):
    # load signalp data
    outfile = os.path.join(outdir,os.path.basename(indir) + ".signal_peptide.csv")
    if os.path.exists(outfile) and not force:
        return
    with open(outfile, "w", newline='') as of:
        writer = csv.writer(of)
        # HMM_Profile	Profile_Length	Gene_ID	Gene_Length	Evalue	Profile_Start	Profile_End	Gene_Start	Gene_End	Coverage
        writer.writerow(['species_prefix','protein_id','peptide_start','peptide_end','probability'])
        for file in os.listdir(indir):
            if file.endswith(".signalp.gff3.gz"):
                with gzip.open(os.path.join(indir,file), "rt") as infh:
                    reader = csv.reader(infh, delimiter='\t')
                    for row in reader:
                        if row[0].startswith('#'):
                            continue
                        id = row[0].split(' ')[0]
                        prefix = id.split('_')[0]                        
                        newrow = [prefix, id, row[3], row[4], row[5]]
                        writer.writerow(newrow)

merops()
cazy_overview()
cazy_hmm()
signalp()