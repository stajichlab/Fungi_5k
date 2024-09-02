#!/usr/bin/env python3
import csv
import sys
import argparse
import pathlib
import re
import os
import time
import gzip
import bsddb3

def main():
    parser = argparse.ArgumentParser(
                    prog='mmseqs2pairwise.py',
                    description='Convert MMseq2 clusters to pairwise orthologs',
                    epilog='Example: mmseqs2pairwise.py [-i input_clusters.tsv] -d outdir')
    parser.add_argument('-i','--input', help='Input MMSeqs cluster file', nargs='?', 
                        type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-o', '--outdir', help='Output directory for pairwise', 
                        required=False, default='pairwise_MMseq_orthologs')
    parser.add_argument('--seqs','--fasta', help='Directory to link gene names to species', 
                        type=pathlib.Path, required=True)
    parser.add_argument('--prefix', help='Cluster Name prefix', default="MMCLUST")
    parser.add_argument('--tmpdir','--tmp', help='Temporary Directory for local DB', default="./tmp")
    parser.add_argument('-v','--debug', help='Debugging output', action='store_true')
    
    args = parser.parse_args()
    
    if os.path.exists(args.outdir):
        print(f"Output directory {args.outdir} already exists. Will overwrite existing files.")
    else:
        os.mkdir(args.outdir)
    if not os.path.exists(args.tmpdir):
        os.mkdir(args.tmpdir)

    # open directory and read seq files
    gene2species = {}
    species_gene_count = {}
    filepat = re.compile(r'(\S+)\.(fa|faa|fasta|seq|aa|pep|cds|nt|dna)')    
    t0 = time.time()
    for file in os.listdir(args.seqs):
        fp = filepat.search(file)
        if fp:
            species = fp.group(1)
            species = species.replace('.final','')   # deal with extra suffixes
            species = species.replace('.proteins','')   # deal with extra suffixes
            # assume this will be okay match
            # t3 = time.time()
            with open(os.path.join(args.seqs,file),'r') as f:
                for line in f:                    
                    if line.startswith('>'):
                        gene = line.split()[0][1:]
                        gene2species[gene] = species
                        species_gene_count[species] = species_gene_count.get(species,0) + 1
            # t4 = time.time()
            # if args.debug:
            #    print(f"Reading {file} took {t4-t3} seconds",file=sys.stderr)
    
    if args.debug:
        t1 = time.time()
        total = t1-t0
        print(f"Reading gene files took {total} seconds",file=sys.stderr)
    
    # Read in the cluster file
    groups = []
    pairwise = {}
    clusterID = 0
    last_seq = ""
    last_cluster = set()
    for line in args.input:
        if line.startswith('#'): continue
        cluster = line.strip().split()
        if len(cluster) != 2: continue
        for gene in cluster:
            if gene not in gene2species:
                print(f"Gene {gene} not found in species list",file=sys.stderr)
                continue
        gene1 = cluster[0]
        gene2 = cluster[1]
        if gene1 != last_seq and last_seq != "":
            if len(last_cluster) > 1:
                #print(f"Cluster {clusterID}:")
                #print(last_cluster)
                names = []
                species_grouping = {}
                for item in last_cluster:
                    species = gene2species[item]
                    if species not in species_grouping:
                        species_grouping[species] = []
                    species_grouping[species].append(item)                
                
                for species1 in species_grouping:
                    if species1 not in pairwise:
                        pairwise[species1] = bsddb3.btopen(os.path.join(args.tmpdir,f'{species1}.bdb'))
                    sp1_orthologs = ", ".join(species_grouping[species1])                    
                    for species2 in species_grouping:
                        if species1 == species2: continue
                        sp2_orthologs = ", ".join(species_grouping[species2])    
                        pairwise[species1].append([f'{args.prefix}{clusterID:0>8}',species2,
                                                sp1_orthologs,sp2_orthologs])
                groups.append(last_cluster)
                clusterID += 1
            last_cluster = set()
        last_cluster.add(gene1)
        last_cluster.add(gene2)
        last_seq = gene1

    print(f"{len(groups)} clusters found",file=sys.stderr)
    for sp in pairwise:
        with open(f"{args.outdir}/{sp}_orthologs.tsv",'w') as f:
            writer = csv.writer(f,delimiter='\t',lineterminator=os.linesep)
            writer.writerow(['Cluster','Species',sp,'Orthologs'])
            for row in pairwise[sp]:
                writer.writerow(row)
        

if __name__ == "__main__":
    main()
