#!/usr/bin/env python3

import sys
import argparse
import pathlib
import re
import os
import time
import gzip
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
    parser.add_argument('-v','--debug', help='Debugging output', action='store_true')
    
    args = parser.parse_args()
    
    if os.path.exists(args.outdir):
        print(f"Output directory {args.outdir} already exists. Will overwrite existing files.",file=sys.stderr)
    else:
        os.mkdir(args.outdir)
    gene2species = {}
    #species_gene_count = {}

    # Read in the gene files and populate the database
    filepat = re.compile(r'(\S+)\.(fa|faa|fasta|seq|aa|pep|cds|nt|dna)')    
    t0 = time.time()
    for file in os.listdir(args.seqs):
        fp = filepat.search(file)
        if fp:
            species = fp.group(1)
            species = species.replace('.final','')   # deal with extra suffixes
            species = species.replace('.proteins','')   # deal with extra suffixes
            with open(os.path.join(args.seqs,file),'r') as f:
                for line in f:                    
                    if line.startswith('>'):
                        gene = line.split()[0][1:]
                        gene2species[gene] = species
                        # species_gene_count[species] = species_gene_count.get(species,0) + 1
    
    if args.debug:
        t1 = time.time()
        total = t1-t0
        print(f"Reading gene files took {total} seconds",file=sys.stderr)
    
    # Read in the cluster file
    pairwise = {}
    clusterID = 0
    last_seq = ""
    last_cluster = set()
    i = 0
    t5 = time.time()
    t4 = t5
    for line in args.input:
        if line.startswith('#'): continue   # speedup if we assume no comment lines?
        cluster = line.strip().split()
        if len(cluster) != 2: continue
        if args.debug and i % 1000000 == 0:
                t6 = time.time()
                print(f"Processing line {i} took {t6-t5} seconds; ClusterID is {clusterID}", file=sys.stderr)
                t5 = t6
        i+=1
        #if args.debug and i > 50000:
        #    break
        for gene in cluster:
            if gene not in gene2species:
                print(f"Gene {gene} not found in species list", file=sys.stderr)
                continue
        (gene1,gene2) = cluster
        if gene1 != last_seq and last_seq != "":
            if len(last_cluster) > 1:
                species_grouping = {}
                for item in last_cluster:
                    species = gene2species[item]                    
                    if species not in species_grouping:
                        species_grouping[species] = []
                    species_grouping[species].append(item)
                
                for species1 in species_grouping:
                    if species1 not in pairwise:
                        fh = gzip.open(f"{args.outdir}/{species1}_orthologs.tsv.gz",'wt')
                        #fh = open(f"{args.outdir}/{species1}_orthologs.tsv",'wt')
                        if fh:
                            #writer = csv.writer(fh,delimiter='\t',lineterminator=os.linesep)
                            print("\t".join(['Cluster','Species',species1,'Orthologs']), file=fh)
                            #writer.writerow(['Cluster','Species',species1,'Orthologs'])
                            pairwise[species1] = fh
                            #pairwise[species1] = {'handle': fh, 'csv': writer }
                        else:
                            print(f"Could not open {args.outdir}/{species1}_orthologs.tsv",file=sys.stderr)
                            sys.exit(1)
                    sp1_orthologs = ", ".join(species_grouping[species1])
                    for species2 in species_grouping:
                        if species1 == species2: continue
                        sp2_orthologs = ", ".join(species_grouping[species2])
                        print("\t".join([f'{args.prefix}{clusterID:0>8}',species2,sp1_orthologs,
                                        sp2_orthologs]), file=pairwise[species1])                
                        #pairwise[species1]['csv'].writerow([f'{args.prefix}{clusterID:0>8}',
                        #species2,sp1_orthologs,
                        #                                    sp2_orthologs])                
                clusterID += 1

            last_cluster = set()
        last_cluster.add(gene1)
        last_cluster.add(gene2)
        last_seq = gene1

    if args.debug:
        t2 = time.time()
        total = t2-t4
        print(f"Processing clusters took {total} seconds",file=sys.stderr)
    for d in pairwise.values():
        #d['handle'].close()
        d.close()

if __name__ == "__main__":
    #cProfile.run('main()')
    main()
