#!/usr/bin/env python3

import sys
import argparse
import pathlib
import re
import os
import time
import gzip
import csv

def main():
    parser = argparse.ArgumentParser(
                    prog='mmseqs2pairwise.py',
                    description='Convert MMseq2 clusters to pairwise orthologs',
                    epilog='Example: mmseqs2pairwise.py [-i input_clusters.tsv] -d outdir')
    parser.add_argument('-i','--input', help='Input MMSeqs cluster file', nargs='?', 
                        type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-o', '--output', help='Output gzip ortholog table like orthofinder', 
                        required=False, default='results/orthogroups_mmseqs.tsv.gz')
    parser.add_argument('-ot', '--output_table', help='Output gzip file for bigquery load', 
                        required=False, default='results/orthogroups_mmseqs_genecount.tsv')
    parser.add_argument('--prefix', help='Cluster Name prefix', default="MMCLUST")
    parser.add_argument('-v','--debug', help='Debugging output', action='store_true')
    parser.add_argument('--force', help='Force overwriting', action='store_true')
    parser.add_argument('-s','--samples', type = str, default = 'samples.csv', help = 'species prefix')

    args = parser.parse_args()

    samples = {}
    species2locus = {}
    species_by_locus = []
    header = []
    with open(args.samples, 'r') as fh:
        sampleinfo = csv.DictReader(fh, delimiter=",")
        for row in sampleinfo:
            samples[row['LOCUSTAG']] = row['SPECIES'].replace(' ', '_')
            if len(row['STRAIN']):
                samples[row['LOCUSTAG']] += f"_{row['STRAIN']}"
            species2locus[ samples[row['LOCUSTAG']] ] = row['LOCUSTAG']
            species_by_locus.append(row['LOCUSTAG'])
            header.append(samples[row['LOCUSTAG']])
        
    # Read in the cluster file
    clusterID = 0
    last_seq = ""
    last_cluster = set()
    i = 0
    t5 = time.time()
    t4 = t5
    with gzip.open(args.output, 'wt') as out, open(args.output_table, 'w') as out_og_table:
        out.write("\t".join(['Orthogroup']+ header) + "\n")
        out_og_table.write("\t".join(['Orthogroup']+ header + ['Total']) + "\n")
        for line in args.input:
            if line.startswith('#'): continue   # speedup if we assume no comment lines?
            cluster = line.strip().split()
            if len(cluster) != 2: continue
            if args.debug and i % 1000000 == 0:
                    t6 = time.time()
                    print(f"Processing line {i} took {t6-t5} seconds; ClusterID is {clusterID}", file=sys.stderr)
                    t5 = t6
            i+=1
            for gene in cluster:
                LOCUSTAG = gene.split("_")[0]
                if LOCUSTAG not in samples:
                    print(f"Gene {gene} Prefix {LOCUSTAG} not found in {args.samples} file", file=sys.stderr)
                    continue
            (gene1,gene2) = cluster
            if gene1 != last_seq and last_seq != "":
                species_grouping = {}
                row = [ f'{args.prefix}{clusterID:0>8}' ]
                countrow = [ f'{args.prefix}{clusterID:0>8}' ]
                for item in last_cluster:
                    locustag = item.split('_')[0]                    
                    if locustag not in species_grouping:
                        species_grouping[locustag] = []
                    species_grouping[locustag].append(item)                
                total = 0
                for locustag in species_by_locus:
                    if locustag not in species_grouping:
                        row.append("")
                        countrow.append("0")
                    else:
                        row.append(", ".join(species_grouping[locustag]))
                        total += len(species_grouping[locustag])
                        countrow.append(f'{len(species_grouping[locustag])}')
                countrow.append(str(total))
                out.write("\t".join(row) + "\n")
                out_og_table.write("\t".join(countrow) + "\n")
                clusterID += 1
                last_cluster = set()
            last_cluster.add(gene1)
            last_cluster.add(gene2)
            last_seq = gene1

        # fence post
        species_grouping = {}
        row = [ f'{args.prefix}{clusterID:0>8}' ]
        countrow = [ f'{args.prefix}{clusterID:0>8}' ]
        for item in last_cluster:
            species = gene2species[item]
            if species not in species_grouping:
                species_grouping[species] = {}
            species_grouping[species].append(item)                
            for species in species_by_locus:
                if species not in species_grouping:
                    row.append("")
                    countrow.append(0)
                else:
                    row.append(", ".join(species_grouping[species]))
                    countrow.append(len(species_grouping[species]))
        
        if args.debug:
            t2 = time.time()
            total = t2-t4
            print(f"Processing clusters took {total} seconds",file=sys.stderr)

if __name__ == "__main__":
    main()
