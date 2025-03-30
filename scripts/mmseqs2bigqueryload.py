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
                    description='Convert MMseq2 clusters to table for orthology',
                    epilog='Example: mmseqs2bigqueryload.py [-i input_clusters.tsv] [-o bigquery/mmseqs_clusters.csv]')
    parser.add_argument('-i','--input', help='Input MMSeqs cluster file', nargs='?', 
                        type=argparse.FileType('r'),
                        default=sys.stdin)
    
    parser.add_argument('-o', '--output', help='Output gzip file for bigquery load', 
                        required=False, default='bigquery/mmseqs_orthogroup_clusters.csv.gz')
    parser.add_argument('--prefix', help='Cluster Name prefix', default="MMCLUST")
    parser.add_argument('-v','--debug', help='Debugging output', action='store_true')
    parser.add_argument('--force', help='Force overwriting', action='store_true')
    args = parser.parse_args()
    if os.path.exists(args.output) and not os.path.exists(args.force):
        print(f"Output file {args.output} already exists. Exiting.", file=sys.stderr)
        sys.exit(1)
    with gzip.open(args.output, 'wt') as out:
        # Read in the cluster file
        print(f"orthogroup,transcript_id", file=out)
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
            if args.debug and i > 0 and i % 1000000 == 0:
                    t6 = time.time()
                    print(f"Processing {i} lines took {t6-t5} seconds; ClusterID is {clusterID}", file=sys.stderr)
                    t5 = t6
            i+=1
            (gene1,gene2) = cluster
            if gene1 != last_seq and last_seq != "":
                clusterID += 1
            print(f"{args.prefix}{clusterID:0>8},{gene2}", file=out)
            last_seq = gene1

if __name__ == "__main__":
    main()
