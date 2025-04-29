#!/usr/bin/env python3

# Process AIupred results to call IDPs
# Definition of IDP region will be a stretch 
# A region that have at least 30 consecutive residues with iupred score over 
# 0.5 would be determined as putative disordered region (IDR).
# Any protein that have at least 1 IDR would be classified as IDP.
import os
import csv
import sys
import re
import time
import gzip
import argparse

def average(lst):
    """
    Calculate the average of a list of numbers.
    """
    if len(lst) == 0:
        return 0
    return sum(lst) / len(lst)


def scores_to_idp_regions(iupred_scores, min_length=30):
    """
    Convert IUPred scores to IDP regions.
    """
    idp_regions = []
    current_start = None
    for i, score in enumerate(iupred_scores):
        if score >= 0.5:
            if current_start is None:
                current_start = i
        else:
            if current_start is not None:
                s = current_start
                e = i - 1
                mean_score = average(iupred_scores[s:e+1])
                # 1 based not zero based
                idp_regions.append((s+1, e+1, mean_score))
                current_start = None
    if current_start is not None:
        s = current_start
        e = len(iupred_scores) - 1
        mean_score = average(iupred_scores[s:e+1])
        # 1 based not 0 based
        idp_regions.append((s+1, e+1, mean_score))
    return idp_regions

def parse_iupred_file(iupred_file):
    """
    Parse IUPred file and return a list of scores.
    """
    iupred_scores = []
    
    with gzip.open(iupred_file, 'rt') as f:
        seqname = None
        iupred_scores = []
        score_set = {}
        for line in f:
            if line.startswith("#>"):
                if seqname is not None:
                    score_set[seqname] = scores_to_idp_regions(iupred_scores)
                    iupred_scores = []
                seqname = line[2:].split()[1]
            elif line.startswith("#"):
                # skip header/comment
                continue
            else:
                parts = line.strip().split()
                if len(parts) < 3:
                    continue
                try:
                    score = float(parts[2])
                    iupred_scores.append(score)
                except ValueError:
                    continue
    return score_set

def main():
    if len(sys.argv) < 2:
        print("Usage: python gather_AIUPred.py iupred.ouput.txt.gz or -d <dir>")
        sys.exit(1)
        
    parser = argparse.ArgumentParser(description="Collect precomputed AIUPred results into a table",
                                    epilog='Example: gather_AIUPred.py')
    parser.add_argument("iupred_file", nargs='*', help="Input IUPred/AIUpred result file(s)")
    parser.add_argument("-d", "--dir", help="Input dir")
    parser.add_argument("-ext", "--ext", default="iupred.txt.gz",help="file extension when reading folder")
    parser.add_argument("-o","--outfile", default='bigquery/idp.csv', 
                        help="Output file")
    parser.add_argument('-v','--debug', help='Debugging output', action='store_true')
    
    args = parser.parse_args()
    with open (args.outfile, "w",newline="") as outfh:
        outwriter = csv.writer(outfh)
        outwriter.writerow(["species_prefix","protein_id","IDP_start","IDP_end", "mean_score"])

        if args.iupred_file:
            for file in args.iupred_file:
                if args.debug:
                    print(f"Processing {file}")
                timestart = time.time() 
                iupred_regions = parse_iupred_file(file)
                timeend = time.time()
                print(f"Time elapsed: {timeend-timestart}")
                
                for seqid, idp_regions in iupred_regions.items():
                    prefix = seqid.split("_")[0]
                    for idp in idp_regions:
                        start, end, mean_score = idp
                        outwriter.writerow([prefix, seqid, start, end, f"{mean_score:.4f}"])
        else:
            print("not ready to do folder processing, usually this is slow")
if __name__ == "__main__":
    main()