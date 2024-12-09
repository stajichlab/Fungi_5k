#!/usr/bin/env python3

import json
import os
import gzip
import csv
import argparse


def load_samples(file_path):
    samples = []
    with open(file_path, 'r') as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            samples.append(row)
    return samples

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Process Samples and look for missing data",
                                    epilog='Example: fix_samples.py')
    parser.add_argument("--samples", default="samples.csv", help="samples.csv file for fungi5k")
    parser.add_argument('-v','--debug', help='Debugging output', action='store_true')

    parser.add_argument("-o","--outfile", default='samples.fixed.csv', 
                        help="Output file")
    
    args = parser.parse_args()
    with open(args.samples,"r") as fh, open(args.outfile, "w",newline="") as outfh:
        csvout = csv.writer(outfh)
        reader = csv.reader(fh)

        header = next(reader)
        header2col = {}
        for i,h in enumerate(header):
            header2col[h] = i     
        csvout.writerow(header)
        
        for row in reader:
            if not row[header2col['GENUS']]:
                row[header2col['GENUS']] = row[header2col['SPECIES']].split(" ")[0]
            csvout.writerow(row)
