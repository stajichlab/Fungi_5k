#!/usr/bin/env python3

import json
import os
import gzip
import csv
import argparse

def load_lineages(fh):
    lineages = {}
    for line in fh:
        
        row = line.strip().split("\t")
        taxid = row[0]
        if len(row) < 3:
            continue
        lineages[taxid] = {}
        names = row[1].split(";")
        ranknames = row[2].split(";")
        for i in range(0, len(names)):
            lineages[taxid][ranknames[i]] = names[i]
    return lineages

def load_samples(fh):
    samples = []
    reader = csv.DictReader(fh)
    for row in reader:
        samples.append(row)
    return samples

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Process Samples and look for missing data",
                                    epilog='Example: update_lineage_missing.py')
    parser.add_argument("--samples", default="samples.csv", help="samples.csv file for fungi5k")
    parser.add_argument("--lineage", default="samples.lineage.txt", 
                        help="samples.lineage.txt file by running ' cut -d, -f5 samples.csv | taxonkit lineage -R > samples.lineage.txt'")

    parser.add_argument('-v','--debug', help='Debugging output', action='store_true')

    parser.add_argument("-o","--outfile", default='samples.lineage_fixed.csv', 
                        help="Output file")
    
    args = parser.parse_args()
    with open(args.samples,"r") as fh,open(args.lineage,"r") as linfh, open(args.outfile, "w",newline="") as outfh:

        lineage = load_lineages(linfh)
        csvout = csv.writer(outfh)
        reader = csv.reader(fh)

        header = next(reader)
        header2col = {}
        for i,h in enumerate(header):
            header2col[h] = i     
        csvout.writerow(header)
        
        for row in reader:
            taxid = row[header2col['NCBI_TAXONID']]
            taxonomyrange = []
            #row[5:13]
            lastname = lineage[taxid]['kingdom']
            print(row[6:14])
            for name in ['phylum','subphylum','class','subclass','order','family','genus','species']:
                if name not in lineage[taxid]:
                    if args.debug:
                        print("Missing",taxid,name)
                    taxonomyrange.append("")
                else:
                    taxonomyrange.append(lineage[taxid][name])
                    lastname = lineage[taxid][name]
            row[6:14] = taxonomyrange
            csvout.writerow(row)
