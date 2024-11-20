#!/usr/bin/env python3

import json
import os
import gzip
import csv
import argparse
import re
import configparser

def load_samples(fh):
    samples = []
    reader = csv.DictReader(fh)
    for row in reader:
        samples.append(row)
    return samples

def main():    
    parser = argparse.ArgumentParser(description="Collect precomputed genome asm stats into a table",
                                    epilog='Example: collect_asm_stats.py')
    parser.add_argument("-d","--genomedir", default="genomes", help="Directory with genomes with pre-computed .stats.txt file")
    parser.add_argument("--samples", default="samples.csv", help="samples.csv file for fungi5k")
    parser.add_argument("-o","--outfile", default="bigquery/asm_stats.csv", 
                        help="output file [bigquery/asm_stats.csv]")
    parser.add_argument('-v','--debug', help='Debugging output', action='store_true')
    
    args = parser.parse_args()
    
    with open(args.samples,"r") as fh, open(args.outfile, "w",newline="") as outfh:
        species = load_samples(fh)
        csvout = csv.writer(outfh)

        # header is 
        statheader = [
            'CONTIG_COUNT', 'TOTAL_LENGTH', 'MIN', 'MAX', 'MEDIAN', 'MEAN',
            'L50', 'N50', 'L90', 'N90', 'GC%', 'T2T_SCAFFOLDS', 
            'TELOMERE_FWD', 'TELOMERE_REV'
        ]
        header = [ 'LOCUSTAG', 'SPECIES', 'FILENAME']
        header.extend(statheader)
        csvout.writerow(header)
        
        for sp in species:
            species_string = f"{sp['SPECIES']}_{sp['STRAIN']}".replace(' ', '_')
            stemname = f"{species_string}.scaffolds.stats.txt"
            if args.debug:
                print(stemname)
            row= [ sp['LOCUSTAG'], sp['SPECIES'], species_string ]
            
            statsfile = os.path.join(args.genomedir,stemname)
            if not os.path.exists(statsfile):
                if args.debug:
                    print("Missing",statsfile)
                continue
            statistics = {}
            
            with open(statsfile,"r") as statsfh:
                for line in statsfh:
                    line = line.strip()
                    if line.startswith("#"):
                        continue
                    if line.startswith("Assembly statistics"):
                        # this is just a sanity check to make sure the filename based
                        # used to open also matches the name that this asm was run on
                        line = re.sub(r'^Assembly statistics for:\s+','',line)
                        asmfile = line.replace('.scaffolds.fa','')
                        if asmfile != species_string:
                            print(f"reading {statsfile} for {species_string} but found {asmfile} in file")
                        continue
                    line = re.sub(r'^\s+','',line)          # remove leading whitespace
                    line = re.sub(r'\s+$','',line)          # remove trailing whitespace                    
                    (name,value) = line.split("=")          # split on the first '='
                    name = re.sub('\s+$','',name)           # remove trailing whitespace
                    value = re.sub('^\s+','',value)           # remove trailing whitespace
                    name = name.replace(' ','_')
                    statistics[name] = value

            for key in statheader:
                if key in statistics:
                    row.append(statistics[key])
                else:
                    print(f"missing statistic {key} in {statsfile}")
                    
            csvout.writerow(row)


if __name__ == "__main__":
    main()