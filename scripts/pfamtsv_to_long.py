#!/usr/bin/env python3

import csv
import argparse
import os
import gzip

def process_file(fh, csvout, multicols):
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        i = 0
        for domain in row['query_name'].split(','):
            outrow = [ row['target_name'], domain ]
            for col in multicols:
                outrow.append(row[col].split(',')[i])
            csvout.writerow(outrow)
            i += 1

def main():    
    parser = argparse.ArgumentParser(description="Convert Pfam TSV from Toronto to long form one line per domain",
                                    epilog='Example: pfamtsv_to_long.py')
    parser.add_argument("-i","--indir", help="Directory with pfam results pre-computed TSV files")
    parser.add_argument("tsv", nargs="*", help="Input TSV file(s)")
    parser.add_argument("-ext","--extension", default="tsv", help="extension of pfam file [tsv]")
    parser.add_argument("-o","--outfile", default="bigquery/pfam.csv", 
                        help="output file [bigquery/pfam.csv]")
    parser.add_argument('-v','--debug', help='Debugging output', action='store_true')
    
    args = parser.parse_args()
    
    with open(args.outfile, "w",newline="") as outfh:
        
        csvout = csv.writer(outfh)
        multicols = ['full_seq_e_value', 'full_seq_score', 'full_seq_bias',
                    'domain_num', 'domain_num_of', 'domain_c_evalue', 'domain_i_evalue',
                    'domain_score', 'domain_bias', 'hmm_from', 'hmm_to', 'ali_from', 'ali_to',
                    'env_from','env_to']
        header = ['protein_id', 'pfam_id']
        header.extend(multicols)
        csvout.writerow(header)
        n = 0
        if args.indir:
            for file in os.listdir(args.indir):
                #if args.debug:
                #    print(f"Processing {file}")
                fh = None
                if file.endswith("." + args.extension):
                    filepath = os.path.join(args.indir, file)
                    fh = open(filepath, "r")
                elif file.endswith("." + args.extension + ".gz"):
                    filepath = os.path.join(args.indir, file)
                    fh = gzip.open(filepath, "rt")
                else:
                    continue
                process_file(fh, csvout, multicols)
                n += 1
                if args.debug and n > 1 and n % 1000 == 0:
                    print(f'Processed {n} files')
                fh.close()
        elif args.tsv:
            for file in args.tsv:
                fh = None
                if file.endswith("." + args.extension):
                    fh = open(file, "r")
                elif file.endswith("." + args.extension + ".gz"):
                    fh = gzip.open(file, "rt")
                else:
                    continue
                process_file(fh, csvout, multicols)
                fh.close()


if __name__ == "__main__":
    main()