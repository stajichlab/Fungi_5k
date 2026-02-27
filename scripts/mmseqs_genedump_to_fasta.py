#!/usr/bin/env python3

import csv
import argparse
import os

def genedump_to_fasta_cluster(input_file, output_dir_in,target_type="SUBPHYLUM",target_taxo_count=2):
    output_dir = f"{output_dir_in}.{target_type.lower()}{target_taxo_count}"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    clusters = {}
    taxo_count = {}
    with open(input_file, 'r') as infile:
        #transcript_id,peptide,orthogroup,group_count,pfam_null_ratio,SPECIES,PHYLUM,SUBPHYLUM,CLASS,ORDER,FAMILY
        reader = csv.DictReader(infile, delimiter=',')        
        for row in reader:           
            gene_id = row['transcript_id']
            sequence = row['peptide']
            orthogroup = row['orthogroup']
            group_count = row['group_count']
            nullratio = float(row['pfam_null_ratio'])
            null_ratio = f'{nullratio:.2f}'
            if not orthogroup in clusters:
                clusters[orthogroup] = []
                taxo_count[orthogroup] = set()
            target = row[target_type.upper()]
            taxo_count[orthogroup].add(target)
            strdesc = f"{gene_id} OG={orthogroup} CLUSTERSIZE={group_count} PFAMRATIO={null_ratio}"
            cols = []
            for other_col in ['PHYLUM','SUBPHYLUM','CLASS','ORDER','FAMILY','SPECIES']:
                cols.append(f"{other_col}={row[other_col]}")
            strdesc += " " + ";".join(cols)
            clusters[orthogroup].append(f">{strdesc}\n{sequence}\n")
    n = 0
    for group, sequences in clusters.items():
        if len(taxo_count[group]) < target_taxo_count:
            continue
        n += 1
        output_file = os.path.join(output_dir, f"{group}.fasta")
        with open(output_file, 'w') as outfile:
            outfile.writelines(sequences)
    print(f"Wrote {n} FASTA files to {output_dir}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert mmseqs genedump CSV to FASTA format.')
    parser.add_argument('-i','--input', dest="input_file",
                        help='Input mmseqs genedump CSV file', 
                        default='results/mmseqs_nopfam_genedump.csv')
    parser.add_argument('-o', '--output', dest='output_dir', 
                        help='Output dir for FASTA files',
                        default = 'results/mmseqs_nopfam_genedump')
    parser.add_argument('-t', '--target_type', dest='target_type', 
                        help='Target taxonomic level for filtering',
                        default = 'SUBPHYLUM')
    parser.add_argument('-c', '--target_taxo_count', dest='target_taxo_count', type=int,
                        help='Minimum number of target taxonomic groups required',
                        default = 2)
    
    args = parser.parse_args()
    genedump_to_fasta_cluster(args.input_file, args.output_dir, args.target_type, args.target_taxo_count)