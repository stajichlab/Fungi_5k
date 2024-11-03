#!/usr/bin/env python3

from Bio import SeqIO

import csv
import os
import sys
from contextlib import ExitStack
import argparse


def parse_gff(gff, dna="", debug=False):
    genedata = {}
    dnadb = None
    if dna:
        dnadb = SeqIO.index(dna)
    with open(gff, 'r') as gff_fh:
        transcript2gene = {}
        for line in gff_fh:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                if debug:
                    print(f"Skipping line with {len(fields)} fields in {gff}")
                continue            
            group_data = {}
            # how expensive is this
            for f in fields[8].split(';'):
                (tag,value) = f.split('=')
                group_data[tag] = value
            (fstart,fend) = sorted(int(fields[3]), int(fields[4]))
            fstrand = -1 if fields[6] == '-' else 1
            if fields[2] == 'gene':
                if 'ID' not in group_data:
                    print(f'Cannot parse groups {groups} does not start with ID in {gff}\n{line}')
                    continue
                gene_id = group_data['ID']
                if gene_id in genedata:
                    print(f"Duplicate gene ID {gene_id} in {gff}")
                    continue                
                genedata[gene_id] = {'chrom': fields[0], 
                                    'start': fstart, 'end': fend,
                                    'strand': fstrand,
                                    'type': 'NULL',
                                    'transcripts': {}, 
                                    }
            elif fields[2] == 'mRNA':
                if not ('ID' in group_data and 'Parent' in group_data):
                    print(f'Cannot parse groups {groups} does not have ID or Parent in {gff}\n{line}')
                    continue
                mrna_id = group_data['ID']
                gene_id = group_data['Parent']
                transcript2gene[mrna_id] = gene_id
                if gene_id not in genedata:
                    print(f"mRNA {mrna_id} has no gene in {gff}")
                    continue
                genedata[gene_id]['type'] = 'protein_coding'
                genedata[gene_id]['transcripts'][mrna_id] = {'chrom': fields[0], 
                                                    'start': fstart, 
                                                    'end': fend, 
                                                    'strand': fstrand, 
                                                    'is_partial': 'NULL',
                                                    'exons': [], 'CDS': [],
                                                    'introns': []}
            elif fields[2] == 'tRNA':
                
                if not ('ID' in group_data and 'Parent' in group_data):
                    print(f'Cannot parse groups {groups} does not have ID or Parent in {gff}\n{line}')
                    continue
                trna_id = group_data['ID']
                gene_id = group_data['Parent']
                transcript2gene[mrna_id] = gene_id
                if gene_id not in genedata:
                    print(f"tRNA {trna_id} has no gene in {gff}")
                    continue
                genedata[gene_id]['type'] = 'tRNA_gene'
                genedata[gene_id]['transcripts'][trna_id] = {'chrom': fields[0], 
                                                    'start': fstart, 
                                                    'end': fend, 
                                                    'strand': fstrand, 
                                                    'is_partial': 'FALSE',
                                                    'exons': [], 
                                                    'introns': []}
            elif fields[2] == 'exon' or fields[2] == "CDS":
                ftype = fields[2]               
                if not ('Parent' in group_data):
                    print(f'Cannot parse groups {groups} does not have at least Parent in {gff}\n{line}')
                    continue
                parent_id = group_data['Parent']
                gene_id = None
                if parent_id not in transcript2gene:
                    print(f"Exon from transcript {parent_id} cannot map to gene id in {gff}\n{line}")
                    continue
                else:
                    gene_id = transcript2gene[parent_id]
                
                if gene_id not in genedata or mrna_id not in genedata[gene_id]['transcripts']:
                    print(f"Exon of {mrna_id} has no gene or mRNA in {gff}\n{line}")
                    continue
                n = len(genedata[gene_id]['transcripts'][parent_id][ftype]) + 1
                if 'ID' in group_data:
                    exon_id = group_data['ID']
                else:
                    exon_id = f'{parent_id}.{ftype}{n}'
                
                genedata[gene_id]['transcripts'][mrna_id][ftype].append({
                    'chrom': fields[0], 
                    'start': fstart,
                    'end': fend, 
                    'strand': fstrand, 
                    'exon_number': None})
    for gene in genedata:
        for transcript in genedata[gene]['transcripts']:
            genedata[gene]['transcripts'][transcript]['exons'] = sorted(
                genedata[gene]['transcripts'][transcript]['exons'], 
                key=lambda x: x['strand'] * x['start'])
            if 'CDS' in genedata[gene]['transcripts'][transcript]:
                genedata[gene]['transcripts'][transcript]['CDS'] = sorted(
                    genedata[gene]['transcripts'][transcript]['CDS'], 
                    key=lambda x: x['strand'] * x['start'])
            

def main():
    if len(sys.argv) < 2:
        print("Usage: python gff_build_big_query.py gff or -d <dir>")
        sys.exit(1)
    
    parser = argparse.ArgumentParser(description="Calculate gene stats from GFF file(s)",
                                    epilog='Example: gff_build_big_query.py gff_file.gff -d dna_dir -p pep_dir')
    parser.add_argument("gff_file", nargs='*', help="Input GFF file(s)")
    parser.add_argument("-g", "--gff_dir", help="GFF files dir")
    parser.add_argument("-d", "--dna_dir", default="genomes", help="DNA files dir")
    parser.add_argument("-p", "--pep_dir", default="input", help="Protein files dir")
    parser.add_argument("-gffext", "--gffext", default="gff3", help="file extension when reading gff folder")
    parser.add_argument("-pepext", "--pepext", default="proteins.fa", help="file extension when reading proteins folder")    
    parser.add_argument("-dnaext", "--dnaext", default=".scaffolds.fa", help="file extension when reading proteins folder")
    
    
    parser.add_argument('-v','--debug', help='Debugging output', action='store_true')
    parser.add_argument("-o","--outdir", default='bigquery', 
                        help="Output folder for gene info, exons, introns, transcripts, proteins")
    
    args = parser.parse_args()
    print(args)
    if args.debug:
        print(f"Reading GFF files from {args.gff_dir} and DNA files from {args.dna_dir}")
    output_files = ['gene_info.csv', 'gene_exons.csv', 'gene_CDS.csv', 'gene_introns.csv', 'gene_transcripts.csv', 'gene_proteins.csv']
    with ExitStack() as stack:
        files = [stack.enter_context(open(f'{args.outdir}/{filename}', 'w', newline='')) for filename in output_files]
        genefile, exonfile, CDSfile, intronsfile, mrnafile, pepfile = files        
        genecsv = csv.writer(genefile)    
        exoncsv = csv.writer(exonfile)
        exoncsv = csv.writer(CDSfile)
        introncsv = csv.writer(intronsfile)
        mrnacsv = csv.writer(mrnafile)
        pepcsv = csv.writer(pepfile)
        genecsv.writerow(['gene_id', 'species_prefix', 'chrom', 'start', 'end', 'strand','gene_type'])
        mrnacsv.writerow(['gene_id', 'mrna_id', 'chrom', 'start', 'end', 'strand', 'is_partial'])
        exoncsv.writerow(['exon_id', 'mrna_id', 'exon_number', 'chrom', 'start', 'end', 'strand','GC_content'])
        introncsv.writerow(['intron_id', 'mrna_id', 'intron_number', 'chrom', 'start', 'end', 'strand',
                        'splice_5', 'splice_3', 'GC_content'])
        pepcsv.writerow(['protein_id', 'mrna_id', 'chrom', 'start', 'end', 'strand', 'length', 'md5checksum'])
        
        IDs_to_process = set()
        filenames = []
        if args.gff_dir and os.path.isdir(args.gff_dir):
            for gff_file in os.listdir(args.gff_dir):
                if gff_file.endswith(args.ext):
                    filename_stem = gff_file.replace(args.gffext,"")
                    filenames.append([filename_stem,os.path.join(args.gff_dir, gff_file)])
        elif args.gff_file:
            for gff_file in args.gff_file:
                filename_stem = gff_file.replace(args.gffext,"")
                filenames.append([filename_stem,args.gff_file])
        for gff_tuple in filenames:
            (fstem, gff_file) = gff_tuple
            dnafile = ""
            if args.dna_dir and os.path.isdir(args.dna_dir):
                dnafile = os.path.join(args.dna_dir, f'{fstem}.{args.dnaext}')
            genedata = parse_gff(gff=gff_file, dna=dnafile, debug=args.debug)
            if genedata:
                IDs_to_process.add(fstem)
            species = None
            for gene in genedata:
                if not species:
                    (species) = gene.split('_')[0]
                genecsv.writerow([gene, species, 
                                genedata[gene]['chrom'], 
                                genedata[gene]['start'], 
                                genedata[gene]['end'], 
                                genedata[gene]['strand'],
                                genedata[gene]['type']])
                # consider saving space by only encoding strand on the gene level
                for transcript in genedata[gene]['transcripts']:
                    m = genedata[gene]['transcripts'][transcript]
                    mrnacsv.writerow([gene, transcript, m['chrom'], m['start'], m['end'], 
                                    m['strand'], m['is_partial']])
                    for exon in m['exons']:
                        e = m['exons'][exon]
                        exoncsv.writerow([exon, transcript, e['exon_number'], 
                                        e['chrom'], 
                                        e['start'], 
                                        e['end'], 
                                        e['strand'],
                                        e['GC_content']])
                    for cds in m['CDS']:
                        c = m['CDS'][cds]
                        CDSfile.writerow([cds, transcript, c['exon_number'], 
                                        c['chrom'], c['start'], c['end'], c['strand'],
                                        c['GC_content']])
                    for intron in m['introns']:
                        i = m['introns'][intron]                                
                        introncsv.writerow([intron, transcript, i['intron_number'], 
                                        i['chrom'], i['start'], i['end'], i['strand'],                                        
                                        i['splice_5'],i['splice_3'], i['GC_content']])
if __name__ == "__main__":
    main()