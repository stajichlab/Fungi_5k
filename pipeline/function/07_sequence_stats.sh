#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 8 --mem 96gb -p short
#SBATCH --job-name=seq_stats
#SBATCH --output=logs/seqstats.log

module load biopython

if [ ! -s bigquery/aa_freq.csv ]; then
    ./scripts/calculate_AA_freq.py -d input -o bigquery/aa_freq.csv
fi
if [ ! -s bigquery/codon_freq.csv ]; then
    ./scripts/calculate_codon_freq.py -d input_cds -o bigquery/codon_freq.csv
fi