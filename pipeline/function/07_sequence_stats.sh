#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 8 --mem 96gb -p short
#SBATCH --job-name=seq_stats
#SBATCH --output=logs/seqstats.log

module load biopython

./scripts/calculate_AA_freq.py -d input