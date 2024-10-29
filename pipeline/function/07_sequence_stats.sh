#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 96 --mem 96gb -p short
#SBATCH --job-name=seq_stats
#SBATCH --output=logs/seqstats.log

CPU=1
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi
module load biopython

EXT=aa_freq.csv
if [ ! -s bigquery/$EXT ]; then
    INDIR=input
    
    ls -U $INDIR | grep \.fa | parallel -J $CPU ./scripts/calculate_AA_freq.py $INDIR/{} -v -o $SCRATCH/{.}.$EXT
    FIRST=$(ls -U $SCRATCH/*.$EXT | head -n 1)
    head -n 1 $FIRST > bigquery/$EXT
    ls -U $SCRATCH | grep -E "$EXT" | xargs -I {} sh -c "tail -n +2 $SCRATCH/{} >> bigquery/$EXT"
fi
EXT=codon_freq.csv
if [ ! -s bigquery/$EXT ]; then
    INDIR=input_cds
    ls -U $INDIR | grep \.fa | parallel -J $CPU ./scripts/calculate_codon_freq.py $INDIR/{} -v -o $SCRATCH/{.}.$EXT
    FIRST=$(ls -U $SCRATCH/*.$EXT | head -n 1)
    head -n 1 $FIRST > bigquery/c$EXT
    ls -U $SCRATCH | grep "$EXT" | xargs -I {} sh -c "tail -n +2 $SCRATCH/{} >> bigquery/$EXT.csv"
#    ./scripts/calculate_codon_freq.py -d input_cds -o bigquery/codon_freq.csv
fi