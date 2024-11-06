#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 96gb -p hpcc_default
#SBATCH --job-name=genegff_stats
#SBATCH --output=logs/genegff_seqstats.log

CPU=1
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi
module load biopython

EXT=gene_info.csv
if [ ! -s bigquery/$EXT ]; then
    INDIR=gff3
    # this should be fast enough
    time ./scripts/build_genestats_bigquery.py --gff_dir gff3
fi
