#!/usr/bin/bash -l
#SBATCH -c 24 --mem 96gb --out logs/filter_aln.%A.log
module load phyling
CPU=${SLURM_CPUS_ON_NODE}
if [ -z $CPU ]; then
    CPU=1
fi

phyling filter -I fungi_phyling_align -t $CPU -o fungi_msa_filter --verbose -n 50
