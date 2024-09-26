#!/usr/bin/bash -l
#SBATCH -p highmem -c 32 -N 1 -n 1 --mem 900gb --out logs/phyling_align.log

module load phyling
CPU=2
if  [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi
phyling align -I input -m fungi_odb10 -t $CPU -o phyling_align -v

