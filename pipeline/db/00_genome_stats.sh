#!/usr/bin/bash -l
#SBATCH -p short -c 64 --mem 64gb --out logs/aaftf_assess.log

module load AAFTF
CPU=64
pushd genomes
parallel -j $CPU AAFTF assess -i {} -r {.}.stats.txt ::: $(ls -U *.fa)

./scripts/collect_asm_stats.py
