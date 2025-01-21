#!/usr/bin/bash -l

#SBATCH --mem 24gb --out split.log

module load bioperl

mkdir -p split
pushd split

bp_dbsplit --size 500000 --prefix fungi_introns  ../all_fungi_introns.fa
