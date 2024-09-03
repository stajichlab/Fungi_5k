#!/usr/bin/bash -l
#SBATCH -p epyc -N 1 -n 48 --mem 384gb  --out logs/mmseqs.%A.log

# Load modules
module load mmseqs2
module load workspace/scratch
PREFIX=Fungi5K
mkdir -p $SCRATCH/mmseqs
DB=$SCRATCH/mmseqs/$PREFIX
mkdir -p results
time mmseqs createdb input/*.proteins.fa $DB --compressed 1 --dbtype 2
mmseqs cluster $DB ${DB}_cluster $SCRATCH --compressed 1 --cluster-reassign 1
mmseqs createtsv $DB $DB ${DB}_cluster results/${PREFIX}_cluster.tsv.gz --compressed 1
