#!/usr/bin/bash -l
#SBATCH -p epyc -N 1 -n 64 --mem 384gb  --out logs/mmseq_linclust.%A.log

# Load modules
module load mmseqs2
module load workspace/scratch
PREFIX=Fungi5K
mkdir -p $SCRATCH/mmseqs
DB=$SCRATCH/mmseqs/$PREFIX
mkdir -p results
time mmseqs createdb input/*.proteins.fa $DB --compressed 1 --dbtype 1
#mmseqs cluster $DB ${DB}_cluster $SCRATCH --compressed 1 --cluster-reassign 1
time mmseqs linclust $DB ${DB}_lincluster $SCRATCH --compressed 1 

mmseqs createtsv $DB $DB ${DB}_lincluster results/${PREFIX}_linclust.tsv
pigz results/${PREFIX}_linclust.tsv
