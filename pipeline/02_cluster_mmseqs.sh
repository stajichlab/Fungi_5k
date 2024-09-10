#!/usr/bin/bash -l
#SBATCH -p epyc -N 1 -n 96 --mem 384gb  --out logs/mmseqs.%A.log

# Load modules
module load mmseqs2
which mmseqs
module load workspace/scratch
PREFIX=Fungi5K
mkdir -p $SCRATCH/mmseqs
DB=$SCRATCH/mmseqs/$PREFIX
#mkdir -p mmseqs
#DB=mmseqs/$PREFIX
mkdir -p results
time mmseqs createdb input/*.proteins.fa $DB --compressed 1 --dbtype 1
# use shorter sequence to compute alignment length but don't limit it to just the alignment part
time mmseqs cluster $DB ${DB}_cluster $SCRATCH --compressed 1 --cluster-reassign 1 --min-seq-id 0.30 --seq-id-mode 1 -c 0.70
mmseqs createtsv $DB $DB ${DB}_cluster results/${PREFIX}_cluster.tsv
pigz results/${PREFIX}_linclust.tsv
