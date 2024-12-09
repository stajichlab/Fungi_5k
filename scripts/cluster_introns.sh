#!/usr/bin/bash -l
#SBATCH -p short -c 96 -C xeon --mem 96gb --out logs/cluster_introns_all.log

module load mmseqs2
DB=fungiIntronDB
#mmseqs createdb all_fungi_introns.fa $DB --dbtype 2
#mmseqs cluster $DB $DB.clustered $SCRATCH
#mmseqs createtsv $DB $DB $DB.clustered ${DB}_clustered.tsv

mmseqs easy-clust --dbtype 2 all_fungi_introns.fa intron_cluster $SCRATCH
