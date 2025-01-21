#!/usr/bin/bash -l
#SBATCH --mem 96gb -N 1 -n 1 -c 24 --out usearch12.log

module load usearch/12
CPU=24
DB=all_fungi_introns.udb
IN=all_fungi_introns.fa
usearch  -usearch_global $IN -db $DB -id 0.8 -strand both -maxaccepts 20 -maxrejects 256 \
	-threads $CPU -output all_fungi.usearch.uc
