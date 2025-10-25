#!/usr/bin/bash -l
#SBATCH -p short -C ryzen --mem 24gb -c 8 --out logs/make_orthogroups.log
PREFIX=CryoendoAnt
pigz -dc results/${PREFIX}_cluster.tsv.gz | ./scripts/mmseqs2bigqueryload.py -v
pigz -dc results/${PREFIX}_cluster.tsv.gz | ./scripts/mmseqs2orthogroups.py -v
