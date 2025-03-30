#!/usr/bin/bash -l
#SBATCH -p short -C ryzen --mem 24gb -c 8 --out logs/make_orthogroups.log

#pigz -dc results/Fungi5K_cluster.tsv.gz | ./scripts/mmseqs2bigqueryload.py -v

pigz -dc results/Fungi5K_cluster.tsv.gz | ./scripts/mmseqs2orthogroups.py -v
