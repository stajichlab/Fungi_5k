#!/usr/bin/bash -l
#SBATCH -p epyc --mem 96g -N 1 -n 1 -c 2 --out logs/mmseqs2pairwise_linclust.%A.log

module load workspace/scratch
module load biopython

filecount=$(ls -U input | wc -l | awk '{print $1}')
filecount=$(($filecount + 100))
if [[ $filecount -gt 1024 ]]; then
	ulimit -n $filecount
fi
echo "There are $filecount files"
pigz -dc results/Fungi5K_linclust.tsv.gz | time python scripts/mmseqs2pairwise.py --seqs input -o results/pairwise_MMseq_linclust -v


