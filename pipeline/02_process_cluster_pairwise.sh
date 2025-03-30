#!/usr/bin/bash -l
#SBATCH -p epyc --mem 384g -N 1 -n 1 -c 4 --out logs/mmseqs2pairwise.%A.log

module load workspace/scratch
module load biopython

filecount=$(ls -U input | wc -l | awk '{print $1}')
filecount=$(($filecount + 100))
if [[ $filecount -gt 1024 ]]; then
	ulimit -n $filecount
fi
pigz -dc results/Fungi5K_cluster.tsv.gz | time python scripts/mmseqs2pairwise.py --seqs input -o $SCRATCH/pairwise_MMseq_cluster -v
parallel -j 8 pigz -p 2 {} ::: $(ls -U $SCRATCH/pairwise_MMseq_cluster/*.tsv)
rsync -a $SCRATCH/pairwise_MMseq_cluster results

