#!/usr/bin/bash -l
#SBATCH -c 24 -N 1 -n 1 --mem 24gb --out logs/fungi_intron.blast.%a.log
CPUS=$SLURM_CPUS_ON_NODE
if [ ! $CPUS ]; then
    CPUS=1
fi


IN=${SLURM_ARRAY_TASK_ID}

if [ -z $IN ]; then
    IN=$1
    if [ -z $IN ]; then
        IN=1
        echo "defaulting to IN value is 1 - specify with --array or cmdline"
    fi
fi

module load ncbi-blast
DB=all_fungi_introns.fa
PID=80
PREF=fungi_introns
OUT=blastn_results
mkdir -p $OUT
if [ ! -f $OUT/all_fungi.self-vs-self.${IN}.BLASTN.tab.gz ]; then
	time blastn -task blastn -perc_identity $PID -num_threads $CPUS -num_alignments 1000000 -reward 1 -penalty -1 -gapopen 2 -gapextend 1 \
	-db $DB -query split/$PREF.${IN} -out $OUT/all_fungi.self-vs-self.${IN}.BLASTN.tab -outfmt "6 qaccver qlen saccver slen pident length mismatch gapopen qstart qend sstart send evalue bitscore"

	pigz $OUT/all_fungi.self-vs-self.${IN}.BLASTN.tab
fi
