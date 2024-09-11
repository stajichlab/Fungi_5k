#!/usr/bin/bash -l
#SBATCH -p gpu --gres=gpu:a100:1 -c 8 --mem 64gb -N 1 -n 1 --out logs/signalp.%a.log
module load signalp/6-gpu

CPU=2
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

WRITECPU=8
if [ $CPU -lt $WRITECPU ]; then
	WRITECPU=$CPU
fi

N=${SLURM_ARRAY_TASK_ID}
if [ -z $N ]; then
    N=$1
    if [ -z $N ]; then
        echo "need to provide a number by --array or cmdline"
        exit
    fi
fi
FILEBATCH=100 # how many files to process at a time

BATCH=100
INDIR=input
OUTDIR=results/function/signalp/
mkdir -p $OUTDIR
sampset=sampleset.txt
if [ ! -s $sampset ]; then
	ls -U $INDIR | sort > $sampset
fi
MAX=$(wc -l $sampset | awk '{print $1}')
START=$(perl -e "print 1 + (($N - 1) * $FILEBATCH)")
END=$(perl -e "print ($N * $FILEBATCH) - 1")
if [ $START -gt $MAX ]; then
	echo "$START too big for $MAX"
	exit
elif [ $END -gt $MAX ]; then
	END=$MAX
fi
echo "running $START - $END"
sed -n ${START},${END}p $sampset | while read INFILE
do
	NAME=$(basename $INFILE .proteins.fa)
	echo "$NAME"
	if [ ! -f $OUTDIR/${NAME}.signalp.results.txt.gz ]; then
		time signalp6 -od $SCRATCH/${NAME} -org euk --mode fast -format txt -fasta $INDIR/$INFILE --write_procs $WRITECPU -bs $BATCH
		pigz -c $SCRATCH/${NAME}/prediction_results.txt > $OUTDIR/${NAME}.signalp.results.txt.gz
		pigz -c $SCRATCH/${NAME}/processed_entries.fasta > $OUTDIR/${NAME}.signalp.processed_entries.fasta.gz
		pigz -c $SCRATCH/${NAME}/output.gff3 > $OUTDIR/${NAME}.signalp.gff3.gz

		rm -rf $SCRATCH/${NAME}
	fi
done
