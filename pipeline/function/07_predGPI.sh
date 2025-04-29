#!/usr/bin/bash -l
#SBATCH -p short -c 96 --mem 96gb -N 1 -n 1 --out logs/predgpi.%a.log
module load predgpi

CPU=1
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

N=${SLURM_ARRAY_TASK_ID}
if [ -z $N ]; then
    N=$1
    if [ -z $N ]; then
        echo "need to provide a number by --array or cmdline"
        exit
    fi
fi
FILEBATCH=500 # how many files to process at a time
INDIR=$(realpath input)
OUTDIR=results/function/predgpi
mkdir -p $OUTDIR
OUTDIR=$(realpath $OUTDIR)
sampset=sampleset.txt
if [ ! -s $sampset ]; then
	ls -U $INDIR | grep -v -P "\.fai$" | sort > $sampset
fi
sampset=$(realpath sampleset.txt)
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

runpredgpi() {
	INFILE=$1
	NAME=$(basename $INFILE .proteins.fa)
	echo "$NAME"
	if [ ! -f $OUTDIR/${NAME}.predgpi.gff3.gz ]; then
		predgpi.py -f $INDIR/$INFILE -m gff3 -o $OUTDIR/${NAME}.predgpi.gff3
		pigz $OUTDIR/${NAME}.predgpi.gff3
	fi

}
export -f runpredgpi
export INDIR OUTDIR
parallel -j $CPU runpredgpi {} ::: $(sed -n ${START},${END}p $sampset)
