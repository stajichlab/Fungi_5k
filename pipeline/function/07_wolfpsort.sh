#!/usr/bin/bash -l
#SBATCH -p short -c 50 --mem 24gb -N 1 -n 1 --out logs/wolfpsort.%a.log
module load wolfpsort

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
OUTDIR=results/function/wolfpsort
mkdir -p $OUTDIR
OUTDIR=$(realpath $OUTDIR)
sampset=sampleset.txt
if [ ! -s $sampset ]; then
	ls -U $INDIR | sort > $sampset
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

runwolfsort() {
	INFILE=$1
	NAME=$(basename $INFILE .proteins.fa)
	echo "$NAME"
	if [ ! -f $OUTDIR/${NAME}.wolfpsort.results.txt.gz ]; then
		cat $INDIR/$INFILE | time runWolfPsortSummary fungi > $OUTDIR/${NAME}.wolfpsort.results.txt
		pigz  $OUTDIR/${NAME}.wolfpsort.results.txt
	fi

}
export -f runwolfsort
export INDIR OUTDIR
parallel -j $CPU runwolfsort {} ::: $(sed -n ${START},${END}p $sampset)
