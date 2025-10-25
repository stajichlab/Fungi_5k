#!/usr/bin/bash -l
#SBATCH -c 48 --mem 32gb -N 1 -n 1 --out logs/targetp.%a.log
module load targetp

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
FILEBATCH=4 # how many files to process at a time
INDIR=$(realpath input)
OUTDIR=results/function/targetP/
mkdir -p $OUTDIR
OUTDIR=$(realpath $OUTDIR)
sampset=sampleset.txt
if [ ! -s $sampset ]; then
	ls -U $INDIR | grep -v -P '\.fai$' | sort > $sampset
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

runtargetp() {
	INFILE=$1
	BATCH=50
	NAME=$(basename $INFILE .proteins.fa)
	echo "$NAME, $INDIR/$INFILE $OUTDIR/${NAME}"
	if [ ! -s $OUTDIR/${NAME}_summary.targetp2.gz ]; then
		time targetp -batch $BATCH -tmp $SCRATCH -format short -fasta $INDIR/$INFILE -org non-pl -prefix $OUTDIR/${NAME}
		pigz -f $OUTDIR/${NAME}_summary.targetp2
	fi

}
export -f runtargetp
export INDIR OUTDIR
parallel -j 4 runtargetp {} ::: $(sed -n ${START},${END}p $sampset)
