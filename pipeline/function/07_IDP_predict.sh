#!/usr/bin/bash -l
#SBATCH -p gpu --gres=gpu:a100:2 -c 4 --mem 100gb -N 1 -n 1 --out logs/IDP_predict.%a.log -a 1-6
module load aiupred

CPU=1
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

# this is really a limitation of how much memory the jobs take
# and the GPU node, seems like 2 max for p100 nodes


N=${SLURM_ARRAY_TASK_ID}
if [ -z $N ]; then
    N=$1
    if [ -z $N ]; then
        echo "need to provide a number by --array or cmdline"
        exit
    fi
fi
FILEBATCH=1000 # how many files to process at a time
INDIR=$(realpath input)
OUTDIR=results/function/aiupred
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

runaiupred() {
	INFILE=$1
	CUDA=$2
	NAME=$(basename $INFILE .proteins.fa)
	echo "$NAME"
	if [ ! -f $OUTDIR/${NAME}.aiupred.txt.gz ]; then
		if [[ $CUDA -lt 0 ]]; then
			echo "running on CPU"
			time aiupred.py -i $INDIR/$INFILE -o $OUTDIR/${NAME}.aiupred.txt
		else
			echo "running on GPU"
			time aiupred.py -g $CUDA -i $INDIR/$INFILE -o $OUTDIR/${NAME}.aiupred.txt
		fi
		pigz $OUTDIR/${NAME}.aiupred.txt
	fi
}
export -f runaiupred 
if [ -z $CUDA_VISIBLE_DEVICES ]; then
	CUDA_VISIBLE_DEVICES=-1
	PARALLELJOBS=$(expr $CPU / 48 + 1)
	echo "running CPU jobs with $PARALLELJOBS"
else
	PARALLELJOBS=$(echo $CUDA_VISIBLE_DEVICES | grep -c "," | awk '{print 1+$1}')
	echo "running GPU jobs with $PARALLELJOBS"

fi
export INDIR OUTDIR 
echo "running $PARALLELJOBS with CUDA=$CUDA_VISIBLE_DEVICES"
parallel --link -j $PARALLELJOBS runaiupred {} ::: $(sed -n ${START},${END}p $sampset) ::: $(echo $CUDA_VISIBLE_DEVICES | perl -p -e 's/,/ /')

