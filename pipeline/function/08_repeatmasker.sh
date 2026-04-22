#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 1 -c 16 --mem 64gb --out logs/repeatmasker.%a.log -a 1-727

CPU=16
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

module load RepeatMasker/4.1.8

FILEBATCH=8
INDIR=$(realpath genomes)
OUTDIR=results/function/repeatmasker
mkdir -p $OUTDIR
OUTDIR=$(realpath $OUTDIR)

sampset=genomeset.txt
if [ ! -s $sampset ]; then
    ls -U $INDIR | grep '\.scaffolds\.fa$' | sort > $sampset
fi
sampset=$(realpath $sampset)

MAX=$(wc -l $sampset | awk '{print $1}')
START=$(perl -e "print 1 + (($N - 1) * $FILEBATCH)")
END=$(perl -e "print ($N * $FILEBATCH)")
if [ $START -gt $MAX ]; then
    echo "$START too big for $MAX"
    exit
elif [ $END -gt $MAX ]; then
    END=$MAX
fi
echo "running $START - $END"

runcmd() {
    INFILE=$1
    INCPU=4
    NAME=$(basename $INFILE .scaffolds.fa)
    SAMPLE_OUTDIR=$OUTDIR/${NAME}
    mkdir -p $SAMPLE_OUTDIR
    MASKED=$SAMPLE_OUTDIR/${NAME}.scaffolds.fa.masked
    if [ ! -s $MASKED ] && [ ! -s ${MASKED}.gz ]; then
        echo "Running RepeatMasker on $NAME"
        RepeatMasker -pa $INCPU -species fungi -gff -xsmall \
            -dir $SAMPLE_OUTDIR $INDIR/$INFILE
        if [ -f $MASKED ]; then
            pigz $SAMPLE_OUTDIR/${NAME}.scaffolds.fa.out
            pigz $MASKED
        fi
    else
        echo "Skipping $NAME (already done)"
    fi
}

export -f runcmd
export INDIR OUTDIR

RUNCPU=$(expr $CPU / 4)
parallel -j $RUNCPU runcmd {} ::: $(sed -n ${START},${END}p $sampset)
