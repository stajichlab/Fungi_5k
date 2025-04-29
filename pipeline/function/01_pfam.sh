#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 1 -c 48 --mem 48gb --out logs/pfam.%a.log -a 1-244

CPU=2
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

module load hmmer/3.4
module load db-pfam
module load workspace/scratch

FILEBATCH=24
INDIR=$(realpath input)
OUTDIR=results/function/pfam_rerun
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

runcmd() {
    INFILE=$1
    INCPU=2
    NAME=$(basename $INFILE .proteins.fa)
    echo "$NAME, $INDIR/$INFILE $OUTDIR/${NAME}"    
    if [ ! -s $OUTDIR/${NAME}.pfam ]; then
        time hmmscan --cut_ga --cpu $INCPU \
            --domtblout $OUTDIR/${NAME}.pfam \
            --tblout $OUTDIR/${NAME}.tblout \
            $SCRATCH/Pfam-A.hmm $INDIR/$INFILE | gzip -c > $OUTDIR/${NAME}.log.gz
    fi
}

export -f runcmd
export INDIR OUTDIR CPU PFAM_DB SCRATCH
RUNCPU=$(expr $CPU / 2)

rsync -a $PFAM_DB/Pfam-A.hmm* $SCRATCH/
parallel -j $RUNCPU runcmd {} ::: $(sed -n ${START},${END}p $sampset)
