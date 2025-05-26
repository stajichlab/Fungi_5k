#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 1 -c 48 --mem 84gb --out logs/iprscan.%a.log -a 1-1938

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

module load interproscan
module load funannotate
module load workspace/scratch

RUNCPU=1
FILEBATCH=3
INDIR=$(realpath input)
OUTDIR=results/function/interpro
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
    INCPU=$2
    NAME=$(basename $INFILE .proteins.fa)
    echo "$NAME, $INDIR/$INFILE $OUTDIR/${NAME}"    
    if [ ! -s $OUTDIR/${NAME}.iprscan.tsv.gz ]; then
        #time interproscan.sh -o $SCRATCH/${NAME}.iprscan.tsv --goterms  \
        #    --iprlookup --tempdir $SCRATCH --cpu $INCPU  \
        #    -f TSV -i $INDIR/$INFILE > $SCRATCH/${NAME}.iprscan.log
	#           
	IPRPATH=$(which interproscan.sh)
	pushd $SCRATCH
	time funannotate iprscan -o ${NAME}.iprscan.xml -i $INDIR/$INFILE -m local -c $INCPU --iprscan_path $IPRPATH > ${NAME}.iprscan.log
	pigz -c ${NAME}.iprscan.xml > $OUTDIR/${NAME}.iprscan.xml.gz
        mv *.log $OUTDIR/
	popd
	
	mkdir -p $SCRATCH/convert
	interproscan.sh -mode convert --tempdir $SCRATCH/convert -f tsv -i $SCRATCH/${NAME}.iprscan.xml -o $OUTDIR/${NAME}.iprscan.tsv
        pigz $OUTDIR/${NAME}.iprscan.tsv 
    fi
}

export -f runcmd
export INDIR OUTDIR CPU PFAM_DB SCRATCH
PERCPU=$(perl -e "print int($CPU / $RUNCPU)")
if [ $PERCPU -lt 1 ]; then
    echo "not enough CPU for $RUNCPU"
    exit
fi
if [ $PERCPU -gt 1 ]; then
    echo "using $PERCPU CPU for $RUNCPU"
fi
parallel -j $RUNCPU runcmd {} ::: $(sed -n ${START},${END}p $sampset) ::: $PERCPU
