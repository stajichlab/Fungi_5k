#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 1 -c 16 --mem 8gb --out logs/cazy.%a.log

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

module load dbcanlight
module load workspace/scratch

INDIR=input
OUTDIR=results/function/cazy/
mkdir -p $OUTDIR

INFILE=$(ls -U $INDIR | sed -n ${N}p)
NAME=$(basename $INFILE .proteins.fa)
echo "To Run CAZY on Sample $NAME ($N): $INFILE"
mkdir -p $OUTDIR/$NAME
if [ ! -f $OUTDIR/$NAME/overview.tsv.gz ]; then
	time dbcanlight search -i $INDIR/$INFILE -m cazyme -o $OUTDIR/$NAME -t $CPU
	time dbcanlight search -i $INDIR/$INFILE -m sub -o $OUTDIR/$NAME -t $CPU
	#time dbcanlight search -i $INDIR/$INFILE -m diamond -o $OUTDIR/$NAME -t $CPU 
	dbcanlight conclude $OUTDIR/$NAME
    pigz -f $OUTDIR/$NAME/cazymes.tsv
    pigz -f $OUTDIR/$NAME/substrates.tsv
    pigz -f $OUTDIR/$NAME/overview.tsv
fi

