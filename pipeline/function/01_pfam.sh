#!/usr/bin/bash -l
#SBATCH -p epyc -N 1 -n 1 -c 2 --mem 2gb --out logs/pfam.%a.log

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

INFILE=db/LsFMGC_AA_95_rep.fasta
TEMP=db/$(basename $INFILE .fasta)__split
PREFIX=LsFMGC
OUTDIR=results/function/pfam
mkdir -p $OUTDIR

IN=$TEMP/${PREFIX}.$N
rsync -a $PFAM_DB/Pfam-A.hmm* $SCRATCH/
time hmmscan --cut_ga --cpu $CPU \
    --domtblout $OUTDIR/${PREFIX}.$N.domtblout \
    --tblout $OUTDIR/${PREFIX}.$N.tblout \
    $SCRATCH/Pfam-A.hmm $IN | gzip -c > $OUTDIR/${PREFIX}.$N.log.gz
