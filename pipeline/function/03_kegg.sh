#!/usr/bin/bash -l
#SBATCH -p epyc -N 1 -n 1 -c 8 --mem 8gb 
#SBATCH --output=logs/domains.kofamscan.%a.log
#SBATCH --array=1-506

module load kofamscan

CPUS=$SLURM_CPUS_ON_NODE
if [ ! $CPUS ]; then
    CPUS=1
fi

IN=${SLURM_ARRAY_TASK_ID}

if [ ! $IN ]; then
    IN=$1
    if [ ! $IN ]; then
        IN=1
        echo "defaulting to IN value is 1 - specify with --array or cmdline"
    fi
fi

ORIG=db/LsFMGC_AA_95_rep.fasta
TEMP=$(realpath db/$(basename $ORIG .fasta)__split)
PREFIX=LsFMGC
OUTDIR=results/function/kofam
mkdir -p $OUTDIR
OUTDIR=$(realpath $OUTDIR)

INFILE=$TEMP/${PREFIX}.$IN
OUT=$OUTDIR/${PREFIX}.$IN.kofam.tsv
OUTRICH=$OUTDIR/${PREFIX}.$IN.kofam.txt
KOLIST=$KOFAM_DB/ko_list
PROFILES=$KOFAM_DB/profiles/prokaryote.hal
# will update to specify profile folder and ko_list file
if [ ! -s $OUT ]; then
	pushd $SCRATCH
	if [ ! -f $OUT ]; then
		time exec_annotation -o $OUT --cpu $CPUS -f mapper -E 0.0001 -k $KOLIST --profile=$PROFILES $INFILE
	fi
	if [ ! -f $OUTRICH.gz ]; then
		time exec_annotation -o $OUTRICH --cpu $CPUS -f detail -E 0.0001 -k $KOLIST --profile=$PROFILES $INFILE
		pigz $OUTRICH
	fi
fi