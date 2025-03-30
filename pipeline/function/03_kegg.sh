#!/usr/bin/bash -l
#SBATCH -p epyc -N 1 -n 1 -c 8 --mem 8gb 
#SBATCH --output=logs/domains.kofamscan.%a.log
#SBATCH --array=1-1215

module load kofamscan

CPUS=$SLURM_CPUS_ON_NODE
if [ ! $CPUS ]; then
    CPUS=1
fi

N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        N=1
        echo "defaulting to IN value is 1 - specify with --array or cmdline"
    fi
fi
INDIR=input
OUTDIR=$(realpath results/function/kofam)
mkdir -p $OUTDIR
sampset=sampleset.txt
if [ ! -s $sampset ]; then
	ls -U $INDIR | sort > $sampset
fi
INFILE=$(sed -n ${N}p $sampset)
PREFIX=$(basename $INFILE .fasta)
INFILE=$(realpath $INDIR/$INFILE)

mkdir -p $OUTDIR/$NAME

OUT=$OUTDIR/${PREFIX}.kofam.tsv
OUTRICH=$OUTDIR/${PREFIX}.kofam.txt
KOLIST=$KOFAM_DB/ko_list
PROFILES=$KOFAM_DB/profiles/eukaryote.hal
# will update to specify profile folder and ko_list file
if [ ! -s $OUT.gz ]; then
	pushd $SCRATCH
	if [ ! -f $OUT ]; then
		time exec_annotation -o $OUT --cpu $CPUS -f mapper -E 0.0001 -k $KOLIST --profile=$PROFILES $INFILE
	fi
	pigz $OUT
	if [ ! -f $OUTRICH.gz ]; then
		time exec_annotation -o $OUTRICH --cpu $CPUS -f detail -E 0.0001 -k $KOLIST --profile=$PROFILES $INFILE
		pigz $OUTRICH
	fi
fi