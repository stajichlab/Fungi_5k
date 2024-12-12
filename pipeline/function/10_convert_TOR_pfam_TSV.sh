#!/usr/bin/bash -l
#SBATCH -c 64 --mem 96gb --out logs/function_convert_pfam_bigquery.log
#SBATCH --job-name=convert_pfam

CPU=1
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

do_query() {
    local START=$1
    local NUM=$2
    local GFOLDER=$3
    local END=$(expr $START + $NUM - 1)
    OUTDIR=$SCRATCH/start_${START}
    echo "st=${START} end=${END}"
    mkdir -p $OUTDIR
    ./scripts/pfamtsv_to_long.py --outfile $OUTDIR/pfam.csv $(ls $GFOLDER | sed -n "${START},${END}p" | xargs -i echo "${GFOLDER}/"{})
}
export -f do_query
EXT=pfam.csv
OUTDIR=bigquery
BINSIZE=100
if [ ! -s $OUTDIR/$EXT.gz ]; then
    INDIR=results/function/pfam
    COUNT=$(ls $INDIR | wc -l)
    echo "starting -->"
    date
    time parallel -j $CPU do_query ::: $(seq 1 $BINSIZE $COUNT) ::: $BINSIZE ::: $INDIR
    for file in $SCRATCH/start_1/pfam.csv; do
        head -n 1 $file > $OUTDIR/$(basename $file)
    done
    for file in $SCRATCH/start_*/pfam.csv;
    do
        echo "processing ${file}"
        tail -n +2 $file >> ${OUTDIR}/$(basename $file)
    done
    date
    echo "--> Finished"
fi
pigz $EXT