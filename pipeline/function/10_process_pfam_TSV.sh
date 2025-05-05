#!/usr/bin/bash -l
#SBATCH -p short -c 48 --mem 96gb --out logs/function_convert_pfam_bigquery.log
#SBATCH --job-name=convert_pfam


# convert the Toronto Pfam results into SQL table loadable csv

CPU=1
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

do_query() {
    local START=$1
    local NUM=$2
    local GFOLDER=$3
    local FILELIST=$4
    echo "do_query, $FILELIST $START $NUM $GFOLDER"
    local END=$(expr $START + $NUM - 1)
    OUTDIR=$SCRATCH/start_${START}
    echo "st=${START} end=${END}"
    mkdir -p $OUTDIR
    ./scripts/pfamtbl_to_long.py --outfile $OUTDIR/pfam.csv $(sed -n "${START},${END}p" $FILELIST | xargs -i echo "${GFOLDER}/"{})
}
export -f do_query
EXT=pfam.csv
FILEEXT=domtblout
OUTDIR=bigquery
BINSIZE=100
if [ ! -s $OUTDIR/$EXT.gz ]; then
    INDIR=results/function/pfam
    TEMPFILE=$SCRATCH/filelist
    ls -U $INDIR | grep -P "\.$FILEEXT\$" > $TEMPFILE
    COUNT=$(wc -l $SCRATCH/filelist | awk '{print $1}')
    echo "count is $COUNT"
    echo "starting -->"
    date
    time parallel -j $CPU do_query ::: $(seq 1 $BINSIZE $COUNT) ::: $BINSIZE ::: $INDIR ::: $TEMPFILE
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
if [ -f $OUTDIR/$EXT ]; then
    echo "compressing $OUTDIR/$EXT"
    pigz $OUTDIR/$EXT
fi