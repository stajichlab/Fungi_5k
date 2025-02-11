#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 64 --mem 96gb -p hpcc_default
#SBATCH --job-name=genedist_stats
#SBATCH --output=logs/genedist_seqstats.log

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
    ./scripts/calculate_intergenic.py --outdir $OUTDIR $(ls $GFOLDER | sed -n "${START},${END}p" | xargs -i echo "${GFOLDER}/"{})
}
export -f do_query
EXT=gene_pairwise_distances.csv
OUTDIR=bigquery
BINSIZE=100
if [ ! -s $OUTDIR/$EXT ]; then
    INDIR=gff3
    COUNT=$(ls $INDIR | wc -l)
    echo "starting -->"
    date
    time parallel -j $CPU do_query ::: $(seq 1 $BINSIZE $COUNT) ::: $BINSIZE ::: $INDIR
    for file in $SCRATCH/start_1/*.csv; do
        head -n 1 $file > $OUTDIR/$(basename $file)
    done
    for file in $SCRATCH/start_*/*.csv;
    do
        echo "processing ${file}"
        tail -n +2 $file >> ${OUTDIR}/$(basename $file)
    done
    date
    echo "--> Finished"
    # this should be fast enough (10 hrs for 5800 files)
    #time ./scripts/calculate_intergenic.py --gff_dir gff3 
fi
