#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 96 --mem 96gb -p short
#SBATCH --job-name=seq_stats
#SBATCH --output=logs/seqstats.log

CPU=1
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi
module load biopython

EXT=aa_freq.csv
if [ ! -s bigquery/$EXT.gz ]; then
    INDIR=input
    
    ls -U $INDIR | grep \.fa | parallel -J $CPU ./scripts/calculate_AA_freq.py $INDIR/{} -v -o $SCRATCH/{.}.$EXT
    FIRST=$(ls -U $SCRATCH/*.$EXT | head -n 1)
    head -n 1 $FIRST > bigquery/$EXT
    ls -U $SCRATCH | grep -E "$EXT" | xargs -I {} sh -c "tail -n +2 $SCRATCH/{} >> bigquery/$EXT"
    pigz -f bigquery/$EXT
fi
EXT=codon_freq.csv
if [ ! -s bigquery/$EXT.gz ]; then
    INDIR=input_cds
    ls -U $INDIR | grep \.fa | parallel -J $CPU ./scripts/calculate_codon_freq.py $INDIR/{} -v -o $SCRATCH/{.}.$EXT
    FIRST=$(ls -U $SCRATCH/*.$EXT | head -n 1)
    head -n 1 $FIRST > bigquery/$EXT
    ls -U $SCRATCH | grep "$EXT" | xargs -I {} sh -c "tail -n +2 $SCRATCH/{} >> bigquery/$EXT.csv"
    pigz -f bigquery/$EXT
fi
EXT=gene_info.csv
if [ ! -s bigquery/$EXT.gz ]; then
    INDIR=gff3
    # do these in parallel
    time ls -U $INDIR | grep \.gff3 | parallel -j $CPU ./scripts/build_genestats_bigquery.py $INDIR/{} --outdir $SCRATCH/{.}.${EXT}    
    FIRST=$(ls -Ud $SCRATCH/*.${EXT} | head -n 1)
    for FILE in $(ls $FIRST/*.csv); do
        FNAME=$(basename $FILE)
        head -n 1 $FILE > bigquery/$FNAME
        ls -U $SCRATCH/ | grep "$EXT" | xargs -I {} sh -c "tail -n +2 $SCRATCH/{}/${FNAME} >> bigquery/${FNAME}"
        pigz -f bigquery/${FNAME}
    done

fi
