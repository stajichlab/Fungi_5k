#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 96 --mem 96gb -p short
#SBATCH --job-name=aiupred_sum
#SBATCH --output=logs/aiupred_summarize.log

CPU=1
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

EXT=idp.csv
EXTSUM=idp_summary.csv
if [ ! -s bigquery/$EXT.gz ]; then
    INDIR=results/function/aiupred
    
    ls -U $INDIR | grep \.aiupred\.txt\.gz | time parallel -J $CPU ./scripts/gather_AIUPred.py $INDIR/{} --outfile $SCRATCH/{.}.$EXT --outfilesum $SCRATCH/{.}.$EXTSUM
    FIRST=$(ls -U $SCRATCH/*.$EXT | head -n 1)
    head -n 1 $FIRST > bigquery/$EXT
    ls -U $SCRATCH | grep -E "$EXT" | xargs -I {} sh -c "tail -n +2 $SCRATCH/{} >> bigquery/$EXT"
    pigz -f bigquery/$EXT

    FIRST=$(ls -U $SCRATCH/*.$EXTSUM | head -n 1)
    head -n 1 $FIRST > bigquery/$EXTSUM
    ls -U $SCRATCH | grep -E "$EXTSUM" | xargs -I {} sh -c "tail -n +2 $SCRATCH/{} >> bigquery/$EXTSUM"
    pigz -f bigquery/$EXTSUM

fi