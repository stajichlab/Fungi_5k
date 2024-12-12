#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 8 --mem 2gb --time 2:00:00
#SBATCH --job-name=MEROPS.domains
#SBATCH --output=logs/merops.%a.log

module load db-merops/124
module load ncbi-blast/2.16.0+

CPU=2
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

WRITECPU=8
if [ $CPU -lt $WRITECPU ]; then
	WRITECPU=$CPU
fi

N=${SLURM_ARRAY_TASK_ID}
if [ -z $N ]; then
    N=$1
    if [ -z $N ]; then
        echo "need to provide a number by --array or cmdline"
        exit
    fi
fi

FILEBATCH=100 # how many genome files to process at a time
INDIR=input
OUTEXT=blasttab
OUTDIR=results/function/merops
MEROPS_CUTOFF=1e-10
MEROPS_MAX_TARGETS=10
mkdir -p $OUTDIR
sampset=sampleset.txt

if [ ! -s $sampset ]; then
	ls -U $INDIR | sort > $sampset
fi
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

if [ ! $MEROPS_DB ]; then
    echo "Need a MEROPS_DB env variable either from config.txt or 'module load db-merops'"
    exit
fi

sed -n ${START},${END}p $sampset | while read INFILE
do
	NAME=$(basename $INFILE .proteins.fa)
	echo "$NAME"
    OUT=$OUTDIR/$NAME.${OUTEXT}
    if [ ! -f ${OUT}.gz ]; then
        time blastp -query $INDIR/$INFILE -db $MEROPS_DB/merops_scan.lib -out ${OUT} \
        -num_threads $CPU -seg yes -soft_masking true \
        -max_target_seqs $MEROPS_MAX_TARGETS \
        -evalue $MEROPS_CUTOFF -outfmt 6 \
        -use_sw_tback
	pigz $OUT
    fi
done
