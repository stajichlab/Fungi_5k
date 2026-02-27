#!/usr/bin/bash -l
#SBATCH --mem 96gb -c 48 --out muscle_run.%a.log

module load muscle
module load hmmer
module load parallel

parallel -j 2 muscle -align {} -output {.}.afa ::: $(ls *.fasta)
parallel -j 24 hmmbuild {.}.hmm {} ::: $(ls *.afa)
