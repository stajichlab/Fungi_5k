#!/usr/bin/bash -l
#SBATCH -p short -c 4 --mem 24gb --out logs/prep_bigquery_function.log

module load biopython

./scripts/prep_for_bigquery_load.py
