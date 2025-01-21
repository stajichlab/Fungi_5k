#!/usr/bin/bash -l
#SBATCH -p short -C xeon --mem 256gb -c 96 --out logs/build_duckdb.log

module load biopython

./scripts/identical_intron_search.py --input "blastn_results/*.tab.gz" --db ./intronblast.db
