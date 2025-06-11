#!/usr/bin/bash -l
#SBATCH -p short -c 2 --mem 24gb

module load duckdb
DB=functionalDB/function.duckdb
duckdb -readonly $DB -c "select SPECIES.PHYLUM, SPECIES.GENUS, SPECIES.SPECIES, startloc, counttargetP from species, (select mean(cleavage_position_start) as startloc,count(*) as counttargetP, species_prefix from targetp group by species_prefix) as cs where SPECIES.LOCUSTAG = cs.species_prefix order by startloc desc;"
