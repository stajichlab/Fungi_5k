#!/usr/bin/bash -l
#SBATCH -p short -C ryzen --mem 128gb -c 96 -N 1 -n 1 --out logs/load_functionDB.log
module load duckdb
DBDIR=functionalDB
DBNAME=function
mkdir -p $DBDIR
# build species table
duckdb -c "CREATE TABLE IF NOT EXISTS species AS SELECT * FROM read_csv_auto('samples.csv')" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_species_locustag ON species(LOCUSTAG)" $DBDIR/$DBNAME.duckdb

# build asm stats table
duckdb -c "CREATE TABLE IF NOT EXISTS asm_stats AS SELECT * FROM read_csv_auto('bigquery/asm_stats.csv')" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_asmstats_locustag ON asm_stats(LOCUSTAG)" $DBDIR/$DBNAME.duckdb

# build proteins
duckdb -c "CREATE TABLE IF NOT EXISTS gene_proteins AS SELECT *, substring(protein_id,1,8) as locustag FROM read_csv_auto('bigquery/gene_proteins.csv.gz')" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_gene_proteins_locustag ON gene_proteins(locustag)" $DBDIR/$DBNAME.duckdb

# build gene info table
duckdb -c "CREATE TABLE IF NOT EXISTS gene_info AS SELECT *, substring(gene_id,1,8) as locustag FROM read_csv_auto('bigquery/gene_info.csv.gz')" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_gene_info_locustag ON gene_info(locustag)" $DBDIR/$DBNAME.duckdb

# build signalp table
duckdb -c "CREATE TABLE IF NOT EXISTS signalp AS SELECT * FROM read_csv_auto('bigquery/signalp.signal_peptide.csv.gz')" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_signalp_locustag ON signalp(species_prefix)" $DBDIR/$DBNAME.duckdb

# build merops table
duckdb -c "CREATE TABLE IF NOT EXISTS merops AS SELECT * FROM read_csv_auto('bigquery/merops.csv.gz')" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_merops_locustag ON merops(species_prefix)" $DBDIR/$DBNAME.duckdb

# build CAZY table
duckdb -c "CREATE TABLE IF NOT EXISTS cazy_overview AS SELECT * FROM read_csv_auto('bigquery/cazy.overview.csv.gz')" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_cazy_overview_locustag ON cazy_overview(species_prefix)" $DBDIR/$DBNAME.duckdb

# build CAZY domains table
duckdb -c "CREATE TABLE IF NOT EXISTS cazy AS SELECT * FROM read_csv_auto('bigquery/cazy.cazymes_hmm.csv.gz')" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_cazy_locustag ON cazy(species_prefix)" $DBDIR/$DBNAME.duckdb

# build Pfam domains table
duckdb -c "CREATE TABLE IF NOT EXISTS pfam AS SELECT *, substring(protein_id,1,8) as species_prefix FROM read_csv_auto('bigquery/pfam.csv.gz')" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_pfam_locustag ON pfam(species_prefix)" $DBDIR/$DBNAME.duckdb
