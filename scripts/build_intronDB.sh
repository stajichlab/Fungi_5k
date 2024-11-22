#!/usr/bin/bash -l
#SBATCH -p short -C ryzen --mem 128gb -c 96 -N 1 -n 1 --out logs/load_intronDB.log
module load duckdb
DBDIR=intronDB
mkdir -p $DBDIR
# build species table
duckdb -c "CREATE TABLE IF NOT EXISTS species AS SELECT * FROM read_csv_auto('samples.csv')" $DBDIR/introns.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_species_locustag ON species(LOCUSTAG)" $DBDIR/introns.duckdb

# build asm stats table
duckdb -c "CREATE TABLE IF NOT EXISTS asm_stats AS SELECT * FROM read_csv_auto('bigquery/asm_stats.csv')" $DBDIR/introns.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_asmstats_locustag ON asm_stats(LOCUSTAG)" $DBDIR/introns.duckdb

# build introns to gene
duckdb -c "CREATE TABLE IF NOT EXISTS gene_introns AS SELECT *, substring(intron_id,1,8) as locustag FROM read_csv_auto('bigquery/gene_introns.csv.gz')" $DBDIR/introns.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_gene_introns_locustag ON gene_introns(locustag)" $DBDIR/introns.duckdb

# build transcript table
duckdb -c "CREATE TABLE IF NOT EXISTS gene_transcripts AS SELECT *, substring(gene_id,1,8) as locustag FROM read_csv_auto('bigquery/gene_transcripts.csv.gz')" $DBDIR/introns.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_gene_transcripts_locustag ON gene_transcripts(locustag)" $DBDIR/introns.duckdb

# build exon table
duckdb -c "CREATE TABLE IF NOT EXISTS gene_exons AS SELECT *, substring(transcript_id,1,8) as locustag FROM read_csv_auto('bigquery/gene_exons.csv.gz')" $DBDIR/introns.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_gene_exon_locustag ON gene_exons(locustag)" $DBDIR/introns.duckdb

# build proteins
duckdb -c "CREATE TABLE IF NOT EXISTS gene_proteins AS SELECT *, substring(protein_id,1,8) as locustag FROM read_csv_auto('bigquery/gene_proteins.csv.gz')" $DBDIR/introns.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_gene_proteins_locustag ON gene_proteins(locustag)" $DBDIR/introns.duckdb

# build gene info table
duckdb -c "CREATE TABLE IF NOT EXISTS gene_info AS SELECT *, substring(gene_id,1,8) as locustag FROM read_csv_auto('bigquery/gene_info.csv.gz')" $DBDIR/introns.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_gene_info_locustag ON gene_info(locustag)" $DBDIR/introns.duckdb




duckdb -c "CREATE TABLE IF NOT EXISTS species AS FROM read_csv_auto('samples.csv')" $DBDIR/intron_matches.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_species_locustag ON species(LOCUSTAG)" $DBDIR/intron_matches.duckdb
duckdb -c "CREATE TABLE IF NOT EXISTS intron_matches AS 
SELECT *, substring(query,1,8) as qlocus, substring(subject,1,8) as slocus  
FROM read_csv('fungi_introns/blastn_results/*.gz',sep='\t',header=false,compression='gzip', 
names=['query', 'qlength', 'subject', 'slength', 'identity', 'alnlen', 'mismatches', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
WHERE identity >= 90 AND query != subject
" $DBDIR/intron_matches.duckdb