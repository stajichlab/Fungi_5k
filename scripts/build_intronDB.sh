#!/usr/bin/bash -l
#SBATCH -p short -C ryzen --mem 96gb -c 24 -N 1 -n 1 --out logs/load_intronDB.log
module load duckdb
DBDIR=intronDB
DBNAME=intron_db
mkdir -p $DBDIR

duckdb -c "CREATE TABLE IF NOT EXISTS species AS SELECT * FROM read_csv_auto('samples.csv')" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE UNIQUE INDEX IF NOT EXISTS idx_species_locustag ON species(LOCUSTAG)" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE UNIQUE INDEX IF NOT EXISTS idx_species_asm ON species(ASMID)" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_species_speciesin ON species(SPECIESIN)" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_species_species ON species(SPECIES,GENUS)" $DBDIR/$DBNAME.duckdb
# build asm stats table
duckdb -c "CREATE TABLE IF NOT EXISTS asm_stats AS SELECT * FROM read_csv_auto('bigquery/asm_stats.csv.gz')" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE UNIQUE INDEX IF NOT EXISTS idx_stats_locustag ON asm_stats(LOCUSTAG)" $DBDIR/$DBNAME.duckdb

# build chrom stats table
duckdb -c "CREATE TABLE IF NOT EXISTS chrom_info AS SELECT * FROM read_csv_auto('bigquery/chrom_info.csv.gz')" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_chrominfo_locustag ON chrom_info(LOCUSTAG)" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE UNIQUE INDEX IF NOT EXISTS idx_chrominfo_locustag_chrom ON chrom_info(LOCUSTAG,chrom_name)" $DBDIR/$DBNAME.duckdb

# build proteins
duckdb -c "CREATE TABLE IF NOT EXISTS gene_proteins AS SELECT *, substring(protein_id,1,8) as locustag FROM read_csv_auto('bigquery/gene_proteins.csv.gz')" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_gene_proteins_locustag ON gene_proteins(locustag)" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_gene_proteins_gene_id ON gene_proteins(gene_id)" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE UNIQUE INDEX IF NOT EXISTS idx_gene_proteins_protein_id ON gene_proteins(protein_id)" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_gene_proteins_tx_id ON gene_proteins(transcript_id)" $DBDIR/$DBNAME.duckdb

# Add tRNA
duckdb -c "CREATE TABLE IF NOT EXISTS gene_trna AS SELECT *, substring(gene_id,1,8) as locustag FROM read_csv_auto('bigquery/gene_trnas.csv.gz')" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_gene_trna_locustag ON gene_trna(locustag)" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE UNIQUE INDEX IF NOT EXISTS idx_gene_trna_gene_id ON gene_trna(gene_id)" $DBDIR/$DBNAME.duckdb

# Add transcripts
duckdb -c "CREATE TABLE IF NOT EXISTS gene_transcripts AS SELECT *, substring(gene_id,1,8) as locustag FROM read_csv_auto('bigquery/gene_transcripts.csv.gz')" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_gene_tx_transcripts_locustag ON gene_transcripts(locustag)" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_gene_tx_transcripts_gene_id ON gene_transcripts(gene_id)" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE UNIQUE INDEX IF NOT EXISTS idx_gene_tx_transcript_id ON gene_transcripts(transcript_id)" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_gene_tx_chrom ON gene_transcripts(chrom)" $DBDIR/$DBNAME.duckdb

# build gene info table
duckdb -c "CREATE TABLE IF NOT EXISTS gene_info AS SELECT *, substring(gene_id,1,8) as locustag FROM read_csv_auto('bigquery/gene_info.csv.gz')" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_gene_info_locustag ON gene_info(locustag)" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE UNIQUE INDEX IF NOT EXISTS idx_gene_info_gene_id ON gene_info(gene_id)" $DBDIR/$DBNAME.duckdb


# build intron tables
duckdb -c "CREATE TABLE IF NOT EXISTS intron_matches AS 
SELECT *, substring(query,1,8) as qlocus, substring(subject,1,8) as slocus  
FROM read_csv('fungi_introns/blastn_results/*.gz',sep='\t',header=false,compression='gzip', 
names=['query', 'qlength', 'subject', 'slength', 'identity', 'alnlen', 'mismatches', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
WHERE identity >= 80 AND query != subject"  $DBDIR/$DBNAME.duckdb

duckdb -c "CREATE INDEX IF NOT EXISTS idx_intron_matches_qlocus ON intron_matches(qlocus,query)" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_intron_matches_slocus ON intron_matches(slocus,subject)"  $DBDIR/$DBNAME.duckdb
