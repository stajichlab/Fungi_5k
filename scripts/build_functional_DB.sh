#!/usr/bin/bash -l
#SBATCH -p short -C ryzen --mem 128gb -c 96 -N 1 -n 1 --out logs/load_functionDB.log
module load duckdb
DBDIR=functionalDB
DBNAME=function
mkdir -p $DBDIR
# build species table
duckdb -c "CREATE TABLE IF NOT EXISTS species AS SELECT * FROM read_csv_auto('samples.csv')" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE UNIQUE INDEX IF NOT EXISTS idx_species_locustag ON species(LOCUSTAG)" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE UNIQUE INDEX IF NOT EXISTS idx_species_asm ON species(ASMID)" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_species_speciesin ON species(SPECIESIN)" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_species_species ON species(SPECIES,GENUS)" $DBDIR/$DBNAME.duckdb
# build asm stats table
duckdb -c "CREATE TABLE IF NOT EXISTS asm_stats AS SELECT * FROM read_csv_auto('bigquery/asm_stats.csv.gz')" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE UNIQUE INDEX IF NOT EXISTS idx_asmstats_locustag ON asm_stats(LOCUSTAG)" $DBDIR/$DBNAME.duckdb

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

# build signalp table
duckdb -c "CREATE TABLE IF NOT EXISTS signalp AS SELECT * FROM read_csv_auto('bigquery/signalp.signal_peptide.csv.gz')" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_signalp_locustag ON signalp(species_prefix)" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_signalp_protein_id ON signalp(protein_id)" $DBDIR/$DBNAME.duckdb

# build merops table
duckdb -c "CREATE TABLE IF NOT EXISTS merops AS SELECT * FROM read_csv_auto('bigquery/merops.csv.gz')" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_merops_locustag ON merops(species_prefix)" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_merops_protein_id ON merops(protein_id)" $DBDIR/$DBNAME.duckdb

# build CAZY table
duckdb -c "CREATE TABLE IF NOT EXISTS cazy_overview AS SELECT * FROM read_csv_auto('bigquery/cazy.overview.csv.gz')" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_cazy_overview_locustag ON cazy_overview(species_prefix)" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_cazy_overview_protein_id ON cazy_overview(protein_id)" $DBDIR/$DBNAME.duckdb

# build CAZY domains table
duckdb -c "CREATE TABLE IF NOT EXISTS cazy AS SELECT * FROM read_csv_auto('bigquery/cazy.cazymes_hmm.csv.gz')" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_cazy_locustag ON cazy(species_prefix)" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_cazy_protein_id ON cazy(protein_id)" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_cazy_HMM ON cazy(HMM_id)" $DBDIR/$DBNAME.duckdb


# build Pfam domains table
duckdb -c "CREATE TABLE IF NOT EXISTS pfam AS SELECT *, substring(protein_id,1,8) as species_prefix FROM read_csv_auto('bigquery/pfam.csv.gz')" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_pfam_locustag ON pfam(species_prefix)" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_pfam_protein_id ON pfam(protein_id)" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_pfam_protein_id ON pfam(pfam_id)" $DBDIR/$DBNAME.duckdb


# Add funguild
duckdb -c "CREATE TABLE IF NOT EXISTS funguild AS SELECT * FROM read_csv_auto('bigquery/species_funguild.csv.gz')" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_funguild_locustag ON funguild(species_prefix)" $DBDIR/$DBNAME.duckdb

# Add codon freq
duckdb -c "CREATE TABLE IF NOT EXISTS codon_frequency AS SELECT * FROM read_csv_auto('bigquery/codon_freq.csv.gz')" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_codonfreq_locustag ON codon_frequency(species_prefix)" $DBDIR/$DBNAME.duckdb

# Add AA freq
duckdb -c "CREATE TABLE IF NOT EXISTS aa_frequency AS SELECT * FROM read_csv_auto('bigquery/aa_freq.csv.gz')" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_aa_locustag ON aa_frequency(species_prefix)" $DBDIR/$DBNAME.duckdb

# Add gene distance
duckdb -c "CREATE TABLE IF NOT EXISTS gene_pairwise_distances AS SELECT * FROM read_csv_auto('bigquery/gene_pairwise_distances.csv.gz')" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_gene_pw_locustag ON gene_pairwise_distances(species_prefix)" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE UNIQUE INDEX IF NOT EXISTS idx_gene_pw_genes ON gene_pairwise_distances(left_gene,right_gene)" $DBDIR/$DBNAME.duckdb

# add tmmhmm
duckdb -c "CREATE TABLE IF NOT EXISTS tmhmm AS SELECT * FROM read_csv('bigquery/tmhmm.csv.gz')" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_tmhmm_locus ON tmhmm(species_prefix)" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_tmhmm_protein_id ON tmhmm(protein_id)" $DBDIR/$DBNAME.duckdb


# add pscan
duckdb -c "CREATE TABLE IF NOT EXISTS prosite AS SELECT * FROM read_csv('bigquery/ps_scan.csv.gz')" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_prosite_locus ON prosite(species_prefix)" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_prosite_protein_ud ON prosite(protein_id)" $DBDIR/$DBNAME.duckdb

# add orthogroups
duckdb -c "CREATE TABLE IF NOT EXISTS mmseqs_orthogroup_clusters AS SELECT * FROM read_csv('bigquery/mmseqs_orthogroup_clusters.csv.gz')" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_og ON mmseqs_orthogroup_clusters(orthogroup)" $DBDIR/$DBNAME.duckdb
duckdb -c "CREATE INDEX IF NOT EXISTS idx_tid ON mmseqs_orthogroup_clusters(transcript_id)" $DBDIR/$DBNAME.duckdb
