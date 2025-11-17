#!/usr/bin/env Rscript

library(tidyverse)
library(duckdb)
library(purrr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

plotdir = "plots/mmseqs_clusters"
dir.create(file.path(plotdir),showWarnings = FALSE)
pdf(file.path(plotdir,"clusters_mmseqs.pdf"))

# to use a database file already created by
con <- dbConnect(duckdb(), dbdir="functionalDB/function.duckdb", read_only = TRUE)

# build a temp table which has the number of non-pfam annotated genes
cluster_create_sql = "
CREATE TEMPORARY TABLE mmseq_cluster_pfam AS SELECT ct.*,cl.transcript_id, gp.length, pfam.pfam_id, pfam_len, full_seq_e_value
FROM mmseqs_orthogroup_cluster_count ct,
mmseqs_orthogroup_clusters cl, 
gene_proteins gp
LEFT JOIN pfam ON pfam.protein_id = gp.transcript_id
WHERE cl.orthogroup = ct.orthogroup AND
      ct.group_count > 100 AND
      cl.transcript_id = gp.transcript_id 
order by group_count"

dbExecute(con, cluster_create_sql)

nulltable_sql = "
CREATE TEMPORARY TABLE mmseq_cluster_pfam_nullcount AS 
SELECT orthogroup, COUNT(*) as pfam_null, group_count, pfam_null / group_count as null_ratio 
FROM mmseq_cluster_pfam where pfam_id IS NULL 
GROUP BY orthogroup, group_count ORDER BY null_ratio DESC,group_count DESC"

dbExecute(con,nulltable_sql)

queryproteins_sql = 
"SELECT gp.transcript_id, gp.peptide, cl.orthogroup, nc.group_count, nc.null_ratio
FROM mmseq_cluster_pfam_nullcount nc,
mmseqs_orthogroup_clusters cl,
gene_proteins gp
WHERE nc.null_ratio > 0.95 AND
cl.transcript_id = gp.transcript_id AND
nc.orthogroup = cl.orthogroup
ORDER BY nc.group_count DESC
"

cluster_info <- dbGetQuery(con,queryproteins_sql)
cluster_info 

# now let us dump out a protein set for each of the 
write_csv(cluster_info,file="results/mmseqs_nopfam_genedump.csv")


dbDisconnect(con, shutdown = TRUE)
