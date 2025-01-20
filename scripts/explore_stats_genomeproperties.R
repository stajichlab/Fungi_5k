#!/usr/bin/env Rscript

library(tidyverse)
library(duckdb)
library(purrr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(paletteer)
library(cowplot)

DBDIR="functionalDB"
DBNAME="function.duckdb"
plotdir = "plots"
dir.create(file.path(plotdir),showWarnings = FALSE)
statsplotdir = file.path(plotdir,"asm_seqstats")
dir.create(statsplotdir, showWarnings = FALSE)

# to use a database file already created by 
con <- dbConnect(duckdb(), dbdir=file.path(DBDIR,DBNAME), read_only = TRUE)

asmstat_sql ="
SELECT sp.*, stats.GC_PERCENT, stats.TOTAL_LENGTH, gene_count, mean_gene_length
FROM 
species as sp,
asm_stats as stats,
(SELECT LOCUSTAG, count(*) as gene_count, MEAN(length) as mean_gene_length
FROM gene_proteins 
GROUP BY LOCUSTAG) as gp
WHERE 
sp.LOCUSTAG = stats.LOCUSTAG and gp.LOCUSTAG = sp.LOCUSTAG"


asmstat_res <- dbGetQuery(con, asmstat_sql)
head(asmstat_res)

# barplot for these GENERA
sumAsm = asmstat_res %>% group_by(ORDER) %>% 
  summarize(
  gc_mean = mean(GC_PERCENT),
  gc_N = n(),
  gc_sd = sd(GC_PERCENT),
  gc_se = gc_sd / sqrt(gc_N),
  len_mean = mean(TOTAL_LENGTH / 1000000),
  len_N = n(),
  len_sd = sd(TOTAL_LENGTH/1000000),
  len_se = len_sd / sqrt(len_N))


len_p <- ggplot(sumAsm) + geom_bar(
                         aes(y=len_mean,
                             x=ORDER),
                         stat="identity", width=0.5) +
  geom_errorbar(
                aes(y=len_mean,
                    x=ORDER,
                    ymin = len_mean - len_se,
                    ymax = len_mean + len_se),
                stat="identity", width=0.5) +
  geom_point(aes(y=len_mean, x=ORDER),size=2) + 
  scale_colour_manual(values = mycolors) +
  scale_fill_manual(values = mycolors) +
  xlab("ORDER") +
  ylab("Genome size (Mb)") +
  theme_cowplot(12) + theme(axis.text.x = element_text(angle = 90,size = 6))

len_p

ggsave(file.path(statsplotdir,"genome_size_order.pdf"),len_p,width=24,height=10)

nb.cols <- length(unique(asmstat_res$SUBPHYLUM))
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)

count_glen_p <- ggplot(asmstat_res) +
  geom_point(aes(x=gene_count, y=mean_gene_length,color=SUBPHYLUM, fill=SUBPHYLUM),size=1.5,alpha=0.5) + 
  scale_colour_manual(values = mycolors) +
  scale_fill_manual(values = mycolors) +
  ylab("Gene Length") +
  xlab("Gene count (Mb)") + ggtitle("Gene Count vs Mean Gene length") +
  theme_cowplot(12) + theme(legend.position="bottom") + scale_x_log10() 
count_glen_p
ggsave(file.path(statsplotdir,"genecount_gene_length.pdf"),count_glen_p,width=15,height=10)

phylum_countglen_p <- ggplot(asmstat_res) +
  geom_point(aes(x=gene_count, y=mean_gene_length,color=SUBPHYLUM, fill=SUBPHYLUM),size=1.5,alpha=0.5) + 
  scale_colour_manual(values = mycolors) +
  scale_fill_manual(values = mycolors) +
  ylab("Gene Length") +
  xlab("Gene count (Mb)") + ggtitle("Gene Count vs Mean Gene length") +
  theme_cowplot(12) + theme(legend.position="bottom") + scale_x_log10() + facet_wrap(~PHYLUM)

phylum_countglen_p
ggsave(file.path(statsplotdir,"genecount_gene_length_facet.pdf"),count_glen_p,width=15,height=10)

basidio <- asmstat_res %>% filter(PHYLUM=="Basidiomycota")
nb.cols <- length(unique(basidio$CLASS))
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)

basidio_countglen_p <- ggplot(basidio) +
  geom_point(aes(x=gene_count, y=mean_gene_length,color=CLASS, fill=CLASS),size=1.5,alpha=0.7) + 
  scale_colour_manual(values = mycolors) +
  scale_fill_manual(values = mycolors) +
  ylab("Gene Length") +
  xlab("Gene count (Mb)") + ggtitle("Gene Count vs Mean Gene length") +
  theme_cowplot(12) + theme(legend.position="bottom") + scale_x_log10() + facet_wrap(~SUBPHYLUM)

basidio_countglen_p
ggsave(file.path(statsplotdir,"genecount_gene_length_basidio.pdf"),count_glen_p,width=15,height=10)

# could make this a function + lapply
for (subphylum in unique(basidio$SUBPHYLUM))
{
  subph <- basidio %>% filter(SUBPHYLUM==subphylum)

  countglen_p <- ggplot(subph) +
  geom_point(aes(x=gene_count, y=mean_gene_length,color=CLASS, fill=CLASS),size=1.5,alpha=0.7) + 
    scale_colour_brewer(palette = "Set1") +
  ylab("Gene Length") +
  xlab("Gene count (Mb)") + ggtitle(sprintf("%s: Gene Count vs Mean Gene length",subphylum)) +
  theme_cowplot(12) + theme(legend.position="bottom") + scale_x_log10()

  countglen_p
  ggsave(file.path(statsplotdir,sprintf("genecount_gene_length_basidio_%s.pdf",subphylum)),countglen_p,width=8,height=8)
}

asco <- asmstat_res %>% filter(PHYLUM=="Ascomycota") %>% filter(SUBPHYLUM != "Ascomycotina" & 
                                                                  SUBPHYLUM != "NA") 
for (subphylum in unique(asco$SUBPHYLUM))
{
  subph <- asco %>% filter(SUBPHYLUM==subphylum)
  nb.cols <- length(unique(subph$CLASS))
  mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)
  countglen_p <- ggplot(subph) +
    geom_point(aes(x=gene_count, y=mean_gene_length,color=CLASS, fill=CLASS),size=1.5,alpha=0.7) + 
    scale_colour_manual(values = mycolors) +
    scale_fill_manual(values = mycolors) +
    ylab("Gene Length") +
    xlab("Gene count (Mb)") + ggtitle(sprintf("%s: Gene Count vs Mean Gene length",subphylum)) +
    theme_cowplot(12) + theme(legend.position="bottom") + scale_x_log10()
  
  countglen_p
  ggsave(file.path(statsplotdir,sprintf("genecount_gene_length_asco_%s.pdf",subphylum)),countglen_p,width=8,height=8)
}

dbDisconnect(con, shutdown = TRUE)


