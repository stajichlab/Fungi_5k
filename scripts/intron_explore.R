#!/usr/bin/env Rscript

library(tidyverse)
library(duckdb)
library(purrr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

library(emojifont)

set.seed(1234)
symbls <- c('fa-github', 'fa-binoculars', 'fa-twitter', 'fa-android', 'fa-coffee', 
            'fa-cube', 'fa-ambulance','fa-check','fa-cutlery','fa-cogs','fa-dot-circle-o','fa-car',
            'fa-building','fa-fire', 'fa-flag','fa-female','fa-gratipay','fa-heart','fa-magnet',
            'fa-lock','fa-map','fa-puzzle-piece','fa-shopping-cart','fa-star','fa-sticky-note',
            'fa-stop-circle-o','fa-volume-down','fa-anchor', 'fa-beer','fa-book','fa-cloud',
            'fa-comment','fa-eject','fa-chrome','fa-child','fa-bomb', 'fa-certificate',
            'fa-desktop','fa-fire-extinguisher','fa-diamond')

idx <- order(symbls)
fontarray <- fontawesome(symbls)
k <- length(fontarray)

# to use a database file already created by 
con <- dbConnect(duckdb(), dbdir="intronDB/introns.duckdb", read_only = TRUE)

speciessql = "SELECT s.LOCUSTAG, PHYLUM, SUBPHYLUM, CLASS, s.ORDER, GENUS, s.SPECIES, s.STRAIN, TOTAL_LENGTH, N50
FROM species as s, asm_stats as stats
WHERE s.LOCUSTAG = stats.LOCUSTAG"

speciesinfo <- dbGetQuery(con, speciessql)

genelensql="
SELECT LOCUSTAG, avg(length) as mean_protein_length, median(length) as median_protein_length, 
                 avg(cds_length) as mean_cds_length, median(cds_length) as median_cds_length
FROM (SELECT substring(transcript_id,1,8) as LOCUSTAG, p.length, p.length * 3 as cds_length FROM gene_proteins p)
GROUP BY LOCUSTAG"
genelen <- dbGetQuery(con, genelensql)

# (SELECT abs(gene_introns.end - gene_introns.start) as intron_length

# rewrite to add introns per KB 
introncountsql="
SELECT LOCUSTAG, avg(intron_count) as mean_intron_ct, median(intron_count) as median_intron_ct
FROM (SELECT substring(gene_introns.transcript_id,1,8) as LOCUSTAG, max(gene_introns.intron_number + 1) as intron_count,
            gene_introns.transcript_id
      FROM gene_introns 
      GROUP BY gene_introns.transcript_id)
GROUP BY LOCUSTAG"

intronct <- dbGetQuery(con, introncountsql)

## test
testsql = "SELECT max(gene_introns.intron_number + 1) as intron_count, gene_introns.transcript_id
      FROM gene_introns WHERE transcript_id LIKE 'FE0E32D4_006065%'
      GROUP BY gene_introns.transcript_id"
testdat <- dbGetQuery(con, testsql)
testdat
testlensql = "SELECT transcript_id, p.length, p.length * 3 as cds_length FROM gene_proteins p
              WHERE transcript_id LIKE 'FE0E32D4_006065%'"
testlendat <- dbGetQuery(con, testlensql)
testlendat
### 

# intron length calculated
intronlensql="
SELECT LOCUSTAG, avg(intron_length) as mean_intron_len, median(intron_length) as median_intron_len
FROM (SELECT substring(transcript_id,1,8) as LOCUSTAG, abs(gene_introns.end - gene_introns.start) as intron_length
      FROM gene_introns)
GROUP BY LOCUSTAG"

intronlens <- dbGetQuery(con, intronlensql)

# intron frequency calculated
intronfreqsql="
SELECT LOCUSTAG, avg(intronsperkb) as mean_intronsperkb, median(intronsperkb) as median_intronsperkb
FROM (SELECT introns.LOCUSTAG, p.transcript_id, 1000 * IFNULL(intron_count,0) / (3 * p.length) as intronsperkb
  FROM gene_proteins AS p LEFT JOIN
  (SELECT substring(gene_introns.transcript_id,1,8) as LOCUSTAG, max(gene_introns.intron_number + 1) as intron_count,
            gene_introns.transcript_id
      FROM gene_introns 
      GROUP BY gene_introns.transcript_id) as introns 
      ON p.transcript_id = introns.transcript_id)
GROUP BY LOCUSTAG
  "

intronfreq <- dbGetQuery(con, intronfreqsql)


intron_data <- intronfreq %>% left_join(intronlens,by="LOCUSTAG") %>% left_join(speciesinfo,by="LOCUSTAG")


# this is where we have to aggregate for all the individual values



sumbysubphylum_freq <- intronfreq %>% left_join(speciesinfo,by="LOCUSTAG") %>%
  group_by(SUBPHYLUM) %>% summarise(freq_mean = mean(intronperkb),
                                    freq_sd = sd(intronperkb),
                                    freq_N = n(),
                                    freq_se = freq_sd / sqrt(freq_N),
                                    freq_upper_limit = freq_mean + freq_se,
                                    freq_lower_limit = freq_mean - freq_se,
                                    freq_total = sum(intronperkb),
                                    freq_median = median(intronperkb))

intronlensql="SELECT substring(transcript_id,1,8) as LOCUSTAG, 
(gene_introns.end - gene_introns.start) as intron_length
FROM gene_introns"

intronlen <- dbGetQuery(con, intronlensql)


intronlen_sum <- intronlen %>% left_join(speciesinfo,by="LOCUSTAG") %>%
  group_by(SUBPHYLUM) %>% summarise(ilen_mean = mean(intron_length),
                                    ilen_sd = sd(intron_length),
                                    ilen_N = n(),
                                    ilen_se = ilen_sd / sqrt(ilen_N),
                                    ilen_upper_limit = ilen_mean + ilen_se,
                                    ilen_lower_limit = ilen_mean - ilen_se,
                                    ilen_total = sum(intron_length),
                                    ilen_median = median(intron_length))

subphylumdata <- sumbysubphylum_freq %>% 
  inner_join(intronlen_sum,by="SUBPHYLUM") %>% 
  inner_join(speciesinfo,by="SUBPHYLUM")


# Define the number of colors you want
nb.cols <- 13
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

subphylumdata$PHYLUM <- factor(subphylumdata$PHYLUM, levels=fontarray[idx])
p <- ggplot(subphylumdata, 
            mapping = aes(x=freq_mean, y=ilen_mean,label=PHYLUM,
                          color=SUBPHYLUM)) +
  geom_point() +
  geom_errorbar(aes(ymin=ilen_lower_limit, 
                    ymax=ilen_upper_limit,
                    xmin=freq_lower_limit,
                    xmax=freq_upper_limit,
                    ), 
                width=.2, 
                
                ) +
  scale_colour_manual(values = mycolors) +
  scale_fill_manual(values = mycolors) +
  xlab("Mean introns per kb") +
  ylab("Mean intron length") +
  theme_cowplot(12) 

ggsave("plots/intron_size_freq.pdf",p,width=15,height=10)

sumbysubphylum_count <- intronct %>% left_join(speciesinfo,by="LOCUSTAG") %>%
  group_by(SUBPHYLUM) %>% summarise(ct_mean = mean(intron_count),
                                    ct_sd = sd(intron_count),
                                    ct_n = n(),
                                    ct_total = sum(intron_count),
                                    ct_median = median(intron_count))


dbDisconnect(con, shutdown = TRUE)


