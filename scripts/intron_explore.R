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

genelensql="SELECT substring(transcript_id,1,8) as LOCUSTAG, p.length as protein_length, p.length * 3 as cds_length, p.transcript_id
FROM gene_proteins p"

genelen <- dbGetQuery(con, genelensql)

intronsql="SELECT substring(transcript_id,1,8) as LOCUSTAG, 
max(gene_introns.intron_number + 1) as intron_count, transcript_id
FROM gene_introns 
GROUP BY transcript_id"

intronct <- dbGetQuery(con, intronsql)

intronfreq = genelen %>% left_join(intronct %>% select(-c(LOCUSTAG)),by="transcript_id") %>% 
  mutate( intronperkb = 1000 * replace_na(intron_count / cds_length,0)) 
  
speciessql = "SELECT s.LOCUSTAG, PHYLUM, SUBPHYLUM, CLASS, TOTAL_LENGTH 
FROM species as s, asm_stats as stats
WHERE s.LOCUSTAG = stats.LOCUSTAG"
speciesinfo <- dbGetQuery(con, speciessql)

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

p

sumbysubphylum_count <- intronct %>% left_join(speciesinfo,by="LOCUSTAG") %>%
  group_by(SUBPHYLUM) %>% summarise(ct_mean = mean(intron_count),
                                    ct_sd = sd(intron_count),
                                    ct_n = n(),
                                    ct_total = sum(intron_count),
                                    ct_median = median(intron_count))


dbDisconnect(con, shutdown = TRUE)


