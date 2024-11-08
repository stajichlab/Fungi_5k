#!/usr/bin/env Rscript

library(tidyverse)
library(vroom)
library(dplyr)
library(ggplot2)
library(cowplot)
library(randomForest)
library(ggfortify)
library(RColorBrewer)
# pak::pak("r-dbi/bigrquery")
library(bigrquery)
library(DBI)
library(collapse)
library(ggpubfigs)

cbPalette <- c( "#E69F00", "#999999","#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00","#CC79A7","#000000")


"019E1F-BADE9B-60248B"
project="ucr-ursa-major-stajich-lab"
con <- dbConnect(
  bigquery(),
  project = "ucr-ursa-major-stajich-lab",
  dataset = "fungi5k",
  billing = "ucr-ursa-major-stajich-lab" 
)

dbListTables(con)

aasql = "SELECT s.PHYLUM, s.SUBPHYLUM, s.CLASS, s.SUBCLASS, s.ORDER, s.FAMILY, s.GENUS, s.SPECIES, s.NCBI_TAXONID, s.LOCUSTAG, aa_freq.amino_acid, aa_freq.frequency
FROM 
 `ucr-ursa-major-stajich-lab.fungi5k.species` as s,
  `ucr-ursa-major-stajich-lab.fungi5k.aa_freq` as aa_freq
  
WHERE
  aa_freq.species_prefix = s.LOCUSTAG"

aafreq <- dbGetQuery(con, aasql)


wide_aa_freq <- aafreq %>% 
  filter( ! is.na(frequency)) %>% 
  pivot_wider(names_from = amino_acid, values_from=frequency) %>% 
  filter(! is.na(PHYLUM) ) %>%
  select(-c(LOCUSTAG))

wide_aa_freq$PHYLUM = factor(wide_aa_freq$PHYLUM) 


# build a random training set?
set.seed(222)
ind <- sample(2, nrow(wide_aa_freq), replace = TRUE, prob = c(0.7, 0.3))
train <- wide_aa_freq[ind==1,]
test <- wide_aa_freq[ind==2,]

train$PHYLUM <- factor(train$PHYLUM)
test$PHYLUM <- factor(test$PHYLUM)


lapply(train, function(x) any(is.na(x)))

rf <- randomForest(PHYLUM~., data=train, proximity=TRUE)

MDSplot(rf, train$PHYLUM, pch=as.numeric(train$PHYLUM))
aacols = colnames(wide_aa_freq %>% select(-c('PHYLUM','SUBPHYLUM','CLASS', 'SUBCLASS','ORDER','FAMILY','GENUS','SPECIES','NCBI_TAXONID')) )

pca <- prcomp(wide_aa_freq %>% select(all_of(aacols)),scale = TRUE)
summary(pca)
apply(wide_aa_freq %>% select(all_of(aacols)), 2, var)
pvar = round(summary(pca)$importance[2, 1:2], 2)

phylum_pcaplot <- autoplot(pca, data = wide_aa_freq, colour = 'PHYLUM') + 
  scale_colour_brewer(palette = "Paired") + theme_cowplot(12) + 
  ggtitle("PCA plot of AA frequency by Phylum") 

phylum_pcaplot
ggsave("pca_plot_AA_phylum_freq.pdf",phylum_pcaplot,width=15,height=15)

ascoAA =  wide_aa_freq %>% filter(PHYLUM == "Ascomycota")

nb.cols <- 24
mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(nb.cols)

ascopca <- prcomp(ascoAA %>% select(all_of(aacols)),scale = TRUE)
asco_pcaplot <- autoplot(ascopca, data = ascoAA, colour = 'CLASS') + 
   theme_cowplot(12) + scale_color_manual(values = mycolors) +
  ggtitle("PCA plot of Ascomycota AA frequency by Class") 

asco_pcaplot
ggsave("pca_plot_AA_ascomycota_class_freq.pdf",asco_pcaplot,width=15,height=15)


noSacchAA =  wide_aa_freq %>% filter(PHYLUM == "Ascomycota" & CLASS != 'Saccharomycetes')

nb.cols <- 20
mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(nb.cols)

pca <- prcomp(noSacchAA %>% select(all_of(aacols)),scale = TRUE)
noSacchAA_pcaplot <- autoplot(pca, data = noSacchAA, colour = 'CLASS') + 
  theme_cowplot(12) + scale_color_manual(values = mycolors) +
  ggtitle("PCA plot of non-Saccharomycetes Ascomycota AA frequency by Class") 

noSacchAA_pcaplot
ggsave("pca_plot_AA_ascomycota_nonsacch_freq.pdf",noSacchAA_pcaplot,width=15,height=15)


DothidAA =  wide_aa_freq %>% filter(CLASS == 'Dothideomycetes')

nb.cols <- 70
mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(nb.cols)

pca <- prcomp(DothidAA %>% select(all_of(aacols)),scale = TRUE)
dothidAA_pcaplot <- autoplot(pca, data = DothidAA, colour = 'FAMILY') + 
  theme_cowplot(12) + scale_color_manual(values = mycolors) +
  ggtitle("PCA plot of Dothideomycetes Ascomycota AA frequency") 

dothidAA_pcaplot
ggsave("pca_plot_AA_Dothideomycetes_freq.pdf",dothidAA_pcaplot,width=15,height=15)


ggplot(sumfreq) + geom_bar( aes(y=PHYLUM, x=MEAN,color=amino_acid,fill=amino_acid), 
                            stat="identity", alpha=0.7) 

p <- ggplot(data = z.points, aes(x = MDS1, y = MDS2)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

p + geom_point() 


aa_freq2 <- aafreq %>% 
  filter( ! is.na(frequency)) %>% filter(amino_acid != "X") %>%
  filter(! is.na(PHYLUM) ) %>% select(-c(LOCUSTAG)) %>% group_by(PHYLUM,amino_acid) %>% select(c(PHYLUM,amino_acid,frequency)) 

means <- aa_freq2 %>% fmean
medians <- aa_freq2 %>% fmedian
SDs <- aa_freq2 %>% fsd

sumAA_table_phyla <- aa_freq2 %>% summarise(
  count = n(),
  mean = mean(frequency, na.rm = TRUE),
  sd = sd(frequency, na.rm = TRUE)
)

mycolors <- friendly_pal("muted_nine")

p <- ggplot(aa_freq2) + geom_boxplot(aes(x = PHYLUM,y=frequency,color=PHYLUM)) + facet_wrap(. ~ amino_acid,nrow = 4) + 
  scale_color_manual(values = mycolors) +
  scale_fill_manual(values = mycolors) + theme_cowplot(12) + 
  ggtitle("AA frequency by Phylum") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p


res.aov <- aov(frequency ~ amino_acid, data = aa_freq2)
summary(res.aov)
TukeyHSD(res.aov)

aa_freq3 <- aafreq %>% 
  filter( ! is.na(frequency)) %>% filter(amino_acid != "X") %>%
  filter(! is.na(CLASS) ) %>% select(-c(LOCUSTAG)) %>% group_by(CLASS,amino_acid) %>% select(c(CLASS,amino_acid,frequency)) %>% filter(amino_acid == "C")

res.aov <- aov(frequency ~ CLASS, data = aa_freq3)
summary(res.aov)
TukeyHSD(res.aov)
