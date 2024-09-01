##### set up ####
library(phyloseq)
library(tidyverse)
library(zoo)

setwd("C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Seagrass_Montague_2021-03-07/galliano_seagrass_2021and2023/imput_data")


uf = readRDS("C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Seagrass_Montague_2021-03-07/galliano_seagrass_2021and2023/imput_data/MH_16S_unfiltered_phyloseq.rds")

## for 18S filtering
uf = subset_taxa(uf, Kingdom !="Unassigned" & 
                   Kingdom !="Eukaryota"&
                   Family!="Mitochondria" &
                   Order !="Chloroplast")



##### remove samples with less than 1000 reads #####
uf@sam_data$rd_uf = sample_sums(uf)

uf = subset_samples(uf, rd_uf>=1000)
uf

##### FILTERING - REMOVE INDIVIDUAL ASVS WITH LESS THAN 0.001% of reads #####
## extract OTU dataframe from phyloseq object
otu.pruned <- as.data.frame(t(as.matrix(otu_table(uf))))

## remove OTU (rows) with less than 100 reads accross whole dataset but keep all samples
## make sure asv sequence is rownames and sample id is column name
otu.pruned$rowsum = rowSums(otu.pruned)

## remove low frequency asvs (less than 0.001%)
total_asvs = sum(otu.pruned$rowsum)
otu.pruned$total = total_asvs
otu.pruned$percent_abundance = otu.pruned$rowsum/otu.pruned$total
otu.pruned = subset(otu.pruned, otu.pruned$percent_abundance>0.00001)

## remove rowsum column
otu = subset(otu.pruned, select=-c(rowsum, total, percent_abundance))


##### remove ASVs found in two samples or less #####
richness = function(x){
  return(sum(x>0))}
## calculate richness on entire dataframe
otu$richness = apply(otu,1,richness) # use all columns of otu dataframe
summary(otu$richness)


## remove OTU (rows) with richness less than 5
otu = subset(otu, otu$richness>2)
## check that it worked (min richness should be 2 or higher)
summary(otu$richness)
## remove richness column
otu = subset(otu, select=-c(richness))
## see which OTUs were lost
#otu.lost = subset(otu.pruned, otu.pruned$richness<=1)
## save file of lost OTUs



##### MAKE ALL THE CELLS IN THE OTU TABLE WITH VALUES 9 OR LESS 0 #####
otu <- mutate_all(otu, funs(ifelse(. < 9, 0, .)))


##### make filtered phyloseq object ####
filt16 = phyloseq(sample_data(uf@sam_data),
                  tax_table(uf@tax_table),
                  otu_table(as.matrix(otu), taxa_are_rows = T))
filt16@sam_data$rd_filt = sample_sums(filt16)
View(filt16@sam_data)


write_rds(filt16, "MH_16S_filtered001percent_notrarefied_phyloseq.rds")

