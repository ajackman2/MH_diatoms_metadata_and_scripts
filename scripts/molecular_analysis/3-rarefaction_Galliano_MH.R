##### set up ####
library(phyloseq)
library(metagMisc)
library(iNEXT)
library(vegan)
library(tidyverse)
library(ggplot2)

#setwd("C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Seagrass_Montague_2021-03-07/galliano_seagrass_2021and2023/imput_data")
setwd("/parfreylab/schenk/seagrass_montague_processing")
  
phylo16 = readRDS("MH_16S_filtered001percent_notrarefied_phyloseq.rds")
phylo18 = readRDS("MH_18S_filtered001percent_notrarefied_phyloseq.rds")
phylorb = readRDS("MHandGalliano_RBCL_filtered_phyloseq.RDS")

##### 16S #####
##### look at minimum, mean, and maximum sample counts 
smin <- min(sample_sums(phylo16)) 
meanreads <-mean(sample_sums(phylo16)) 
smax <- max(sample_sums(phylo16))
totalreads <- sum(sample_sums(phylo16))

##### USE iNEXT TO FIND SAMPLE COVERAGE
### check if taxa are rows in the phyloseq object
taxa_are_rows(phylo16)
## mine is TRUE

### prepare otu table for iNEXT function
x <- metagMisc::prepare_inext(
  as.data.frame(otu_table(phylo16)),
  correct_singletons = T)

## calculate coverage with iNEXT Char.Ind
SC <- plyr::llply(.data = x, .fun = function(z){ try( iNEXT:::Chat.Ind(z, sum(z)) ) })
plyr::ldply(.data = SC, .fun = class)

## get dataframe to sort coverage and decide what to pick 
coverage = as.data.frame(t(as.matrix(as.data.frame(SC))))

#### run coveraged-based  rarefaction 
COVERAGE_RAREF_1000 <- phyloseq_coverage_raref(physeq=phylo16,
                                               coverage = 0.8,
                                               iter = 1000,
                                               replace = F, 
                                               correct_singletons = T,
                                               drop_lowcoverage = F,
                                               multithread = T) 


write_rds(COVERAGE_RAREF_1000, "COVERAGE_RAREF_MH_phylo16_phyloseq.RDS")

### select to to the entire 1000 phyloseq objects
all_1000_phyloseq <- COVERAGE_RAREF_1000[c(1:1000)]

### first, extract otu tables from the 1000 phyloseq objects
otu_tables_1000 <- lapply(all_1000_phyloseq, function(z) as.data.frame(t(phyloseq::otu_table(z))))

### average all matrices to get the mean abundance across all iterations
average_otu_tables_1000 <- Reduce("+",otu_tables_1000)/length(otu_tables_1000)

### IMPORTANT! NEED TO ROUND IT OTHERWISE iNEXT WILL NOT WORK!
average_otu_tables_1000_round <- average_otu_tables_1000 %>% mutate_at(vars(starts_with("ASV")), funs(round(., 0)))

### add SampleID column back
average_otu_tables_1000_round$sample_id<- rownames(average_otu_tables_1000) 
average_otu_tables_1000_round <- average_otu_tables_1000_round %>% 
  select(sample_id, everything())

### save the averaged phyloseq objects 
write.csv(average_otu_tables_1000_round, "phylo16_average_otu_tables.csv", quote=F, row.names=F)

##### SAVE RAREFIED PHYLOSEQ OBJECT 
oturare = read.csv("phylo16_average_otu_tables.csv") |> 
  column_to_rownames(var="sample_id") |> 
  as.matrix()

rare = phyloseq(sample_data(phylo16@sam_data),
                otu_table(oturare, taxa_are_rows = F),
                tax_table(phylo16@tax_table))
rare

write_rds(rare, "MH_phylo16_coverage_rarefied_2021to2023.RDS")


##### 18S #####
##### look at minimum, mean, and maximum sample counts 
smin <- min(sample_sums(phylo18)) 
meanreads <-mean(sample_sums(phylo18)) 
smax <- max(sample_sums(phylo18))
totalreads <- sum(sample_sums(phylo18))

##### USE iNEXT TO FIND SAMPLE COVERAGE
### check if taxa are rows in the phyloseq object
taxa_are_rows(phylo18)
## mine is TRUE

### prepare otu table for iNEXT function
x <- metagMisc::prepare_inext(
  as.data.frame(otu_table(phylo18)),
  correct_singletons = T)

## calculate coverage with iNEXT Char.Ind
SC <- plyr::llply(.data = x, .fun = function(z){ try( iNEXT:::Chat.Ind(z, sum(z)) ) })
plyr::ldply(.data = SC, .fun = class)

## get dataframe to sort coverage and decide what to pick 
coverage = as.data.frame(t(as.matrix(as.data.frame(SC))))

#### run coveraged-based  rarefaction 
COVERAGE_RAREF_1000 <- phyloseq_coverage_raref(physeq=phylo18,
                                               coverage = 0.8,
                                               iter = 1000,
                                               replace = F, 
                                               correct_singletons = T,
                                               drop_lowcoverage = F,
                                               multithread = T) 


write_rds(COVERAGE_RAREF_1000, "COVERAGE_RAREF_MH_phylo18_phyloseq.RDS")
COVERAGE_RAREF_1000 = readRDS("COVERAGE_RAREF_MH_phylo18_phyloseq.RDS")


### select to to the entire 1000 phyloseq objects
all_1000_phyloseq <- COVERAGE_RAREF_1000[c(1:1000)]

### first, extract otu tables from the 1000 phyloseq objects
otu_tables_1000 <- lapply(all_1000_phyloseq, function(z) as.data.frame(t(phyloseq::otu_table(z))))

### average all matrices to get the mean abundance across all iterations
average_otu_tables_1000 <- Reduce("+",otu_tables_1000)/length(otu_tables_1000)

### IMPORTANT! NEED TO ROUND IT OTHERWISE iNEXT WILL NOT WORK!
average_otu_tables_1000_round <- average_otu_tables_1000 %>% mutate_at(vars(starts_with("ASV")), funs(round(., 0)))

### add SampleID column back
average_otu_tables_1000_round$sample_id<- rownames(average_otu_tables_1000) 
average_otu_tables_1000_round <- average_otu_tables_1000_round %>% 
  select(sample_id, everything())

### save the averaged phyloseq objects 
write.csv(average_otu_tables_1000_round, "phylo18_average_otu_tables.csv", quote=F, row.names=F)

##### SAVE RAREFIED PHYLOSEQ OBJECT 
oturare = read.csv("phylo18_average_otu_tables.csv") |> 
  column_to_rownames(var="sample_id") |> 
  as.matrix()

rare = phyloseq(sample_data(phylo18@sam_data),
                otu_table(oturare, taxa_are_rows = F),
                tax_table(phylo18@tax_table))
rare

write_rds(rare, "MH_phylo18_coverage_rarefied_2021to2023.RDS")

##### RBCL #####
##### look at minimum, mean, and maximum sample counts 
smin <- min(sample_sums(phylorb)) 
meanreads <-mean(sample_sums(phylorb)) 
smax <- max(sample_sums(phylorb))
totalreads <- sum(sample_sums(phylorb))

##### USE iNEXT TO FIND SAMPLE COVERAGE
### check if taxa are rows in the phyloseq object
taxa_are_rows(phylorb)
## mine is TRUE

### prepare otu table for iNEXT function
x <- metagMisc::prepare_inext(
  as.data.frame(otu_table(phylorb)),
  correct_singletons = T)

## calculate coverage with iNEXT Char.Ind
SC <- plyr::llply(.data = x, .fun = function(z){ try( iNEXT:::Chat.Ind(z, sum(z)) ) })
plyr::ldply(.data = SC, .fun = class)

## get dataframe to sort coverage and decide what to pick 
coverage = as.data.frame(t(as.matrix(as.data.frame(SC))))

#### run coveraged-based  rarefaction 
COVERAGE_RAREF_1000 <- phyloseq_coverage_raref(physeq=phylorb,
                                               coverage = 0.8,
                                               iter = 1000,
                                               replace = F, 
                                               correct_singletons = T,
                                               drop_lowcoverage = F,
                                               multithread = T) 


write_rds(COVERAGE_RAREF_1000, "COVERAGE_RAREF_MH_phylorb_phyloseq.RDS")

### select to to the entire 1000 phyloseq objects
COVERAGE_RAREF_1000=readRDS("COVERAGE_RAREF_MH_phylorb_phyloseq.RDS")
all_1000_phyloseq <- COVERAGE_RAREF_1000[c(1:1000)]

### first, extract otu tables from the 1000 phyloseq objects
otu_tables_1000 <- lapply(all_1000_phyloseq, function(z) as.data.frame(t(phyloseq::otu_table(z))))

### average all matrices to get the mean abundance across all iterations
average_otu_tables_1000 <- Reduce("+",otu_tables_1000)/length(otu_tables_1000)

### IMPORTANT! NEED TO ROUND IT OTHERWISE iNEXT WILL NOT WORK!
average_otu_tables_1000_round <- average_otu_tables_1000 %>% mutate_at(vars(starts_with("ASV")), funs(round(., 0)))

### add SampleID column back
average_otu_tables_1000_round$sample_id<- rownames(average_otu_tables_1000) 
average_otu_tables_1000_round <- average_otu_tables_1000_round %>% 
  select(sample_id, everything())

### save the averaged phyloseq objects 
write.csv(average_otu_tables_1000_round, "phylorb_average_otu_tables.csv", quote=F, row.names=F)

##### SAVE RAREFIED PHYLOSEQ OBJECT 
oturare = read.csv("phylorb_average_otu_tables.csv") |> 
  column_to_rownames(var="sample_id") |> 
  as.matrix()

rare = phyloseq(sample_data(phylorb@sam_data),
                otu_table(oturare, taxa_are_rows = F),
                tax_table(phylorb@tax_table))
rare

write_rds(rare, "MH_phyloRBCL_coverage_rarefied_2021to2023.RDS")

## end