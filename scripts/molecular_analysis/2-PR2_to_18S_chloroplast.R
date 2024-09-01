## load packages
library(phyloseq)
library(tidyverse)
library(plyr)
library(dada2)
library(data.table)
library(zoo)

##### CLUSTER -  ASSIGN TAXONOMY #####

setwd("/parfreylab/schenk/seagrass_montague_processing/")
chl18 = readRDS("chloroplast_18S.RDS")



## SILVA needs to be downloaded in the directory to run
seqtab = as.matrix(as.data.frame(as.matrix(chl18@otu_table)))

taxa <- assignTaxonomy(seqtab, "/parfreylab/schenk/seagrass_montague_processing/pr2_version_4.14.0_SSU_dada2.fasta", 
                       multithread=TRUE, tryRC=TRUE,
                       taxLevels = c("Kingdom","Supergroup","Division","Class","Order","Family","Genus","Species"))

# Save taxonomy table
write_rds(taxa, "MH_pr2_taxonomy_18S.RDS")

##### DESKTOP #####
setwd("C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Seagrass_Montague_2021-03-07/galliano_seagrass_2021and2023/imput_data")

chl18 = readRDS("chloroplast_18S.RDS")
tax = readRDS("MH_pr2_taxonomy_18S.RDS")
