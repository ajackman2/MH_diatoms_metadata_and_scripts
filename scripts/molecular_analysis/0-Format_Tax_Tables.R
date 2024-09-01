###### SET UP ######
library(tidyverse)
library(phyloseq)

setwd("C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Seagrass_Montague_2021-03-07/galliano_seagrass_2021and2023/galliano_tax_tables")

diatoms = read.csv("diatoms_DiatBarcodeV10.csv")
euk = read.csv("eukaryotes_Silva132.csv")

galbb = readRDS("C:/Users/siobh/OneDrive - The University Of British Columbia/Project - LP MW AS RBCL Galliano Seq/laura_galiano_2023bb_filtered.RDS")
mh18 = readRDS("C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Seagrass_Montague_2021-03-07/galliano_seagrass_2021and2023/imput_data/MH_18S_filtered001percent_notrarefied_phyloseq.rds")
mhrb = readRDS("C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Seagrass_Montague_2021-03-07/galliano_seagrass_2021and2023/imput_data/MH_RBCL_filtered001percent_notrarefied_phyloseq.rds")
mwdiatoms =readRDS("C:/Users/siobh/OneDrive - The University Of British Columbia/Project - LP MW AS RBCL Galliano Seq/mark_diatoms_2023_filtered.RDS")
mw18 = readRDS("C:/Users/siobh/OneDrive - The University Of British Columbia/Project - LP MW AS RBCL Galliano Seq/mark_18S_2023_filtered.RDS")

##### MAKE LIST ######
## galliano bio blitz
galbb = as.data.frame(galbb@tax_table)
galbb = as.data.frame(c(unique(galbb$Species)))
names(galbb) <- c("Species")
galbb$InBioBlitz <- "y"

## seagrass diatoms 
mhrb = as.data.frame(mhrb@tax_table)
mhrb = as.data.frame(c(unique(mhrb$Species)))
names(mhrb) <- c("Species")
mhrb$InSeagrass <-"y"

## Mark diatoms
mwdiatoms = as.data.frame(mwdiatoms@tax_table)
mwdiatoms = as.data.frame(c(unique(mwdiatoms$Species)))
names(mwdiatoms) <-c("Species")
mwdiatoms$InMark <-"y"

## seagrass 18S
mh18 = as.data.frame(mh18@tax_table)
mh18 = as.data.frame(c(unique(mh18$species)))
names(mh18) <-c("Species")
mh18$InSeagrass <-"y"

## Mark 18S
mw18 = as.data.frame(mw18@tax_table)
mw18 = as.data.frame(c(unique(mw18$Species)))
names(mw18) <-c("Species")
mw18$InMark <-"y"

###### ADD VARIABLE INFO - DIATOMS ######
## merge diatoms datasets
detects = merge(merge(galbb, mhrb, all=T), mwdiatoms, all=T)
detects[is.na(detects)]<-"n"

withloc = full_join(diatoms, detects)
write.csv(withloc, 'galliano_diatoms.csv')

##### MERGE 18S TOGETHER #####
detects = merge(mh18, mw18, all=T)
detects[is.na(detects)]<-"n"

withloc = full_join(euk, detects)
write.csv(withloc, 'galliano_eukaryotes.csv')


