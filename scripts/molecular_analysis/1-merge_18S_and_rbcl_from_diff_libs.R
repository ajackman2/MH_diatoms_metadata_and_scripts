#### SET UP #####
## load packages
library(plyr)
library(tidyverse)
library(phyloseq)
library(dada2)
library(data.table)
library(zoo)


## set workign directory 
setwd("/parfreylab/schenk/seagrass_montague_processing")

## read in files
st18.1 = as.matrix(readRDS("seqtab_nochim_18S_20231027_fw_only.RDS"))
st18.2 = as.matrix(readRDS("seqtab_nochim_18S_MH_2021.RDS")) |>  
  t() |> 
  as.data.frame() |> 
  rownames_to_column(var="asv_seq")

strb.1 = as.matrix(readRDS("seqtab_nochim_LP_AS_MW_MHSGredo.rds"))
strb.2 = as.matrix(readRDS("seqtab_nochim_RBCL_MH_2021.RDS"))

#### trim 18S to be the same lengths ####
# 200 bp "C:\Users\siobh\OneDrive - The University Of British Columbia\Project - LP MW AS RBCL Galliano Seq\1-18S-AS_MW_redo_initial_processing_fwonly.R"

# trim to 175 bp to match 2022 redo fw only 
st18.2$asv_seq = substr(st18.2$asv_seq, 1, 200)

## 
otutab = ncol(st18.2)

# turn the seqtab 
st18.2long = pivot_longer(st18.2, 
                             cols=c(2:otutab),
                             names_to = "sample_id",
                             values_to = "asv_count")

# summarize to merge the sequences that ar the same
st18.2long = ddply(st18.2long, c("asv_seq","sample_id"),
                          summarise,
                          asv.n = sum(asv_count))

# pivot back to wide format and format for analysis
st18.2form = pivot_wider(st18.2long,
                            names_from = "asv_seq",
                            values_from = "asv.n") |> 
  column_to_rownames(var="sample_id") |> 
  as.matrix()

##### MERGE 18S SEQTBAS TOGETHER #####
st18 <- mergeSequenceTables(st18.1, st18.2form, tryRC=TRUE)
## remove bimeras
seqtab.nochim18 <- removeBimeraDenovo(st18, method="consensus", multithread=TRUE,verbose = FALSE)
dim(seqtab.nochim18)

# Collapse ASVs that only differ by end base pair
seqtab.nochim18<-collapseNoMismatch(seqtab.nochim18, verbose=T) # Check Issue 716 from dada2 github
# https://github.com/benjjneb/dada2/issues/626
dim(seqtab.nochim18)

## save seqtab file as RDS
write_rds(seqtab.nochim18, "MH_2021to2023_18S_seqtab_nochim.RDS")

##### MERGE RBCL SEQTBAS TOGETHER #####
strb <- mergeSequenceTables(strb.1, strb.2, tryRC=TRUE)
## remove bimeras
seqtab.nochimrb <- removeBimeraDenovo(strb, method="consensus", multithread=TRUE,verbose = FALSE)
dim(seqtab.nochimrb)

# Collapse ASVs that only differ by end base pair
seqtab.nochimrb<-collapseNoMismatch(seqtab.nochimrb, verbose=T) # Check Issue 716 from dada2 github
# https://github.com/benjjneb/dada2/issues/626
dim(seqtab.nochimrb)

## save seqtab file as RDS
write_rds(seqtab.nochimrb, "MHandGalliano_2021to2023_RBCL_seqtab_nochim.RDS")

##### SILVA 132 #####

## SILVA needs to be downloaded in the directory to run
taxa <- assignTaxonomy(seqtab.nochim18, "/parfreylab/schenk/taxonomy_databases/Silva138/silva_nr_v132_train_set.fa", 
                       multithread=TRUE, tryRC=TRUE)

# Save taxonomy table
write.table(data.frame("row_names"=rownames(taxa),taxa),"MH2021to2023_tax_SILVAv132_18s.txt", 
            row.names=FALSE, quote=F, sep="\t")

## add species information to taxonomy data
taxa_sp <- addSpecies(taxa, "/parfreylab/schenk/taxonomy_databases/Silva138/silva_species_assignment_v138.fa")

#NA taxa are hard to separate later if they have no label. apply "Unassigned" label here now.
unique(taxa_sp[,1]) #possible labels here: eukaryotic, archaeal, bacterial, and "NA" taxa. 
NAs <- is.na(taxa_sp[,1]) #test for NA
NAs <- which(NAs == TRUE) #get indices of NA values
taxa_sp[NAs,1] <- "Unassigned" #apply new label to identified indices
colnames(taxa_sp) <- c("domain", "phylum", "class", "order", "family", "genus", "species") 


# Save taxonomy table
write.table(data.frame("row_names"=rownames(taxa_sp),taxa_sp),"MH2021to2023_tax_noNA_SILVAv132_18s.txt", 
            row.names=FALSE, quote=F, sep="\t")

##### DIATBARCODE #####
## SILVA needs to be downloaded in the directory to run
taxa.rb <- assignTaxonomy(seqtab.nochimrb, "/parfreylab/schenk/taxonomy_databases/DiatBarcodeV10/diat_barcode_v10_tax_assign_dada2.fa", 
                       multithread=TRUE)

# Save taxonomy table
write.table(data.frame("row_names"=rownames(taxa.rb),taxa.rb),"MHandGalliano_taxonomy_diatv10_2021to2023.txt", 
            row.names=FALSE, quote=F, sep="\t")

## add species information to taxonomy data
taxa_sp_rb <- addSpecies(taxa.rb, "/parfreylab/schenk/taxonomy_databases/DiatBarcodeV10/diat_barcode_v10_sp_assign_dada2.fa")

#NA taxa are hard to separate later if they have no label. apply "Unassigned" label here now.
unique(taxa_sp_rb[,1]) #possible labels here: eukaryotic, archaeal, bacterial, and "NA" taxa. 
NAs <- is.na(taxa_sp_rb[,1]) #test for NA
NAs <- which(NAs == TRUE) #get indices of NA values
taxa_sp[NAs,1] <- "Unassigned" #apply new label to identified indices
colnames(taxa_sp_rb) <- c("Domain","Kingdom","Subkingdom","Phylum","Class","Order","Empty","Genus","Species","subsp") 



# Save taxonomy table
write.table(data.frame("row_names"=rownames(taxa_sp_rb),taxa_sp_rb),"MHandGalliano_taxonomy_noNA_diatv10_2021to2023.txt", 
            row.names=FALSE, quote=F, sep="\t")

# mank 01 2023-11-13 eneving 
## back to pc