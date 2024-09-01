##### set up ####
library(phyloseq)
library(tidyverse)
library(zoo)

setwd("C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Seagrass_Montague_2021-03-07/galliano_seagrass_2021and2023/imput_data")


seq = readRDS("MH_2021to2023_18S_seqtab_nochim.RDS")
tax = read.table("MH2021to2023_tax_noNA_SILVAv132_18s.txt", header=T)


#### propagates taxonomy from left ####
tax_propagated = tax %>%
  t() %>% #transpose (moves taxonomy from column names to row names)
  na.locf() %>% #fill the NAs with the values from the cell to the left (higher taxonomic rank)
  t() |> # transpose back to have column names be taxonomy and row names be ASV
  as.data.frame() |> 
  rownames_to_column(var="asv_id") |> 
  column_to_rownames(var="row_names")

ntax = ncol(tax_propagated)

## reorder columns
tax_propagated = tax_propagated[,c(2:ntax,1)]

## add asv id
tax_propagated$asv_id = paste0("asv", tax_propagated$asv_id)

##### make into unfiltered phyloseq object ####
meta = as.data.frame(rownames(seq))
meta$sample_id = meta$`rownames(seq)`
meta = meta |> column_to_rownames(var="sample_id")

#### fix metadata to remvoe non MH samples #####

meta = separate(meta, col=`rownames(seq)`,
                into=c("leaf_number","leaf_section", "target_region", "illumina_number"),
                sep="-")

meta = subset(meta, meta$leaf_number %in% c("Sterivex", "1", "2", "3", "4", "5", 
                                            "6", "7", "8", "9", "10","11", "blank", "pcr_blank_S358","extraction_blank_18S_S102"))

#### make phyloseq #####
uf = phyloseq(sample_data(meta),
              otu_table(as.matrix(seq), taxa_are_rows = F),
              tax_table(as.matrix(tax_propagated)))
uf



##### remove off target taxa #####

## keep for PR2
chl18S = subset_taxa(uf, order=='Chloroplast')
View(chl18S@tax_table)

write_rds(chl18S, "chloroplast_18S.RDS")

## for 18S filtering
uf = subset_taxa(uf, domain !="Bacteria" & 
                   domain!="Archaea"&
                   domain !="Unassigned" & genus!="Zostera")

#View(uf@tax_table)

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
filt18 = phyloseq(sample_data(uf@sam_data),
                  tax_table(uf@tax_table),
                  otu_table(as.matrix(otu), taxa_are_rows = T))
filt18@sam_data$rd_filt = sample_sums(filt18)
View(filt18@sam_data)


filt18 = subset_samples(filt18, rd_filt>0)

write_rds(filt18, "MH_18S_filtered001percent_notrarefied_phyloseq.rds")

