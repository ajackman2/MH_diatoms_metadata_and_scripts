# incorporating the BLAST results into our data

# load libraries
library(dplyr)
library(plyr)
library(tidyverse)
library(tidyr)
library(stringr)
library(phyloseq)
library(qualpalr)
library(ggpubr)
library(ggh4x)
library(ggplot2); theme_set(theme_bw()+
                              theme(panel.grid = element_blank(),
                                    strip.background = element_rect(fill="white"),
                                    axis.text.y = element_text(size = 12, colour = "black"),
                                    axis.title = element_text(size=15, face="bold"),
                                    strip.text = element_text(color="black", size=10),
                                    legend.text=element_text(size=10),
                                    axis.line = element_line(colour = "black")))


## load functions
dephyloseq = function(phylo_obj){
  
  ## get the metadata
  meta = as.data.frame(as.matrix(phylo_obj@sam_data))
  
  ## how many metadata columns you have 
  metacols = ncol(meta)+1
  
  ## get out the otu table 
  ## if your metadta is empty after running this, you need to use 
  otu = as.data.frame(t(as.matrix(phylo_obj@otu_table)))
  #otu = as.data.frame(as.matrix(phylo_obj@otu_table))
  
  ## merge the metadata and otu table by the rownames (sample ids from the Illumina sequencing data)
  mo = merge(meta, otu, by=0)
  
  ## get out the taxonomy file 
  tax = as.data.frame(phylo_obj@tax_table)
  
  ## get the ASV ID out. This the matches the placeholder ASV ID in the OTU table
  tax = tax %>% rownames_to_column(var="ASVid")
  
  ## pivot longer to be able to match the ASVs in the OTU table to the taxonomy table 
  mo = mo %>% pivot_longer(cols = -c(1:metacols), names_to = "ASVid", values_to="asv_abundance")
  
  ## Join the metadata and otu table with the taoxnomy table 
  mot = full_join(mo, tax)
  
  ## Specify the output for the dephyloseq funciton 
  output = mot
}

setwd("/Users/andreajackman/R_Stuff/Parfrey/MH_diatoms_metadata_and_scripts")

# read in the data
diatom_data <- read.csv("morphological_data/master_diatom_data.csv")
mol_18S = readRDS("molecular_data/MH_18S_phyloseq.rds")
mol_RBCL = readRDS("molecular_data/MH_RBCL_phyloseq.rds")

##### FORMAT SEM COUNTS ######
diatom_data = diatom_data |>  pivot_longer(cols=c(9:72),
                                           names_to = "Genus",
                                           values_to = "diatom_count")

diatom_data = ddply(diatom_data, c("Genus"),
                    summarise,
                    sum.sem.counts = sum(diatom_count),
                    n.sem.samples = length(unique(box)))

diatom_data$mean.sem.count = diatom_data$sum.sem.counts/diatom_data$n.sem.samples

diatom_data$in_counts = ifelse(diatom_data$mean.sem.count, "yes", "excluded")

# fixing one entry
diatom_data <- diatom_data %>%
  mutate(Genus = ifelse(Genus == "Paribellus", "Parlibellus", Genus))

##### FORMAT MOLECULAR DATA #####
## per sample read depth
mol_RBCL@sam_data$rd_filt = sample_sums(mol_RBCL)
mol_18S@sam_data$rd_filt = sample_sums(mol_18S)

## remove non-diatoms
mol_RBCL = subset_taxa(mol_RBCL, Phylum !="Eukaryota" &
                         Kingdom!="Sar") #Heterosigma akashiwo chloroplast DNA by BLAST

mol_18S = subset_taxa(mol_18S, phylum %in% c("Ochrophyta") &
                        !(class %in% c("Chrysophyceae", "Ochrophyta", "Phaeophyceae", "Pelagophyceae", "Raphidophyceae")))

## get out of phyloseq
mol_RBCL = dephyloseq(mol_RBCL)
mol_18S = dephyloseq(mol_18S) 

# read in the files with the new phylogenetic tree and blast results
blast_18S <- read.csv("output/diatoms_in_18S_data.csv")
blast_rbcl <- read.csv("output/diatoms_in_rbcl_data_2.csv")

#### RBCL data ####
# reassign the genera
merged_data <- merge(mol_RBCL,blast_rbcl, by = "asv_id", suffixes = c("_mol", "_blast"))
merged_data$Genus_mol <- merged_data$tree_corrected
merged_data$Order_mol <- merged_data$tree_order
merged_data <- merged_data |>
  select(asv_id:subsp_mol)

# edit mol_RBCL
for (i in 1:nrow(mol_RBCL)){
  for (j in 1:nrow(merged_data)){
  if (mol_RBCL$asv_id[i] == merged_data$asv_id[j]) {
      mol_RBCL$Genus[i] <- merged_data$Genus_mol[j]
      mol_RBCL$Order[i] <- merged_data$Order_mol[j]
    }
  }
}

#### 18S data ####
# reassign the genera
merged_data <- merge(mol_18S, blast_18S, by = "asv_id", suffixes = c("_mol", "_blast"))
merged_data$genus_mol <- merged_data$tree_corrected
merged_data$order_mol <- merged_data$tree_order
merged_data <- merged_data |>
  select(asv_id:genus_mol)


# edit mol_18S
for (i in 1:nrow(mol_18S)){
  for (j in 1:nrow(merged_data)){
    if (mol_18S$asv_id[i] == merged_data$asv_id[j]) {
      mol_18S$genus[i] <- merged_data$genus_mol[j]
      mol_18S$order[i] <- merged_data$order_mol[j]
    }
  }
}

# save the new dataframes
write.csv(mol_18S, "counting_data/tree_corrected_mol_18S")
write.csv(mol_RBCL, "counting_data/tree_corrected_mol_RBCL")
