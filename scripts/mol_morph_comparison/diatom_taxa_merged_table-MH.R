##### set up #####
# taxa merging
# MJ Herrin (from A. Jackman scripts)
# Dec 1 2023

# libraries
library(dplyr)
library(tidyverse)
library(plyr)
library(tidyr)
library(stringr)
library(phyloseq)
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


# read in the data
diatom_data <- read.csv("counting_data/master_diatom_data.csv")
mol_18S = readRDS("imput_data/MH_18S_filtered001percent_notrarefied_phyloseq.rds")
mol_RBCL = readRDS("imput_data/MH_RBCL_filtered001percent_notrarefied_phyloseq.RDS")

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

## get per sample relative abundance
mol_RBCL$ra = as.numeric(mol_RBCL$asv_abundance)/as.numeric(mol_RBCL$rd_filt)
mol_18S$ra = as.numeric(mol_18S$asv_abundance)/as.numeric(mol_18S$rd_filt)

## rename
names(mol_18S)[names(mol_18S)=="genus"]<-"Genus"

inRBCL = c(unique(mol_RBCL$Genus))
in18S = c(unique(mol_18S$Genus))

mol = full_join(mol_RBCL, mol_18S)

mol = ddply(mol, c("Genus"),
            summarise,
            sum.mol.ra = sum(ra),
            mean.mol.ra = mean(ra),
            n.mol.samples =length(unique(Row.names)))

mol$in18S = ifelse(mol$Genus %in% c(in18S), "yes", "no")
mol$inRBCL = ifelse(mol$Genus %in% c(inRBCL), "yes", "no")


###### combine data #####
all = full_join(mol, diatom_data)


write.csv(all, "counting_data/combined_sem_illumina_counts_new.csv")
