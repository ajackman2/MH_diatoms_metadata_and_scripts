###### SET UP #######
## load packages
library(phyloseq)
library(tidyverse)
library(plyr)
library(qualpalr)
library(ggpubr)
library(ggh4x)
library(indicspecies)
library(data.table)
library(ggplot2); theme_set(theme_bw()+
                              theme(panel.grid = element_blank(),
                                    strip.background = element_rect(fill="white"),
                                    axis.text.y = element_text(size = 10, colour = "black"),
                                    axis.title = element_text(size=10, face="bold"),
                                    strip.text = element_text(color="black", size=10),
                                    legend.text=element_text(size=10),
                                    axis.line = element_line(colour = "black"),
                                    axis.text.x = element_blank()))

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

setwd("C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Seagrass_Montague_2021-03-07/galliano_seagrass_2021and2023/output")

t18 = readRDS("C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Seagrass_Montague_2021-03-07/galliano_seagrass_2021and2023/imput_data/MH_18S_filtered001percent_notrarefied_phyloseq.rds")
rbcl = readRDS("C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Seagrass_Montague_2021-03-07/galliano_seagrass_2021and2023/imput_data/MH_RBCL_filtered001percent_notrarefied_phyloseq.rds")

t18@sam_data$rd_filt = sample_sums(t18)
rbcl@sam_data$rd_filt = sample_sums(rbcl)


##### only keep diatoms #####
rbcl.phyla.list = as.data.frame(as.matrix(rbcl@tax_table))
rbcl.phyla.list = c(unique(rbcl.phyla.list$Phylum))
rbcl = subset_taxa(rbcl, Phylum !="Ochrophyta")

t18.phyla.list = as.data.frame(as.matrix(t18@tax_table))
t18.phyla.list = c(unique(t18.phyla.list$phylum))
print(t18.phyla.list)
t18.phyla.list = c("Ochrophyta", "Eukaryota")
t18 = subset_taxa(t18, phylum %in% c(t18.phyla.list))
#View(t18@tax_table)
t18.order.list = as.data.frame(as.matrix(t18@tax_table))
t18.order.list = c(unique(t18.order.list$order))
print(t18.order.list)
t18.order.list= c("Bacillariophytina", "Coscinodiscophytina","Eukaryota","Ochrophyta","Diatomea")
t18 = subset_taxa(t18, order %in% c(t18.order.list))

##### get taxa tables as data frames #####
rbcl.df = as.data.frame(as.matrix(rbcl@tax_table))
t18.df = as.data.frame(as.matrix(t18@tax_table))


write.csv(rbcl.df, "diatoms_in_rbcl_data.csv")
write.csv(t18.df, "diatoms_in_18S_data.csv")
getwd()
