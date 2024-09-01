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


# read in the data
diatom_data <- read.csv("counting_data/master_diatom_data.csv")
mol_18S = readRDS("imput_data/MH_18S_filtered001percent_notrarefied_phyloseq.rds")
mol_RBCL = readRDS("imput_data/MHandGalliano_RBCL_filtered_phyloseq.RDS")

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

# re-assign ASVs - only ones with confidence from Mark
# extract Achnanthales

only_achnan <- mol_RBCL |>
  filter(Species == "Achnanthales")

# write a fasta file
fasta_entries <- character(nrow(only_achnan))

# Loop through each row of the dataframe
for (i in 1:nrow(only_achnan)) {
  # Extract ASVid and Domain:Species information
  asv_id <- only_achnan$ASVid[i]
  domain_species <- paste(only_achnan[i, "Domain"], only_achnan[i, "Phylum"], only_achnan[i, "Class"], 
                          only_achnan[i, "Order"], only_achnan[i, "Genus"], only_achnan[i, "Species"], sep = ":")
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

fasta_entries <- unique(fasta_entries)
# Write the FASTA entries to a file
writeLines(fasta_entries, "blast/achnanthales_diat.fasta")

mol_RBCL$Genus[mol_RBCL$asv_id == 'asv237'] <- "Planothidium"
mol_RBCL$Genus[mol_RBCL$asv_id == 'asv482'] <- "Achnanthidium"

# extract Bacillariales

only_bacill <- mol_RBCL |>
  filter(Species == "Bacillariales")

# write a fasta file
fasta_entries <- character(nrow(only_bacill))

# Loop through each row of the dataframe
for (i in 1:nrow(only_bacill)) {
  # Extract ASVid and Domain:Species information
  asv_id <- only_bacill$ASVid[i]
  domain_species <- paste(only_bacill[i, "Domain"], only_bacill[i, "Phylum"], only_bacill[i, "Class"], 
                          only_bacill[i, "Order"], only_bacill[i, "Genus"], only_bacill[i, "Species"], sep = ":")
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

fasta_entries <- unique(fasta_entries)
# Write the FASTA entries to a file
writeLines(fasta_entries, "blast/bacillariales_diat.fasta")

mol_RBCL$Genus[mol_RBCL$asv_id == 'asv631'] <- "Petroneis" # 0 for petroneis and 10 mismatches for cylindrotheca
# mol_RBCL$Genus[mol_RBCL$asv_id == 'asv510'] <- "Cylindrotheca" # 12 mismatches
mol_RBCL$Genus[mol_RBCL$asv_id == 'asv154'] <- "Nitzschia" # 0 mismatches
mol_RBCL$Genus[mol_RBCL$asv_id == 'asv234'] <- "Cylindrotheca" # 0 mismatches

# extract Bacillariophyceae

only_bacillario <- mol_18S |>
  filter(genus == "Bacillariophyceae")

# write a fasta file
fasta_entries <- character(nrow(only_bacillario))

# Loop through each row of the dataframe
for (i in 1:nrow(only_bacillario)) {
  # Extract ASVid and Domain:Species information
  asv_id <- only_bacillario$ASVid[i]
  domain_species <- paste(only_bacillario[i, "domain"], only_bacillario[i, "phylum"], only_bacillario[i, "class"], 
                          only_bacillario[i, "order"], only_bacillario[i, "genus"],  sep = ":")
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

fasta_entries <- unique(fasta_entries)
# Write the FASTA entries to a file
writeLines(fasta_entries, "blast/bacillariophyceae_silva.fasta")

mol_18S$genus[mol_18S$asv_id == 'asv178'] <- "Proschkinia" # 0 mismatches
mol_18S$genus[mol_18S$asv_id == 'asv204'] <- "Achnanthidium" # 0 mismatches
mol_18S$genus[mol_18S$asv_id == 'asv557'] <- "Cocconeis" # 4 mismatches
mol_18S$genus[mol_18S$asv_id == 'asv793'] <- "Cocconeis" # 4 mismatches
mol_18S$genus[mol_18S$asv_id == 'asv998'] <- "Cocconeis" # 0 mismatches


# bacillariop
only_bacill <- mol_RBCL |>
  filter(Species == "Bacillariophyceae")

# write a fasta file
fasta_entries <- character(nrow(only_bacill))

# Loop through each row of the dataframe
for (i in 1:nrow(only_bacill)) {
  # Extract ASVid and Domain:Species information
  asv_id <- only_bacill$ASVid[i]
  domain_species <- paste(only_bacill[i, "Domain"], only_bacill[i, "Phylum"], only_bacill[i, "Class"], 
                          only_bacill[i, "Order"], only_bacill[i, "Genus"], only_bacill[i, "Species"], sep = ":")
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

fasta_entries <- unique(fasta_entries)
# Write the FASTA entries to a file
writeLines(fasta_entries, "blast/bacillariophyceae_diat.fasta")

# bacillariophyta
only_bacill <- mol_RBCL |>
  filter(Species == "Bacillariophyta")

# write a fasta file
fasta_entries <- character(nrow(only_bacill))

# Loop through each row of the dataframe
for (i in 1:nrow(only_bacill)) {
  # Extract ASVid and Domain:Species information
  asv_id <- only_bacill$ASVid[i]
  domain_species <- paste(only_bacill[i, "Domain"], only_bacill[i, "Phylum"], only_bacill[i, "Class"], 
                          only_bacill[i, "Order"], only_bacill[i, "Genus"], only_bacill[i, "Species"], sep = ":")
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

fasta_entries <- unique(fasta_entries)
# Write the FASTA entries to a file
writeLines(fasta_entries, "blast/Bacillariophyta_diat.fasta")

mol_RBCL$Genus[mol_RBCL$asv_id == 'asv101'] <- "Podosira" # 0 mismatches

# diatomea
only_diatomea <- mol_18S |>
  filter(species == "Diatomea")

# write a fasta file
fasta_entries <- character(nrow(only_diatomea))

# Loop through each row of the dataframe
for (i in 1:nrow(only_diatomea)) {
  # Extract ASVid and Domain:Species information
  asv_id <- only_diatomea$ASVid[i]
  domain_species <- paste(only_diatomea[i, "domain"], only_diatomea[i, "phylum"], only_diatomea[i, "class"], 
                          only_diatomea[i, "order"], only_diatomea[i, "genus"],  sep = ":")
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

fasta_entries <- unique(fasta_entries)
# Write the FASTA entries to a file
writeLines(fasta_entries, "blast/diatomea_silva.fasta")

mol_18S$genus[mol_18S$asv_id == 'asv182'] <- "Melosira" # 6 mismatches

# Fragilariales
only_frag <- mol_RBCL |>
  filter(Species == "Fragilariales")

# write a fasta file
fasta_entries <- character(nrow(only_frag))

# Loop through each row of the dataframe
for (i in 1:nrow(only_frag)) {
  # Extract ASVid and Domain:Species information
  asv_id <- only_frag$ASVid[i]
  domain_species <- paste(only_frag[i, "Domain"], only_frag[i, "Phylum"], only_frag[i, "Class"], 
                          only_frag[i, "Order"], only_frag[i, "Genus"], only_frag[i, "Species"], sep = ":")
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

fasta_entries <- unique(fasta_entries)
# Write the FASTA entries to a file
writeLines(fasta_entries, "blast/fragiliarles_diat.fasta")

mol_RBCL$Genus[mol_RBCL$asv_id == 'asv604'] <- "Gedaniella" # 0 mismatches

#### Fragiliarles ####
# diatomea
only_frag <- mol_18S |>
  filter(species == "Fragilariales")

# write a fasta file
fasta_entries <- character(nrow(only_frag))

# Loop through each row of the dataframe
for (i in 1:nrow(only_frag)) {
  # Extract ASVid and Domain:Species information
  asv_id <- only_frag$ASVid[i]
  domain_species <- paste(only_frag[i, "domain"], only_frag[i, "phylum"], only_frag[i, "class"], 
                          only_frag[i, "order"], only_frag[i, "genus"],  sep = ":")
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

fasta_entries <- unique(fasta_entries)
# Write the FASTA entries to a file
writeLines(fasta_entries, "blast/frag_silva.fasta")

mol_18S$genus[mol_18S$asv_id == 'asv182'] <- "Melosira" # 6 mismatches

#### Fragilariophyceae #####
# Fragilariales
only_frag <- mol_RBCL |>
  filter(Species == "Fragilariophyceae")

# write a fasta file
fasta_entries <- character(nrow(only_frag))

# Loop through each row of the dataframe
for (i in 1:nrow(only_frag)) {
  # Extract ASVid and Domain:Species information
  asv_id <- only_frag$ASVid[i]
  domain_species <- paste(only_frag[i, "Domain"], only_frag[i, "Phylum"], only_frag[i, "Class"], 
                          only_frag[i, "Order"], only_frag[i, "Genus"], only_frag[i, "Species"], sep = ":")
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

fasta_entries <- unique(fasta_entries)
# Write the FASTA entries to a file
writeLines(fasta_entries, "blast/fragilariophyceae_diat.fasta")

mol_RBCL$Genus[mol_RBCL$asv_id == 'asv49'] <- "Licmorphora" # 0 mismatches

#### Mediophyceae ####

only_medio <- mol_RBCL |>
  filter(Species == "Mediophyceae")

# write a fasta file
fasta_entries <- character(nrow(only_medio))

# Loop through each row of the dataframe
for (i in 1:nrow(only_medio)) {
  # Extract ASVid and Domain:Species information
  asv_id <- only_medio$ASVid[i]
  domain_species <- paste(only_medio[i, "Domain"], only_medio[i, "Phylum"], only_medio[i, "Class"], 
                          only_medio[i, "Order"], only_medio[i, "Genus"], only_medio[i, "Species"], sep = ":")
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

fasta_entries <- unique(fasta_entries)
# Write the FASTA entries to a file
writeLines(fasta_entries, "blast/mediophyceae_diat.fasta")

mol_RBCL$Genus[mol_RBCL$asv_id == 'asv364'] <- "Attheya" # 0 mismatches

#### Mediophyceae #####
# diatomea
only_medio <- mol_18S |>
  filter(species == "Mediophyceae")

# write a fasta file
fasta_entries <- character(nrow(only_medio))

# Loop through each row of the dataframe
for (i in 1:nrow(only_medio)) {
  # Extract ASVid and Domain:Species information
  asv_id <- only_medio$ASVid[i]
  domain_species <- paste(only_medio[i, "domain"], only_medio[i, "phylum"], only_medio[i, "class"], 
                          only_medio[i, "order"], only_medio[i, "genus"],  sep = ":")
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

fasta_entries <- unique(fasta_entries)
# Write the FASTA entries to a file
writeLines(fasta_entries, "blast/medio_silva.fasta")

mol_18S$genus[mol_18S$asv_id == 'asv182'] <- "Melosira" # 6 mismatches

# write.csv(mol_18S, "counting_data/blast_corrected_mol_18S")
# write.csv(mol_RBCL, "counting_data/blast_corrected_mol_RBCL")