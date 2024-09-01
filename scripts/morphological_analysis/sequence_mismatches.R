# counting primer mismatches
# Andrea Jackman

# first need the sequences from the taxa I want to look into
# primer mismatches for cocconeis and tabularia
# mismatch for the 8 morpho not in molecular plus 21 from Mark

# read in the data fram diat barcode
# load libraries
library(dplyr)
library(plyr)
library(tidyverse)
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
mol_RBCL = readRDS("imput_data/MH_RBCL_filtered001percent_notrarefied_phyloseq.rds")

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

## get per sample relative abundance
mol_RBCL$ra = as.numeric(mol_RBCL$asv_abundance)/as.numeric(mol_RBCL$rd_filt)
mol_18S$ra = as.numeric(mol_18S$asv_abundance)/as.numeric(mol_18S$rd_filt)

## rename
names(mol_18S)[names(mol_18S)=="genus"]<-"Genus"

# inRBCL = c(unique(mol_RBCL$Genus))
# in18S = c(unique(mol_18S$Genus))
# 
# mol = full_join(mol_RBCL, mol_18S)
# 
# mol = ddply(mol, c("Genus"),
#             summarise,
#             sum.mol.ra = sum(ra),
#             mean.mol.ra = mean(ra),
#             n.mol.samples =length(unique(Row.names)))
# 
# mol$in18S = ifelse(mol$Genus %in% c(in18S), "yes", "no")
# mol$inRBCL = ifelse(mol$Genus %in% c(inRBCL), "yes", "no")

# remove not found in march from morphological data
diatom_data <- diatom_data %>%
  filter(!(Genus %in% c('Undatella', 'Actinoptychus', 'Trachyneis', 'Rhopalodia',
                        'Pseudonitzschia', 'Parlibellus', 'Hobaniella', 'Hanzschia',
                        'Gyrosigma', 'Cyclotella')))

# find the top 5 in 18S and RBCL
library(seqinr)

# filter for the top 5 diatoms
primer_mismatch_18S <- mol_18S |>
  filter(Genus %in% c('Navicula', 'Cylindrotheca', 'Nitzschia', 'Pseudogonphonema', 'Minidiscus'))

primer_mismatch_RBCL <- mol_RBCL |>
  filter(Genus %in% c('Navicula', 'Cylindrotheca', 'Nitzschia', 'Pseudogonphonema', 'Minidiscus'))

# make fasta files of the sequences
fasta_entries <- character(nrow(primer_mismatch_18S))

# Loop through each row of the dataframe
for (i in 1:nrow(primer_mismatch_18S)) {
  # Extract ASVid and Domain:Species information
  asv_id <- primer_mismatch_18S$ASVid[i]
  domain_species <- paste(primer_mismatch_18S[i, "domain"], primer_mismatch_18S[i, "phylum"], primer_mismatch_18S[i, "class"], 
                          primer_mismatch_18S[i, "order"], primer_mismatch_18S[i, "Genus"], primer_mismatch_18S[i, "species"], sep = ":")
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

# Write the FASTA entries to a file
writeLines(fasta_entries, "imput_data/primer_mismatch_18S.fasta")

# make fasta files of the sequences
fasta_entries <- character(nrow(primer_mismatch_RBCL))

# Loop through each row of the dataframe
for (i in 1:nrow(primer_mismatch_RBCL)) {
  # Extract ASVid and Domain:Species information
  asv_id <- primer_mismatch_RBCL$ASVid[i]
  domain_species <- paste(primer_mismatch_RBCL[i, "Domain"], primer_mismatch_RBCL[i, "Phylum"], primer_mismatch_RBCL[i, "Class"], 
                          primer_mismatch_RBCL[i, "Order"], primer_mismatch_RBCL[i, "Genus"], primer_mismatch_RBCL[i, "Species"], sep = ":")
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

# Write the FASTA entries to a file
writeLines(fasta_entries, "imput_data/primer_mismatch_RBCL.fasta")

# Load required libraries
library(Biostrings)

# Function to read sequences from a FASTA file
read_fasta <- function(file_path) {
  fasta_data <- readDNAStringSet(file_path)
  return(fasta_data)
}

# Read sequences from primer_mismatch_RBCL.fasta
primer_sequences <- read_fasta("imput_data/primer_mismatch_RBCL.fasta")

# Read sequences from diat_barcode_v10_tax_assign_dada2.fa
master_sequences <- read_fasta("imput_data/diat_barcode_v10_tax_assign_dada2.fa")


# make new dataframes to check the matches
# Initialize an empty dataframe
RBCL_data <- data.frame(ASVid = character(), domain_species = character(), stringsAsFactors = FALSE)

# Loop through each row of primer_mismatch_RBCL dataframe
for (i in 1:nrow(primer_mismatch_RBCL)) {
  # Extract ASVid and Domain:Species information
  asv_id <- primer_mismatch_RBCL$ASVid[i]
  domain_species <- paste(primer_mismatch_RBCL[i, "Domain"], primer_mismatch_RBCL[i, "Phylum"], primer_mismatch_RBCL[i, "Class"], 
                          primer_mismatch_RBCL[i, "Order"], primer_mismatch_RBCL[i, "Genus"], primer_mismatch_RBCL[i, "Species"], sep = ":")
  
  # Add the extracted information to the result dataframe
  RBCL_data <- rbind(RBCL_data, data.frame(ASVid = asv_id, domain_species = domain_species))
}

# Remove the first row which is empty
RBCL_data <- RBCL_data[-1, ]

data_18S <- data.frame(ASVid = character(), domain_species = character(), stringsAsFactors = FALSE)

# Loop through each row of primer_mismatch_RBCL dataframe
for (i in 1:nrow(primer_mismatch_18S)) {
  # Extract ASVid and Domain:Species information
  asv_id <- primer_mismatch_18S$ASVid[i]
  domain_species <- paste(primer_mismatch_18S[i, "domain"], primer_mismatch_18S[i, "phylum"], primer_mismatch_18S[i, "class"], 
                          primer_mismatch_18S[i, "order"], primer_mismatch_18S[i, "Genus"], primer_mismatch_18S[i, "species"], sep = ":")
  
  # Add the extracted information to the result dataframe
  data_18S <- rbind(data_18S, data.frame(ASVid = asv_id, domain_species = domain_species))
}

# Remove the first row which is empty
data_18S <- data_18S[-1, ]

#### check matches ####
# read master data in
master_sequences <- read.csv("imput_data/v10_correspondence.csv")
# Function to check if a sequence from primer_mismatch_RBCL.fasta matches any sequence in diat_barcode_v10_tax_assign_dada2.fa
check_matches <- function(RBCL_data, master_seqs) {
  matches <- vmatchPattern(as.character(RBCL_data), as.character(master_seqs))
  return(length(matches) > 0)
}

results_df <- data.frame(
  Sequence_ID = character(),
  Match_Found = logical(),
  stringsAsFactors = FALSE
)

# Check matches for each sequence in primer_mismatch_RBCL.fasta
for (i in 1:length(RBCL_data$ASVid)) {
  match_found <- check_matches(RBCL_data[i], master_sequences$sequence)
  results_df <- rbind(results_df, data.frame(Sequence_ID = names(primer_sequences)[i], Match_Found = match_found))
}

View(results_df) 

# check for the 18S database
primer_sequences <- read_fasta("imput_data/primer_mismatch_18S.fasta")

# Read sequences from diat_barcode_v10_tax_assign_dada2.fa
master_sequences <- read_fasta("imput_data/diat_barcode_v10_tax_assign_dada2.fa")

# Function to check if a sequence from primer_mismatch_RBCL.fasta matches any sequence in diat_barcode_v10_tax_assign_dada2.fa
check_matches <- function(primer_seq, master_seqs) {
  matches <- vmatchPattern(as.character(primer_seq), master_seqs)
  return(length(matches) > 0)
}

results_df <- data.frame(
  Sequence_ID = character(),
  Match_Found = logical()
)

# Check matches for each sequence in primer_mismatch_RBCL.fasta
for (i in 1:length(primer_sequences)) {
  match_found <- check_matches(primer_sequences[i], master_sequences)
  results_df <- rbind(results_df, data.frame(Sequence_ID = names(primer_sequences)[i], Match_Found = match_found))
}

View(results_df) # all have matches

# Function to check if a sequence ASVid from RBCL_data matches any sequence in master_seqs
check_matches <- function(RBCL_data, master_seqs) {
  asvid <- as.character(RBCL_data$ASVid)
  match_found <- sapply(master_seqs$sequence, function(seq) grepl(asvid, seq))
  return(any(match_found))
}

# Initialize a dataframe to store the results
results_df <- data.frame(
  Sequence_ID = character(),
  Match_Found = logical(),
  stringsAsFactors = FALSE
)

# Check matches for each sequence in RBCL_data
for (i in 1:nrow(RBCL_data)) {
  match_found <- check_matches(RBCL_data[i, ], master_sequences)
  results_df <- rbind(results_df, data.frame(Sequence_ID = RBCL_data$ASVid[i], Match_Found = match_found))
}

# Remove the first row which is empty
results_df <- results_df[-1, ]

# Print the resulting dataframe
View(results_df)

#### filter for the diatoms not detected well by illumina #####
# cocconeis, tabularia, fragilariopsis, gomphonemopsis, rhoicosphenia
# make fasta files of the sequences

# filter for cocconeis
cocconeis_RBCL <- mol_RBCL |>
  filter(Genus == 'Cocconeis')

cocconeis_RBCL <- cocconeis_RBCL[!duplicated(cocconeis_RBCL$ASVid), ]

fasta_entries <- character(nrow(cocconeis_RBCL))

# Loop through each row of the dataframe
for (i in 1:nrow(cocconeis_RBCL)) {
  # Extract ASVid and Domain:Species information
  asv_id <- cocconeis_RBCL$ASVid[i]
  domain_species <- paste(cocconeis_RBCL[i, "Domain"], cocconeis_RBCL[i, "Phylum"], cocconeis_RBCL[i, "Class"], 
                          cocconeis_RBCL[i, "Order"], cocconeis_RBCL[i, "Genus"], cocconeis_RBCL[i, "Species"], sep = ":")
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

# Write the FASTA entries to a file
writeLines(fasta_entries, "primer_mismatches/cocconeis_RBCL.fasta")

# get just cocconeis diat barcode
diat_barcode <- read.csv("imput_data/v10_correspondence.csv")

only_cocconeis <- diat_barcode |>
  filter(genus == "Cocconeis")

unique(only_cocconeis$sequence)

# write a fasta file
fasta_entries <- character(nrow(only_cocconeis))

# Loop through each row of the dataframe
for (i in 1:nrow(only_cocconeis)) {
  # Extract ASVid and Domain:Species information
  asv_id <- only_cocconeis$sequence[i]
  domain_species <- paste(only_cocconeis[i, "domain"], only_cocconeis[i, "phylum"], only_cocconeis[i, "class"], 
                          only_cocconeis[i, "order"], only_cocconeis[i, "genus"], only_cocconeis[i, "species_compromise"], sep = ":")
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

# Write the FASTA entries to a file
writeLines(fasta_entries, "primer_mismatches/cocconeis_diat.fasta")

#### Filter for most common ones ####
# filter for cocconeis
navicula_mol <- mol_RBCL |>
  filter(Genus == 'Navicula')

# navicula_18S <- mol_18S |>
#   filter(Genus == 'Navicula')

#navicula_mol <- full_join(navicula_RBCL, navicula_18S)
navicula_mol <- navicula_mol[!duplicated(navicula_mol$ASVid), ]

fasta_entries <- character(nrow(navicula_mol))

# Loop through each row of the dataframe
for (i in 1:nrow(navicula_mol)) {
  # Extract ASVid and Domain:Species information
  asv_id <- navicula_mol$ASVid[i]
  domain_species <- paste(navicula_mol[i, "Domain"], navicula_mol[i, "Phylum"], navicula_mol[i, "Class"], 
                          navicula_mol[i, "Order"], navicula_mol[i, "Genus"], navicula_mol[i, "Species"], sep = ":")
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

# Write the FASTA entries to a file
writeLines(fasta_entries, "primer_mismatches/navicula_RBCL.fasta")

# get just cocconeis diat barcode
diat_barcode <- read.csv("imput_data/v10_correspondence.csv")

only_navicula <- diat_barcode |>
  filter(genus == "Navicula")

# write a fasta file
fasta_entries <- character(nrow(only_navicula))

# Loop through each row of the dataframe
for (i in 1:nrow(only_navicula)) {
  # Extract ASVid and Domain:Species information
  asv_id <- only_navicula$sequence[i]
  domain_species <- paste(only_navicula[i, "domain"], only_navicula[i, "phylum"], only_navicula[i, "class"], 
                          only_navicula[i, "order"], only_navicula[i, "genus"], only_navicula[i, "species_compromise"], sep = ":")
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

# Write the FASTA entries to a file
writeLines(fasta_entries, "primer_mismatches/navicula_diat.fasta")

#### Tabularia ####
# filter for cocconeis
tabularia_mol <- mol_RBCL |>
  filter(Genus == 'Tabularia')

# navicula_18S <- mol_18S |>
#   filter(Genus == 'Navicula')

#navicula_mol <- full_join(navicula_RBCL, navicula_18S)
tabularia_mol <- tabularia_mol[!duplicated(tabularia_mol$ASVid), ]

fasta_entries <- character(nrow(tabularia_mol))

# Loop through each row of the dataframe
for (i in 1:nrow(tabularia_mol)) {
  # Extract ASVid and Domain:Species information
  asv_id <- tabularia_mol$ASVid[i]
  domain_species <- paste(tabularia_mol[i, "Domain"], tabularia_mol[i, "Phylum"], tabularia_mol[i, "Class"], 
                          tabularia_mol[i, "Order"], tabularia_mol[i, "Genus"], tabularia_mol[i, "Species"], sep = ":")
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

# Write the FASTA entries to a file
writeLines(fasta_entries, "primer_mismatches/tabularia_RBCL.fasta")

# get just cocconeis diat barcode
diat_barcode <- read.csv("imput_data/v10_correspondence.csv")

only_tabularia <- diat_barcode |>
  filter(genus == "Tabularia")

# write a fasta file
fasta_entries <- character(nrow(only_tabularia))

# Loop through each row of the dataframe
for (i in 1:nrow(only_tabularia)) {
  # Extract ASVid and Domain:Species information
  asv_id <- only_tabularia$sequence[i]
  domain_species <- paste(only_tabularia[i, "domain"], only_tabularia[i, "phylum"], only_tabularia[i, "class"], 
                          only_tabularia[i, "order"], only_tabularia[i, "genus"], only_tabularia[i, "species_compromise"], sep = ":")
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

# Write the FASTA entries to a file
writeLines(fasta_entries, "primer_mismatches/tabularia_diat.fasta")




#### Gomphonemopsis ####
gomphonemopsis_mol <- mol_RBCL |>
  filter(Genus == 'Gomphonemopsis')

# navicula_18S <- mol_18S |>
#   filter(Genus == 'Navicula')

#navicula_mol <- full_join(navicula_RBCL, navicula_18S)
gomphonemopsis_mol <- gomphonemopsis_mol[!duplicated(gomphonemopsis_mol$ASVid), ]

fasta_entries <- character(nrow(gomphonemopsis_mol))

# Loop through each row of the dataframe
for (i in 1:nrow(gomphonemopsis_mol)) {
  # Extract ASVid and Domain:Species information
  asv_id <- gomphonemopsis_mol$ASVid[i]
  domain_species <- paste(gomphonemopsis_mol[i, "Domain"], gomphonemopsis_mol[i, "Phylum"], gomphonemopsis_mol[i, "Class"], 
                          gomphonemopsis_mol[i, "Order"], gomphonemopsis_mol[i, "Genus"], gomphonemopsis_mol[i, "Species"], sep = ":")
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

# Write the FASTA entries to a file
writeLines(fasta_entries, "primer_mismatches/gomphonemopsis_RBCL.fasta")

# get just cocconeis diat barcode
diat_barcode <- read.csv("imput_data/v10_correspondence.csv")

only_gomp <- diat_barcode |>
  filter(genus == "Gomphonemopsis")

# write a fasta file
fasta_entries <- character(nrow(only_gomp))

# Loop through each row of the dataframe
for (i in 1:nrow(only_gomp)) {
  # Extract ASVid and Domain:Species information
  asv_id <- only_gomp$sequence[i]
  domain_species <- paste(only_gomp[i, "domain"], only_gomp[i, "phylum"], only_gomp[i, "class"], 
                          only_gomp[i, "order"], only_gomp[i, "genus"], only_gomp[i, "species_compromise"], sep = ":")
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

# Write the FASTA entries to a file
writeLines(fasta_entries, "primer_mismatches/gomphonemopsis_diat.fasta")


#### Skeletonema ####
skeletonema_mol <- mol_RBCL |>
  filter(Genus == 'Skeletonema')

# navicula_18S <- mol_18S |>
#   filter(Genus == 'Navicula')

#navicula_mol <- full_join(navicula_RBCL, navicula_18S)
skeletonema_mol <- skeletonema_mol[!duplicated(skeletonema_mol$ASVid), ]

fasta_entries <- character(nrow(skeletonema_mol))

# Loop through each row of the dataframe
for (i in 1:nrow(skeletonema_mol)) {
  # Extract ASVid and Domain:Species information
  asv_id <- skeletonema_mol$ASVid[i]
  domain_species <- paste(skeletonema_mol[i, "Domain"], skeletonema_mol[i, "Phylum"], skeletonema_mol[i, "Class"], 
                          skeletonema_mol[i, "Order"], skeletonema_mol[i, "Genus"], skeletonema_mol[i, "Species"], sep = ":")
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

# Write the FASTA entries to a file
writeLines(fasta_entries, "primer_mismatches/skeletonema_RBCL.fasta")

# get just cocconeis diat barcode
diat_barcode <- read.csv("imput_data/v10_correspondence.csv")

only_skel <- diat_barcode |>
  filter(genus == "Skeletonema")

# write a fasta file
fasta_entries <- character(nrow(only_skel))

# Loop through each row of the dataframe
for (i in 1:nrow(only_skel)) {
  # Extract ASVid and Domain:Species information
  asv_id <- only_skel$sequence[i]
  domain_species <- paste(only_skel[i, "domain"], only_skel[i, "phylum"], only_skel[i, "class"], 
                          only_skel[i, "order"], only_skel[i, "genus"], only_skel[i, "species_compromise"], sep = ":")
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

# Write the FASTA entries to a file
writeLines(fasta_entries, "primer_mismatches/skeletonema_diat.fasta")


#### Planothidium ####
planothidium_mol <- mol_RBCL |>
  filter(Genus == 'Planothidium')

# navicula_18S <- mol_18S |>
#   filter(Genus == 'Navicula')

#navicula_mol <- full_join(navicula_RBCL, navicula_18S)
planothidium_mol <- planothidium_mol[!duplicated(planothidium_mol$ASVid), ]

fasta_entries <- character(nrow(planothidium_mol))

# Loop through each row of the dataframe
for (i in 1:nrow(planothidium_mol)) {
  # Extract ASVid and Domain:Species information
  asv_id <- planothidium_mol$ASVid[i]
  domain_species <- paste(planothidium_mol[i, "Domain"], planothidium_mol[i, "Phylum"], planothidium_mol[i, "Class"], 
                          planothidium_mol[i, "Order"], planothidium_mol[i, "Genus"], planothidium_mol[i, "Species"], sep = ":")
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

# Write the FASTA entries to a file
writeLines(fasta_entries, "primer_mismatches/planothidium_RBCL.fasta")

# get just cocconeis diat barcode
diat_barcode <- read.csv("imput_data/v10_correspondence.csv")

only_plan <- diat_barcode |>
  filter(genus == "Planothidium")

# write a fasta file
fasta_entries <- character(nrow(only_plan))

# Loop through each row of the dataframe
for (i in 1:nrow(only_plan)) {
  # Extract ASVid and Domain:Species information
  asv_id <- only_plan$sequence[i]
  domain_species <- paste(only_plan[i, "domain"], only_plan[i, "phylum"], only_plan[i, "class"], 
                          only_plan[i, "order"], only_plan[i, "genus"], only_plan[i, "species_compromise"], sep = ":")
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

# Write the FASTA entries to a file
writeLines(fasta_entries, "primer_mismatches/planothidium_diat.fasta")


#### Cylindrotheca ####
cylindrotheca_mol <- mol_RBCL |>
  filter(Genus == 'Cylindrotheca')

# navicula_18S <- mol_18S |>
#   filter(Genus == 'Navicula')

#navicula_mol <- full_join(navicula_RBCL, navicula_18S)
cylindrotheca_mol <- cylindrotheca_mol[!duplicated(cylindrotheca_mol$ASVid), ]

fasta_entries <- character(nrow(cylindrotheca_mol))

# Loop through each row of the dataframe
for (i in 1:nrow(cylindrotheca_mol)) {
  # Extract ASVid and Domain:Species information
  asv_id <- cylindrotheca_mol$ASVid[i]
  domain_species <- paste(cylindrotheca_mol[i, "Domain"], cylindrotheca_mol[i, "Phylum"], cylindrotheca_mol[i, "Class"], 
                          cylindrotheca_mol[i, "Order"], cylindrotheca_mol[i, "Genus"], cylindrotheca_mol[i, "Species"], sep = ":")
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

# Write the FASTA entries to a file
writeLines(fasta_entries, "primer_mismatches/cylindrotheca_RBCL.fasta")

# get just cocconeis diat barcode
diat_barcode <- read.csv("imput_data/v10_correspondence.csv")

only_cyln <- diat_barcode |>
  filter(genus == "Cylindrotheca")

# write a fasta file
fasta_entries <- character(nrow(only_cyln))

# Loop through each row of the dataframe
for (i in 1:nrow(only_cyln)) {
  # Extract ASVid and Domain:Species information
  asv_id <- only_cyln$sequence[i]
  domain_species <- paste(only_cyln[i, "domain"], only_cyln[i, "phylum"], only_cyln[i, "class"], 
                          only_cyln[i, "order"], only_cyln[i, "genus"], only_cyln[i, "species_compromise"], sep = ":")
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

# Write the FASTA entries to a file
writeLines(fasta_entries, "primer_mismatches/cylindrotheca_diat.fasta")

#### Nitzschia ####
nitzschia_mol <- mol_RBCL |>
  filter(Genus == 'Nitzschia')

# navicula_18S <- mol_18S |>
#   filter(Genus == 'Navicula')

#navicula_mol <- full_join(navicula_RBCL, navicula_18S)
nitzschia_mol <- nitzschia_mol[!duplicated(nitzschia_mol$ASVid), ]

fasta_entries <- character(nrow(nitzschia_mol))

# Loop through each row of the dataframe
for (i in 1:nrow(nitzschia_mol)) {
  # Extract ASVid and Domain:Species information
  asv_id <- nitzschia_mol$ASVid[i]
  domain_species <- paste(nitzschia_mol[i, "Domain"], nitzschia_mol[i, "Phylum"], nitzschia_mol[i, "Class"], 
                          nitzschia_mol[i, "Order"], nitzschia_mol[i, "Genus"], nitzschia_mol[i, "Species"], sep = ":")
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

# Write the FASTA entries to a file
writeLines(fasta_entries, "primer_mismatches/nitzschia_RBCL.fasta")

# get just cocconeis diat barcode
diat_barcode <- read.csv("imput_data/v10_correspondence.csv")

only_nit <- diat_barcode |>
  filter(genus == "Nitzschia")

# write a fasta file
fasta_entries <- character(nrow(only_nit))

# Loop through each row of the dataframe
for (i in 1:nrow(only_nit)) {
  # Extract ASVid and Domain:Species information
  asv_id <- only_nit$sequence[i]
  domain_species <- paste(only_nit[i, "domain"], only_nit[i, "phylum"], only_nit[i, "class"], 
                          only_nit[i, "order"], only_nit[i, "genus"], only_nit[i, "species_compromise"], sep = ":")
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

# Write the FASTA entries to a file
writeLines(fasta_entries, "primer_mismatches/nitzschia_diat.fasta")


#### Pseudogonphonema ####
pseudogonphonema_mol <- mol_RBCL |>
  filter(Genus == 'Pseudogomphonema')

# navicula_18S <- mol_18S |>
#   filter(Genus == 'Navicula')

#navicula_mol <- full_join(navicula_RBCL, navicula_18S)
pseudogonphonema_mol <- pseudogonphonema_mol[!duplicated(pseudogonphonema_mol$ASVid), ]

fasta_entries <- character(nrow(pseudogonphonema_mol))

# Loop through each row of the dataframe
for (i in 1:nrow(pseudogonphonema_mol)) {
  # Extract ASVid and Domain:Species information
  asv_id <- pseudogonphonema_mol$ASVid[i]
  domain_species <- paste(pseudogonphonema_mol[i, "Domain"], pseudogonphonema_mol[i, "Phylum"], pseudogonphonema_mol[i, "Class"], 
                          pseudogonphonema_mol[i, "Order"], pseudogonphonema_mol[i, "Genus"], pseudogonphonema_mol[i, "Species"], sep = ":")
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

# Write the FASTA entries to a file
writeLines(fasta_entries, "primer_mismatches/pseudogonphonema_RBCL.fasta")

# get just cocconeis diat barcode
diat_barcode <- read.csv("imput_data/v10_correspondence.csv")

only_pseudogon <- diat_barcode |>
  filter(genus == "Pseudogomphonema")

# write a fasta file
fasta_entries <- character(nrow(only_pseudogon))

# Loop through each row of the dataframe
for (i in 1:nrow(only_pseudogon)) {
  # Extract ASVid and Domain:Species information
  asv_id <- only_pseudogon$sequence[i]
  domain_species <- paste(only_pseudogon[i, "domain"], only_pseudogon[i, "phylum"], only_pseudogon[i, "class"], 
                          only_pseudogon[i, "order"], only_pseudogon[i, "genus"], only_pseudogon[i, "species_compromise"], sep = ":")
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

# Write the FASTA entries to a file
writeLines(fasta_entries, "primer_mismatches/pseudogonphonema_diat.fasta")


#### Minidiscus ####
minidiscus_mol <- mol_RBCL |>
  filter(Genus == 'Minidiscus')

# navicula_18S <- mol_18S |>
#   filter(Genus == 'Navicula')

#navicula_mol <- full_join(navicula_RBCL, navicula_18S)
minidiscus_mol <- minidiscus_mol[!duplicated(minidiscus_mol$ASVid), ]

fasta_entries <- character(nrow(minidiscus_mol))

# Loop through each row of the dataframe
for (i in 1:nrow(minidiscus_mol)) {
  # Extract ASVid and Domain:Species information
  asv_id <- minidiscus_mol$ASVid[i]
  domain_species <- paste(minidiscus_mol[i, "Domain"], minidiscus_mol[i, "Phylum"], minidiscus_mol[i, "Class"], 
                          minidiscus_mol[i, "Order"], minidiscus_mol[i, "Genus"], minidiscus_mol[i, "Species"], sep = ":")
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

# Write the FASTA entries to a file
writeLines(fasta_entries, "primer_mismatches/minidiscus_RBCL.fasta")

# get just cocconeis diat barcode
diat_barcode <- read.csv("imput_data/v10_correspondence.csv")

only_mini <- diat_barcode |>
  filter(genus == "Minidiscus")

# write a fasta file
fasta_entries <- character(nrow(only_mini))

# Loop through each row of the dataframe
for (i in 1:nrow(only_mini)) {
  # Extract ASVid and Domain:Species information
  asv_id <- only_mini$sequence[i]
  domain_species <- paste(only_mini[i, "domain"], only_mini[i, "phylum"], only_mini[i, "class"], 
                          only_mini[i, "order"], only_mini[i, "genus"], only_mini[i, "species_compromise"], sep = ":")
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

# Write the FASTA entries to a file
writeLines(fasta_entries, "primer_mismatches/minidiscus_diat.fasta")


#### Halamphora ####
halamphora_mol <- mol_RBCL |>
  filter(Genus == 'Halamphora')

# navicula_18S <- mol_18S |>
#   filter(Genus == 'Navicula')

#navicula_mol <- full_join(navicula_RBCL, navicula_18S)
halamphora_mol <- halamphora_mol[!duplicated(halamphora_mol$ASVid), ]

fasta_entries <- character(nrow(halamphora_mol))

# Loop through each row of the dataframe
for (i in 1:nrow(halamphora_mol)) {
  # Extract ASVid and Domain:Species information
  asv_id <- halamphora_mol$ASVid[i]
  domain_species <- paste(halamphora_mol[i, "Domain"], halamphora_mol[i, "Phylum"], halamphora_mol[i, "Class"], 
                          halamphora_mol[i, "Order"], halamphora_mol[i, "Genus"], halamphora_mol[i, "Species"], sep = ":")
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

# Write the FASTA entries to a file
writeLines(fasta_entries, "primer_mismatches/halamphora_RBCL.fasta")

# get just cocconeis diat barcode
diat_barcode <- read.csv("imput_data/v10_correspondence.csv")

only_hal <- diat_barcode |>
  filter(genus == "Halamphora")

# write a fasta file
fasta_entries <- character(nrow(only_hal))

# Loop through each row of the dataframe
for (i in 1:nrow(only_hal)) {
  # Extract ASVid and Domain:Species information
  asv_id <- only_hal$sequence[i]
  domain_species <- paste(only_hal[i, "domain"], only_hal[i, "phylum"], only_hal[i, "class"], 
                          only_hal[i, "order"], only_hal[i, "genus"], only_hal[i, "species_compromise"], sep = ":")
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

# Write the FASTA entries to a file
writeLines(fasta_entries, "primer_mismatches/halamphora_diat.fasta")

library(ggpubr)

# read in the primer mismatch csv
primer_mismatch <- read.csv("sequence_mismatches/sequence_mismatches.csv")

# do a box plot
seq_mismatch <- ggplot(primer_mismatch, aes(x = category, y = avg_mismatches)) +
  geom_boxplot(alpha = 0.2, color = "gray", outlier.alpha = 0) +  
  geom_jitter(width = 0.1, height = 0, alpha = 0.8, color = "black") +
  xlab("Category (high in the Illumina data vs under-represented)") +
  ylab("Average number of sequence mismatches per genus") +
  theme_classic() +
  annotate("text", x = c(1, 2), y = c(10.75, 14), label = c("a", "b"), size = 4, color = 'darkred') +
  theme(axis.text.x = element_text(size = 12),   # Adjust x-axis label size
        axis.text.y = element_text(size = 12))
 
seq_mismatch 

# t-test
high_illumina_data <- primer_mismatch |>
  filter(category == "high_illumina")

under_represented_data <- primer_mismatch |>
  filter(category == "under_represented")

# Perform t-test
t_test_result <- t.test(high_illumina_data$avg_mismatches, under_represented_data$avg_mismatches)
t_test_result # 0.01
