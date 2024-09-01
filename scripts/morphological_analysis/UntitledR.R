# get just cocconeis diat barcode
diat_barcode <- read.csv("imput_data/v10_correspondence.csv")

only_hal <- diat_barcode |>
  filter(genus == "Entomoneis")

only_hal <- only_hal |>
  unique(only_hal$ID_seq_original)

# write a fasta file
fasta_entries <- character(nrow(only_hal))

# Loop through each row of the dataframe
for (i in 1:nrow(only_hal)) {
  # Extract ASVid and Domain:Species information
  asv_id <- only_hal$sequence[i]
  domain_species <- only_hal$ID_seq_original[i]
  
  # Create the FASTA entry
  fasta_entry <- paste0(">", domain_species, "\n", asv_id, "\n")
  
  # Add the FASTA entry to the vector
  fasta_entries[i] <- fasta_entry
}

# Write the FASTA entries to a file
writeLines(fasta_entries, "sequence_mismatches/ento_diat2.fasta")
