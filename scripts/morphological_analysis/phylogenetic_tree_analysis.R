library(Biostrings)

# read in the data
data <- read.csv("tree_stuff/2024-06-17_example_tree.csv", stringsAsFactors = FALSE)

# make sure columns are in the right format
str(data)

# combine the two columns into one
data$Species <- gsub(" ", "_", data$Species)
data$Identifier <- paste(data$Sequence.ID, data$Order, data$Species, sep = "_")

clean_sequences <- function(seqs) {
  cleaned_seqs <- gsub("\n", "", seqs)  # Remove newline characters
  return(cleaned_seqs)
}

# Create DNAStringSet object from the sequences
data$sequence <- clean_sequences(data$Sequence)

# Create DNAStringSet object
sequences <- DNAStringSet(data$sequence)

# Assign names to the sequences
names(sequences) <- data$Identifier

# Name output file
output_file <- "tree_stuff/2024-06-17_example_tree.fasta"  

# Make FASTA file
writeXStringSet(sequences, filepath = output_file)

