# Difference in percent dissimilarity between well and under represented taxa
# Andrea Jackman


library(ggpubr)

# read in the seq mismatch csv
primer_mismatch <- read.csv("molecular_data/primer_mismatches/sequence_mismatches.csv")

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
