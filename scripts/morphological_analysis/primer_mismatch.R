# read in libraries
library(ggpubr)
library(dplyr)
library(ggplot2)

# read in the primer mismatch csv
f_primer_mismatch <- read.csv("counting_data/Primer_mismatches_forward_avg.csv")
r_primer_mismatch <- read.csv("counting_data/Primer_mismatches_reverse_avg.csv")

# do a box plot
plot <- ggplot(f_primer_mismatch, aes(x = well_under_represented, y = avg_r_mismatches, label = Genus)) +
  geom_boxplot(alpha = 0.2, color = "gray") +  
  geom_point(size = 1) +
  geom_text(size = 3, position = position_dodge(width = 0.75), vjust = -0.75) +
  xlab("Category (high in the Illumina data vs under-represented)") +
  ylab("Average number of mismatches per Genus") 
  # annotate("text", x = c(1, 2), y = c(10.75, 14), label = c("a", "b"), size = 4, color = 'darkred')

plot 


# t-test
high_illumina_data <- f_primer_mismatch |>
  filter(well_under_represented == 'high_illumina')

under_represented_data <- f_primer_mismatch |>
  filter(well_under_represented == "under_represented")

# Perform t-test
t_test_result <- t.test(high_illumina_data$avg_r_mismatches, under_represented_data$avg_r_mismatches)
t_test_result # 0.794

# do a box plot
plot <- ggplot(r_primer_mismatch, aes(x = well_under_represented, y = avg_r_mismatches, label = Genus)) +
  geom_boxplot(alpha = 0.2, color = "gray") +  
  geom_point(size = 1) +
  geom_text(size = 3, position = position_dodge(width = 0.75), vjust = -0.75) +
  xlab("Category (high in the Illumina data vs under-represented)") +
  ylab("Average number of mismatches per Genus") 
# annotate("text", x = c(1, 2), y = c(10.75, 14), label = c("a", "b"), size = 4, color = 'darkred')

plot 


# t-test
high_illumina_data <- r_primer_mismatch |>
  filter(well_under_represented == 'high_illumina')

under_represented_data <- r_primer_mismatch |>
  filter(well_under_represented == "under_represented")

# Perform t-test
t_test_result <- t.test(high_illumina_data$avg_r_mismatches, under_represented_data$avg_r_mismatches)
t_test_result # 0.6298


r_mismatches <- ggplot(r_primer_mismatch, aes(x = well_under_represented, y = avg_r_mismatches)) +
  geom_boxplot(alpha = 0.2, color = "gray", outlier.alpha = 0) +  # Transparent boxplot
  geom_jitter(width = 0.1, height = 0, alpha = 0.8, color = "black") +
  xlab("Category (high in the Illumina data vs under-represented)") +
  ylab("Average number of reverse primer mismatches per genus") +
  annotate("text", x = c(1, 2), y = c(1.5, 1.1), label = c("a", "a"), size = 5, color = 'darkred') +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14),   
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17))

r_mismatches

f_mismatches <- ggplot(f_primer_mismatch, aes(x = well_under_represented, y = avg_r_mismatches)) +
  geom_boxplot(alpha = 0.2, color = "gray", outlier.alpha = 0) +  # Transparent boxplot
  geom_jitter(width = 0.1, height = 0, alpha = 0.8, color = "black") +
  xlab("Category (high in the Illumina data vs under-represented)") +
  ylab("Average number of forward primer mismatches per genus") +
  theme_classic() +
  annotate("text", x = c(1, 2), y = c(2.1, 3.1), label = c("a", "a"), size = 5, color = 'darkred') +
  theme(axis.text.x = element_text(size = 14),   
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17))

f_mismatches

primer_mismatch <- read.csv("sequence_mismatches/primer_mismatches.csv")

# do a box plot
seq_mismatch <- ggplot(primer_mismatch, aes(x = category, y = avg_mismatches)) +
  geom_boxplot(alpha = 0.2, color = "gray", outlier.alpha = 0) +  
  geom_jitter(width = 0.1, height = 0, alpha = 0.8, color = "black") +
  xlab("Category (high in the Illumina data vs under-represented)") +
  ylab("Average number of sequence mismatches per genus") +
  theme_classic() +
  annotate("text", x = c(1, 2), y = c(10.75, 14), label = c("a", "b"), size = 5, color = 'darkred') +
  theme(axis.text.x = element_text(size = 14),   
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17))

seq_mismatch 

ggarrange(f_mismatches,r_mismatches, seq_mismatch, ncol=3, nrow=1, labels = c("A", "B", "C"))
