# plots for the overlap between molecular and morphological data

# load libraries
library(dplyr)
library(plyr)
library(tidyverse)
library(tidyr)
library(stringr)
library(ggplot2); theme_set(theme_bw()+
                              theme(panel.grid = element_blank(),
                                    strip.background = element_rect(fill="white"),
                                    axis.text.y = element_text(size = 12, colour = "black"),
                                    axis.title = element_text(size=15, face="bold"),
                                    strip.text = element_text(color="black", size=10),
                                    legend.text=element_text(size=10),
                                    axis.line = element_line(colour = "black")))


# read the data in
overlap_data <- read_csv("counting_data/combined_sem_illumina_counts2.csv")

# make a new dataframe and remove higher level classifications
plot_data <- overlap_data %>%
  filter(!(Genus %in% c('Achnanthales', 'Bacillariales', 'Bacillariophyceae', 'Bacillariophytina', 'Bacillariophyta', 'Coscinodiscophytina', 
                        'Coscinodiscophyceae', 'Cymbellales', 'Diatomea', 'Fragilariales', 'Fragilariophyceae', 'Mediophyceae',
                        'Melosirids', 'Naviculales', 'Thalassiophysales', 'Thalassiosirales', 'Ochrophyta')))

# remove not found in march from morphological data
plot_data <- plot_data %>%
  filter(!(Genus %in% c('Undatella', 'Actinoptychus', 'Trachyneis', 'Rhopalodia',
                        'Hobaniella', 'Hanzschia',
                        'Gyrosigma', 'Cyclotella')))

# collapse this
plot_data <- plot_data %>%
  select(-...1) %>%
  mutate(Genus = ifelse(Genus == "Paribellus", "Parlibellus", Genus))

new_values <- tibble(
  Genus = "Parlibellus",
  sum.mol.ra = 9.187296e-02,
  mean.mol.ra = 2.262881e-04,
  n.mol.samples = 29,
  in18S = 'no',
  inRBCL = 'yes',
  sum.sem.counts = 0,
  n.sem.samples = 17,
  mean.sem.count = 0,
  in_counts = 'excluded'
  
)

# Filter out the original Parlibellus rows
filtered_data <- plot_data %>%
  filter(Genus != "Parlibellus")

# Combine the filtered data with the new values for Parlibellus
plot_data <- bind_rows(filtered_data, new_values)
  

# fix Asteromphalus
plot_data <- plot_data |>
  mutate(in_counts = ifelse(Genus == 'Asteromphalus', 'excluded', in_counts))

# replace all of the zeros
plot_data <- plot_data %>%
  mutate(across(everything(), ~replace(., is.na(.), 0)))

# make new cols for color classification
plot_data$in18S <- ifelse(is.na(plot_data$in18S), "no", plot_data$in18S)
plot_data$inRBCL <- ifelse(is.na(plot_data$inRBCL), "no", plot_data$inRBCL)
plot_data$in_molecular <- ifelse(plot_data$in18S == "yes" | plot_data$inRBCL == "yes", "yes", "no")
plot_data$in_counts <- ifelse(is.na(plot_data$in_counts), "no", plot_data$in_counts)
plot_data$in_counts <- ifelse(plot_data$in_counts == "excluded" | plot_data$in_counts == "yes", "yes", "no")

plot_data$where_found <- NA
plot_data$where_found <- ifelse(plot_data$in_counts == "yes" & plot_data$in_molecular == "yes", "in_both",
                                ifelse(plot_data$in_counts == "yes" & plot_data$in_molecular == "no", "in_morpho",
                                       ifelse(plot_data$in_counts == "no" & plot_data$in_molecular == "yes", "in_molecular", NA)))

# correct the abundances
total_morphological_sum <- sum(plot_data$mean.sem.count, na.rm = TRUE)


# calculate per-sample relative abundance
plot_data$sem.ra <- plot_data$mean.sem.count / total_morphological_sum

# replace NA with zero
plot_data$sem.ra[is.na(plot_data$sem.ra)] <- 0
plot_data$sum.mol.ra[is.na(plot_data$sum.mol.ra)] <- 0

filtered_plot_data <- plot_data[plot_data$sum.mol.ra != 0 & plot_data$sem.ra != 0,]

# Create the scatter plot with the filtered data and add the trendline
scatter_plot3 <- ggplot(plot_data, aes(x = sem.ra, y = sum.mol.ra, color = where_found)) +
  geom_point(data = filtered_plot_data) +
  geom_smooth(data = filtered_plot_data, method = 'lm') +
  geom_point(data = plot_data[plot_data$sum.mol.ra == 0 | plot_data$sem.ra == 0,]) +
  labs(title = "Comparison of Molecular and Morphological Abundance",
       x = "Morphological Relative Abundance", y = "Molecular Relative Abundance")

scatter_plot3
# plot the data
scatter_plot <- ggplot(plot_data, aes(x = sem.ra, y = sum.mol.ra, label = Genus, color = where_found)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(title = "Comparison of Molecular and Morphological Abundance",
       x = "Morphological Relative Abundance", y = "Molecular Relative Abundance") 
  #coord_cartesian(xlim = c(0, 0.03), ylim = c(0, 2))

scatter_plot

lm <- lm(sum.mol.ra~sem.ra, data = filtered_plot_data)
summary(lm) # not significant

# do an ANOVA between the groups
anova_result <- aov(sum.mol.ra ~ where_found, data = plot_data)

# summary of ANOVA
summary(anova_result) # 0.111
TukeyHSD(anova_result)

# do an ANOVA between the groups
anova_result <- aov(sem.ra ~ where_found, data = plot_data)

# summary of ANOVA
summary(anova_result) # 0.13

# make a scatterplot just for in_botj
plot_data2 <- plot_data |>
  filter(where_found == 'in_both')

scatter_plot2 <- ggplot(plot_data2, aes(x = sem.ra, y = sum.mol.ra)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(title = "Comparison of Molecular and Morphological Abundance",
       x = "Morphological Relative Abundance", y = "Molecular Relative Abundance")

scatter_plot2

lm <- lm(sum.mol.ra ~ sem.ra, data = plot_data2)
summary(lm) # totally not significant 0.441

# try with relative rank instead

plot_data2$mol_ra_rank <- rank(-plot_data2$sum.mol.ra, ties.method = "random")
plot_data2$sem_ra_rank <- rank(-plot_data2$sem.ra, ties.method = "random")

scatter_plot3 <- ggplot(plot_data2, aes(x = sem_ra_rank, y = mol_ra_rank)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(title = "Comparison of Molecular and Morphological Abundance",
       x = "Morphological Rank", y = "Molecular Rank") +
  annotate("text", x = min(plot_data2$sem_ra_rank) +7, y = max(plot_data2$mol_ra_rank)-1.4,
           label = paste("p-value: 0.006**", "\ny = 0.40x + 9.0296"), 
           hjust = 1, vjust = 1, size = 4)

scatter_plot3

lm <- lm(mol_ra_rank ~ sem_ra_rank, data = plot_data2)
summary(lm) # 0.008356

# try with the whole data

plot_data$mol_ra_rank <- rank(plot_data$sum.mol.ra, ties.method = "random")
plot_data$sem_ra_rank <- rank(plot_data$sem.ra, ties.method = "random")

scatter_plot3 <- ggplot(plot_data, aes(x = sem_ra_rank, y = mol_ra_rank)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(title = "Comparison of Molecular and Morphological Abundance",
       x = "Morphological Rank", y = "Molecular Rank") +
  annotate("text", x = min(plot_data$sem_ra_rank) +15, y = max(plot_data$mol_ra_rank)-1,
           label = paste("p-value: 0.001**", "\ny = 0.32x + 32.92"), 
           hjust = 1, vjust = 1, size = 4)

scatter_plot3

lm <- lm(mol_ra_rank ~ sem_ra_rank, data = plot_data)
summary(lm) #0.007618

# try with the ranks and diatoms not found as zero
plot_data <- plot_data %>%
  dplyr::mutate(across(everything(), ~replace(., is.na(.), 0)))

plot_data$mol_ra_rank <- ifelse(plot_data$n.mol.samples == 0, 0, rank(-plot_data$sum.mol.ra, ties.method = "random"))
plot_data$sem_ra_rank <- ifelse(plot_data$sum.sem.counts == 0 & plot_data$where_found != "in_both", 0, rank(-plot_data$sem.ra, ties.method = "random"))

filtered_plot_data <- plot_data |>
  filter(where_found == "in_both")

# Create the scatter plot with the filtered data and add the trendline
scatter_plot3 <- ggplot(plot_data, aes(x = sem_ra_rank, y = mol_ra_rank, color = where_found)) +
  geom_point(data = filtered_plot_data) +
  geom_smooth(data = filtered_plot_data, method = 'lm') +
  geom_point(data = plot_data[plot_data$sum.mol.ra == 0 | plot_data$sem.ra == 0,]) +
  labs(x = "Morphological Rank", y = "Molecular Rank") +
  annotate("text", x = min(plot_data$sem_ra_rank) +20, y = max(plot_data$mol_ra_rank)-1,
           label = paste("p-value: 0.07", "\ny = 0.24x + 15.44"), 
           hjust = 1, vjust = 1, size = 4) +
  scale_color_manual(labels = c("Both", "Molecular only", "Morphological only"),
                     values = c('purple3', 'red3', 'blue3')) +
  guides(color=guide_legend(title="Found in Morphological or Molecular Data?"))

scatter_plot3

lm <- lm(mol_ra_rank ~ sem_ra_rank, data = filtered_plot_data)
summary(lm) # 0.07

#### make a venn diagram
# read the data in
overlap_data <- read_csv("counting_data/combined_sem_illumina_counts2.csv")

# make a new dataframe and remove higher level classifications
# plot_data <- overlap_data %>%
#   filter(!(Genus %in% c('Achnanthales', 'Bacillariales', 'Bacillariophyceae', 'Bacillariophytina', 'Bacillariophyta', 'Coscinodiscophytina', 
#                         'Coscinodiscophyceae', 'Cymbellales', 'Diatomea', 'Fragilariales', 'Fragilariophyceae', 'Mediophyceae',
#                         'Melosirids', 'Naviculales', 'Thalassiophysales', 'Thalassiosirales', 'Ochrophyta')))

# remove not found in march from morphological data
plot_data <- overlap_data %>%
  filter(!(Genus %in% c('Undatella', 'Actinoptychus', 'Trachyneis', 'Rhopalodia',
                        'Hobaniella', 'Hanzschia',
                        'Gyrosigma', 'Cyclotella')))

# collapse this
plot_data <- plot_data %>%
  select(-...1) %>%
  mutate(Genus = ifelse(Genus == "Paribellus", "Parlibellus", Genus))

new_values <- tibble(
  Genus = "Parlibellus",
  sum.mol.ra = 9.187296e-02,
  mean.mol.ra = 2.262881e-04,
  n.mol.samples = 29,
  in18S = 'no',
  inRBCL = 'yes',
  sum.sem.counts = 0,
  n.sem.samples = 17,
  mean.sem.count = 0,
  in_counts = 'excluded'
  
)

# Filter out the original Parlibellus rows
filtered_data <- plot_data %>%
  filter(Genus != "Parlibellus")

# Combine the filtered data with the new values for Parlibellus
plot_data <- bind_rows(filtered_data, new_values)


# fix Asteromphalus
plot_data <- plot_data |>
  mutate(in_counts = ifelse(Genus == 'Asteromphalus', 'excluded', in_counts))

# include Mark's data
library(VennDiagram)

ven_data2 <- plot_data
ven_data2$in18S <- ifelse(is.na(ven_data2$in18S), "no", ven_data2$in18S)
ven_data2$inRBCL <- ifelse(is.na(ven_data2$inRBCL), "no", ven_data2$inRBCL)
ven_data2$in_counts <- ifelse(is.na(ven_data2$in_counts), "no", ven_data2$in_counts)
ven_data2$in_marks <- ifelse(ven_data2$in_counts %in% c("excluded", "yes"), "yes", "no")
ven_data2$in_counts <- gsub("excluded", "no", ven_data2$in_counts)

# make a new col called in_molecular
ven_data2$in_molecular <- ifelse(ven_data2$in18S == "yes" | ven_data2$inRBCL == "yes", "yes", "no")

# make a new vec with the in_molecular and in_counts
molecular_genera_vec2 <- c()
morpho_genera_vec2 <- c()
marks_genera_vec <- c()

for (i in 1:nrow(ven_data2)) {
  if (ven_data2$in_molecular[i] == "yes") {
    molecular_genera_vec2 <- c(molecular_genera_vec2, ven_data2$Genus[i])
  }
}

for (i in 1:nrow(ven_data2)) {
  if (ven_data2$in_counts[i] == "yes") {
    morpho_genera_vec2 <- c( morpho_genera_vec2, ven_data2$Genus[i])
  }
}

for (i in 1:nrow(ven_data2)) {
  if (ven_data2$in_marks[i] == "yes") {
    marks_genera_vec <- c(marks_genera_vec, ven_data2$Genus[i])
  }
}

venn.plot2 <- venn.diagram(
  x = list(Molecular = molecular_genera_vec2, Morphological = morpho_genera_vec2, Mark = marks_genera_vec),
  category.names = c("Molecular", "Morphological Count", "Mark's Count"),
  filename = NULL
)


grid.draw(venn.plot2)

