###### SET UP #######
## load packages
library(tidyverse)
library(dplyr)
library(qualpalr)
library(ggpubr)
library(ggh4x)
library(indicspecies)
library(data.table)
library(reshape2)
library(ggplot2); theme_set(theme_bw()+
                              theme(panel.grid = element_blank(),
                                    strip.background = element_rect(fill="white"),
                                    axis.text.y = element_text(size = 10, colour = "black"),
                                    axis.title = element_text(size=10, face="bold"),
                                    strip.text = element_text(color="black", size=10),
                                    legend.text=element_text(size=10),
                                    axis.line = element_line(colour = "black"),
                                    axis.text.x = element_blank()))

## tell R whereto get data
setwd("/Users/andreajackman/R_Stuff/Parfrey/diatoms_new/counting_data")

## read in data
df <- read.csv("diatom_sum_per_box.csv")

# Define columns to analyze
columns_to_analyze <- c("Navicula", "Tabularia", "Cylindrotheca", "Pseudonitzschia", "Gomphonemopsis", "Hyalodiscus", "Coscinodiscus",
                        "Plagiotropis", "Nitzschia", "Rhoicosphenia", "Pseudogomphonema", "Thalassiosira", "Diploneis", "Pariraphis",
                        "Planothidium", "Halamphora", "Achnanthes", "Amphora", "Bacillaria", "Ellerbeckia", "Minidiscus", "Skeletonema",
                        "Gomphoseptatum", "Cocconeis")

#columns_to_analyze <- c("Navicula", "Tabularia", "Cylindrotheca", "Gomphonemopsis", "Nitzschia", "Rhoicosphenia", "Pseudogomphonema", "Planothidium", "Cocconeis")

calculate_proportion <- function(df, columns) {
  proportions <- sapply(columns, function(col) mean(df[[col]] > 0))
  return(proportions)
}

# Calculate proportions for each blade position
blade_positions <- unique(df$blade_position)
proportions_list <- lapply(blade_positions, function(bp) {
  df_bp <- subset(df, blade_position == bp)
  proportions <- calculate_proportion(df_bp, columns_to_analyze)
  proportions_df <- data.frame(blade_position = bp, species = names(proportions), proportion = proportions)
  return(proportions_df)
})

# Combine results
result <- do.call(rbind, proportions_list)

# Determine the blade section with the highest proportion for each species
max_proportion <- result %>%
  group_by(species) %>%
  summarize(max_proportion = max(proportion))

max_proportion_bp <- result %>%
  group_by(species) %>%
  filter(proportion == max(proportion)) %>%
  select(species, blade_position) %>%
  distinct()

# Add the max blade position to the result data
result <- merge(result, max_proportion_bp, by = "species", suffixes = c("", "_max"))

species_to_change <- c("Cocconeis", "Navicula", "Tabularia", "Cylindrotheca")
result$blade_position_max[result$species %in% species_to_change] <- 'T'

# Order species within each facet by prevalence
result <- result %>%
  #group_by(blade_position_max) %>%
  arrange(proportion) %>%
  mutate(species = factor(species, levels = unique(species)))

# Reshape data for heatmap
result_melt <- melt(result[, c("blade_position", "species", "proportion", "blade_position_max")], id.vars = c("blade_position", "species", "blade_position_max"), variable.name = "variable", value.name = "value")

# Extract the values at blade_position = "T"
value_at_T <- result_melt %>%
  filter(blade_position == "T") %>%
  select(species, value) %>%
  arrange(value)

# Order species by the value at blade_position = "T"
unique_species_order <- unique(value_at_T$species)
result_melt$species <- factor(result_melt$species, levels = unique_species_order)

# Create heatmap sorted by value at blade_position = "T"
ggplot(result_melt, aes(x = blade_position, y = species, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "mediumorchid4") +
  labs(x = "Blade Position", y = "Genus", fill = "Proportion of \nSamples Present") +
  theme_classic()

# Calculate mean prevalence for ordering facets
mean_prevalence <- result %>%
  group_by(blade_position_max) %>%
  summarize(mean_prevalence = mean(proportion))

# Order blade_position_max by mean_prevalence
#result_melt$blade_position_max <- factor(result_melt$blade_position_max, levels = mean_prevalence$blade_position_max[order(mean_prevalence$mean_prevalence, decreasing = FALSE)])

custom_labels <- c(B = "Base", M = "Middle", T = "Tip")

# Create heatmap with faceting by blade section of highest proportion
ggplot(result_melt, aes(x = blade_position, y = species, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "mediumorchid4") +
  labs(x = "Blade Position", y = "Genus", fill = "Proportion of \nSamples Present") +
  theme_classic()
  #facet_wrap(~ blade_position_max, scales = "free_y", labeller = labeller(blade_position_max = custom_labels))


# Calculate the count of each genus found per blade section
genus_individuals_per_section <- df %>%
  group_by(blade_position) %>%
  dplyr::summarize(across(all_of(columns_to_analyze), sum)) %>%
  pivot_longer(cols = -blade_position, names_to = "species", values_to = "individuals")

result <- merge(result, genus_individuals_per_section, by = c("blade_position", "species"))

# Reshape data for heatmap
result_melt <- melt(result[, c("blade_position", "species", "proportion", "blade_position_max", "individuals")], id.vars = c("blade_position", "species", "blade_position_max", "individuals"), variable.name = "variable", value.name = "value")

# Calculate mean prevalence for ordering facets
mean_prevalence <- result %>%
  group_by(blade_position_max) %>%
  summarize(mean_prevalence = mean(proportion))

# Order blade_position_max by mean_prevalence
result_melt$blade_position_max <- factor(result_melt$blade_position_max, levels = mean_prevalence$blade_position_max[order(mean_prevalence$mean_prevalence, decreasing = FALSE)])

custom_labels <- c(B = "Base", M = "Middle", T = "Tip")

# Create heatmap with faceting by blade section of highest proportion
ggplot(result_melt, aes(x = blade_position, y = species, fill = value)) +
  geom_tile()+
  scale_fill_gradient(low = "white", high = "mediumorchid4") +
  labs(x = "Blade Position", y = "Genus", fill = "Proportion of \nSamples Present") +
  theme_classic() +
  geom_text(aes(label = individuals), color = "black", size = 3) +
  facet_wrap(~ blade_position_max, scales = "free_y")

# Calculate mean proportion at the Tip (T) position
mean_proportion_T <- result %>%
  filter(blade_position_max == "T") %>%
  group_by(species) %>%
  summarize(mean_proportion_T = mean(proportion)) %>%
  arrange(desc(mean_proportion_T))

# Order species by mean proportion at the Tip
result_melt$species <- factor(result_melt$species, levels = mean_proportion_T$species)

# Create heatmap with faceting by blade section of highest proportion
ggplot(result_melt, aes(x = blade_position, y = species, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "mediumorchid4") +
  labs(x = "Blade Position", y = "Genus", fill = "Proportion of \nSamples Present") +
  theme_classic() +
  geom_text(aes(label = individuals), color = "black", size = 3) +
  facet_wrap(~ blade_position_max, scales = "free_y")



