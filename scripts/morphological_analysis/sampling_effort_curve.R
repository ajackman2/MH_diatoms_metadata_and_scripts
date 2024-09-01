# import libraries
library(ggplot2)
library(dplyr)  
library(vegan)

# read in data
diatom_data <- read.csv("counting_data/master_diatom_data.csv")

# filter for box 9 nad 0_19_7
diatom_filtered <- diatom_data |>
  filter(box %in% c('9', '0_19_7')) |>
  select(-grid_position, -under_review, -blade_position, -folder, - suggested_id, -file_name)

# make image_num numeric
diatom_filtered$image_num <- as.numeric(diatom_filtered$image_num)
diatom_filtered <- mutate_at(diatom_filtered, vars(Navicula:Undatella), as.numeric)

# make a species matrix
species_matrix <- as.matrix(select(diatom_filtered, Navicula:Undatella) > 0)

# create species accumulation curve
sac <- specaccum(species_matrix)

# plot
plot(sac, xlab = "Number of Samples", ylab = "Number of Species",
     main = "Species Accumulation Curve")

species_matrix <- as.matrix(select(diatom_filtered, Navicula:Undatella) > 0)

# Calculate species accumulation curves for each box category
sac_list <- diatom_filtered %>%
  group_by(box) %>%
  summarize(species_accumulation_curve = list(specaccum(species_matrix)))

# Combine the species accumulation curves into a data frame
sac_df <- data.frame(do.call(rbind, sac_list$species_accumulation_curve))
sac_df$Box <- rep(sac_list$box, each = nrow(sac_df))

# Plot species accumulation curves with confidence intervals
ggplot(sac_df, aes(x = sites, y = specimens, color = box, group = box)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Box), alpha = 0.3) +
  labs(x = "Number of Samples", y = "Number of Species",
       title = "Species Accumulation Curve with 95% Confidence Intervals") +
  theme_minimal()
