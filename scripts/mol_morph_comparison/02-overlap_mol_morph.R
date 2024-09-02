# look at the blast overlap

# libraries
library(dplyr)
library(plyr)
library(tidyverse)
library(tidyr)
library(stringr)
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

# set the working directory
setwd("/Users/andreajackman/R_Stuff/Parfrey/MH_diatoms_metadata_and_scripts")

# read the data in
mol_18S <- read.csv("molecular_data/tree_corrected_mol_18S")
mol_RBCL <- read.csv("molecular_data/tree_corrected_mol_rbcl")
diatom_data <- read.csv("morphological_data/master_diatom_data.csv")

# fix the diatom data
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

# add a row with synedra
new_row <- data.frame(
  Genus = "Synedra",
  sum.sem.counts = 0,
  n.sem.samples = 17,
  mean.sem.count = 0,
  in_counts = "excluded"
)

# Adding the new row
diatom_data <- rbind(diatom_data, new_row)

diatom_data$Genus[diatom_data$Genus == "Pseudonitzschia"] <- "Pseudo-nitzschia"

## get per sample relative abundance
mol_RBCL <- mol_RBCL |>
  filter(!(Genus %in% c('Aulacodiscus', 'Bolidomonas', 'Biddulphiopsis', 'Caloneis', 'Cerataulina', 'Cymatosira',
                        'Cyclotella', 'Detonula', 'Diatoma', 'Discostella', 'Gomphonema', 'Guinardia', 'Lauderia',
                        'Pauliella', 'Placoneis', 'Porosira', 'Psammogramma', 'Rhoikoneis', 'Shionodiscus',
                        'Staurosira', 'Stephanodiscus', 'Hyalosira', 'Nanofrustulum', 'Pteridomonas')))

mol_18S <- mol_18S |>
  filter(!(genus %in% c('Aulacodiscus', 'Bolidomonas', 'Biddulphiopsis', 'Caloneis', 'Cerataulina', 'Cymatosira',
                        'Cyclotella', 'Detonula', 'Diatoma', 'Discostella', 'Gomphonema', 'Guinardia', 'Lauderia',
                        'Pauliella', 'Placoneis', 'Porosira', 'Psammogramma', 'Rhoikoneis', 'Shionodiscus',
                        'Staurosira', 'Stephanodiscus', 'Hyalosira', 'Nanofrustulum', 'Pteridomonas')))

mol_RBCL$ra = as.numeric(mol_RBCL$asv_abundance)/as.numeric(mol_RBCL$rd_filt)
mol_18S$ra = as.numeric(mol_18S$asv_abundance)/as.numeric(mol_18S$rd_filt)

# change plant_id to character
mol_18S$plant_id = as.character(mol_18S$leaf_number)
mol_RBCL$plant_id = as.character(mol_RBCL$leafnumber)

## rename
names(mol_18S)[names(mol_18S)=="genus"]<-"Genus"

inRBCL = c(unique(mol_RBCL$Genus))
in18S = c(unique(mol_18S$Genus))

mol = full_join(mol_RBCL, mol_18S)

mol$ra[is.na(mol$ra)] <- 0

mol = ddply(mol, c("Genus"),
            summarise,
            sum.mol.ra = sum(ra),
            mean.mol.ra = mean(ra),
            n.mol.samples =length(unique(Row.names)))

mol$in18S = ifelse(mol$Genus %in% c(in18S), "yes", "no")
mol$inRBCL = ifelse(mol$Genus %in% c(inRBCL), "yes", "no")

# remove the genera that we are not confident are there

in_one <- c('ASteromphalus', 'Coscinodiscus', 'Dimeregramma', 'Ditylum', 'Ellerbeckia', 'Entomoneis', 'Fragilariopsis', 'Fragilaria', 
            'Gomphonemopsis', 'Grammatophora', 'Leptocylindrus', 'Lyrella', 'Paralia', 'Plagiogramma',
            'Psammodictyon', 'Pseudo-nitzschia', 'Rhizosolenia', 'Trachyneis', 'Opephora', 'Fogedia')

mol$one_sample <- NA
mol <- mol %>%
  mutate(one_sample = ifelse(Genus %in% in_one, "yes", "no"))

###### combine data #####
all = full_join(mol, diatom_data)


#write.csv(all, "counting_data/combined_sem_illumina_counts_blast.csv")

#### make the venn diagram ####
# remove not found in march from morphological data
plot_data <- all %>%
  filter(!(Genus %in% c('Undatella', 'Actinoptychus', 'Rhopalodia',
                        'Hobaniella','Cyclotella', 'Donkinia')))

# add proskinia
new_values <- c(
  sum.sem.counts = 0,
  n.sem.samples = 17,
  mean.sem.count = 0,
  in_counts = "excluded"
)

# Find the row where Genus is 'Proschkinia'
row_index <- which(plot_data$Genus == "Proschkinia")

# Modify the values in that row
plot_data[row_index, c('sum.sem.counts', 'n.sem.samples', 'mean.sem.count', 'in_counts')] <- new_values

# collapse this
plot_data <- plot_data %>%
  mutate(Genus = ifelse(Genus == "Paribellus", "Parlibellus", Genus))

new_values <- tibble(
  Genus = "Parlibellus",
  sum.mol.ra = 9.181825e-02,
  mean.mol.ra =1.334568e-04,
  n.mol.samples = 43,
  in18S = 'no',
  inRBCL = 'yes',
  one_sample = 'no',
  sum.sem.counts = 0,
  n.sem.samples = 17,
  mean.sem.count = 0,
  in_counts = 'excluded'
  
)

# Filter out the original Parlibellus rows
filtered_data <- plot_data %>%
  filter(Genus != "Parlibellus")

# Combine the filtered data with the new values for Parlibellus
filtered_data$mean.sem.count = as.numeric(filtered_data$mean.sem.count)
filtered_data$sum.sem.counts = as.numeric(filtered_data$sum.sem.counts)
filtered_data$n.sem.samples = as.numeric(filtered_data$n.sem.samples)
plot_data <- bind_rows(filtered_data, new_values)


# fix Asteromphalus
plot_data <- plot_data |>
  mutate(in_counts = ifelse(Genus == 'Asteromphalus', 'excluded', in_counts))

plot_data <- plot_data |>
  filter(Genus != 'Paribellus')

# add new genera
new_values <- tibble(
  Genus = "Asteromphalus",
  sum.mol.ra = 4.822754e-03,
  mean.mol.ra =1.121571e-04,
  n.mol.samples = 43,
  in18S = 'no',
  inRBCL = 'yes',
  one_sample = 'yes',
  sum.sem.counts = 0,
  n.sem.samples = 17,
  mean.sem.count = 0,
  in_counts = 'excluded'
  
)

plot_data <- bind_rows(plot_data, new_values)

new_values <- tibble(
  Genus = "Eunotia",
  sum.mol.ra = 0,
  mean.mol.ra =0,
  n.mol.samples = 0,
  in18S = 'no',
  inRBCL = 'no',
  one_sample = 'no',
  sum.sem.counts = 0,
  n.sem.samples = 17,
  mean.sem.count = 0,
  in_counts = 'excluded'
  
)

plot_data <- bind_rows(plot_data, new_values)

plot_data <- plot_data %>%
  filter(Genus != "Fragilaria")

new_values <- tibble(
  Genus = "Fragilaria",
  sum.mol.ra = 0,
  mean.mol.ra =0,
  n.mol.samples = 0,
  in18S = 'yes',
  inRBCL = 'no',
  one_sample = 'no',
  sum.sem.counts = 0,
  n.sem.samples = 17,
  mean.sem.count = 0,
  in_counts = 'excluded'
  
)

plot_data <- bind_rows(plot_data, new_values)


#write.csv(plot_data, "counting_data/mol_morph_overlap_blast_updated.csv")

# include Mark's data
library(VennDiagram)

ven_data2 <- plot_data
ven_data2$in18S <- ifelse(is.na(ven_data2$in18S), "no", ven_data2$in18S)
ven_data2$inRBCL <- ifelse(is.na(ven_data2$inRBCL), "no", ven_data2$inRBCL)
ven_data2$in_counts <- ifelse(is.na(ven_data2$in_counts), "no", ven_data2$in_counts)
ven_data2$in_marks <- ifelse(ven_data2$in_counts %in% c("excluded", "yes"), "yes", "no")
ven_data2$in_counts <- gsub("excluded", "no", ven_data2$in_counts)

# make a new col called in_molecular
ven_data2$in_18S <- ifelse(ven_data2$in18S == "yes" , "yes", "no")
ven_data2$in_rbcl <- ifelse(ven_data2$inRBCL == "yes" , "yes", "no")

# make a new vec with the in_molecular and in_counts
molecular_genera_vec2 <- c()
molecular_genera_vec3 <- c()
morpho_genera_vec2 <- c()
marks_genera_vec <- c()

for (i in 1:nrow(ven_data2)) {
  if (ven_data2$in_18S[i] == "yes") {
    molecular_genera_vec2 <- c(molecular_genera_vec2, ven_data2$Genus[i])
  }
}

for (i in 1:nrow(ven_data2)) {
  if (ven_data2$in_rbcl[i] == "yes") {
    molecular_genera_vec3 <- c(molecular_genera_vec3, ven_data2$Genus[i])
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
  x = list(mol = molecular_genera_vec2, rbcL = molecular_genera_vec3, Morphological = morpho_genera_vec2, Mark = marks_genera_vec),
  category.names = c("mol", "rbcL", "Morphological Count", "Mark's Count"),
  filename = NULL
)

#grid.draw(venn.plot2)

#### Abundance plot ####
# make a new dataframe and remove higher level classifications
plot_data <- plot_data %>%
  filter(!(Genus %in% c('Achnanthales', 'Bacillariales', 'Bacillariophyceae', 'Bacillariophytina', 'Bacillariophyta', 'Coscinodiscophytina', 
                        'Coscinodiscophyceae', 'Cymbellales', 'Diatomea', 'Fragilariales', 'Fragilariophyceae', 'Mediophyceae',
                        'Melosirids', 'Naviculales', 'Thalassiophysales', 'Thalassiosirales', 'Ochrophyta', 'Stauroneidaceae', 'Surirellales')))


# replace all of the zeros
plot_data <- plot_data %>%
  mutate(across(everything(), ~replace(., is.na(.), 0)))

# make new cols for color classification
plot_data$in18S <- ifelse(is.na(plot_data$in18S), "no", plot_data$in18S)
plot_data$inRBCL <- ifelse(is.na(plot_data$inRBCL), "no", plot_data$inRBCL)
plot_data$in_molecular <- ifelse(plot_data$in18S == "yes" | plot_data$inRBCL == "yes", "yes", "no")
plot_data$in_counts <- ifelse(is.na(plot_data$in_counts), "no", plot_data$in_counts)
#plot_data$in_counts <- ifelse(plot_data$in_counts == "excluded" | plot_data$in_counts == "yes", "yes", "no")

plot_data$where_found <- NA
plot_data$where_found <- ifelse(
  (plot_data$in_counts == "yes" | plot_data$in_counts == 'excluded') & plot_data$in_molecular == "yes" & plot_data$one_sample != 'yes', "in_both",
  ifelse(
    plot_data$in_counts == "yes" & plot_data$in_molecular == "no" , "in_morpho",
    ifelse(
      plot_data$in_counts == "no" & plot_data$in_molecular == "yes", "in_molecular",
      ifelse(
        plot_data$in_counts == "excluded" & plot_data$in_molecular == "no", "in_marks",
      ifelse(
        plot_data$one_sample == 'yes', 'in_one', NA
    )
  )
)
))


# correct the abundances
total_morphological_sum <- sum(plot_data$mean.sem.count, na.rm = TRUE)


# calculate per-sample relative abundance
plot_data$sem.ra <- plot_data$mean.sem.count / total_morphological_sum

ggplot(plot_data, aes(x = sem.ra, y = mean.mol.ra, color = where_found)) +
  geom_point()

# replace NA with zero
plot_data$sem.ra[is.na(plot_data$sem.ra)] <- 0
plot_data$sum.mol.ra[is.na(plot_data$sum.mol.ra)] <- 0


# # try with the ranks and diatoms not found as zero
# plot_data <- plot_data %>%
#   dplyr::mutate(across(everything(), ~replace(., is.na(.), 0)))


# plot_data$sem_ra_rank <- ifelse(plot_data$sum.sem.counts == 0, 0, rank(-plot_data$sem.ra, ties.method = "random"))
# plot_data$mol_ra_rank <- ifelse(plot_data$n.mol.samples == 0, 0, rank(-plot_data$sum.mol.ra, ties.method = "random"))

plot_data$sem_ra_rank <- ifelse(plot_data$sum.sem.counts == 0, 0, rank(-plot_data$sem.ra, ties.method = "min"))
plot_data$mol_ra_rank <- ifelse(plot_data$n.mol.samples == 0, 0, rank(-plot_data$sum.mol.ra, ties.method = "min"))

filtered_plot_data <- plot_data |>
  filter(where_found == "in_both" & sum.sem.counts > 0)

# redefine the basline
plot_data$mol_ra_rank[is.na(plot_data$mol_ra_rank) | plot_data$mol_ra_rank == 0] <- 60
plot_data$sem_ra_rank[is.na(plot_data$sem_ra_rank) | plot_data$sem_ra_rank == 0] <- 29

# Create the scatter plot with the filtered data and add the trendline
scatter_plot3 <- ggplot(plot_data, aes(x = sem_ra_rank, y = mol_ra_rank, color = where_found)) +
  geom_point(data = filtered_plot_data, size = 2.5, position = position_jitter(width = 0.3, height = 0)) +
  geom_smooth(data = filtered_plot_data, method = 'lm') +
  geom_point(data = plot_data[plot_data$sum.mol.ra == 0 | plot_data$sem.ra == 0 | plot_data$where_found == 'in_one' | plot_data$where_found == 'in_marks',], size = 2.5) +
  labs(x = "Morphological Rank", y = "Molecular Rank") +
  annotate("text", x = min(plot_data$sem_ra_rank) +5, y = max(plot_data$mol_ra_rank)-1,
           label = paste("p-value: 0.04", "\ny = 0.64x + 5.58"),
           hjust = 1, vjust = 1, size = 4) +
  scale_color_manual(labels = c("Both", "Molecular only", "Morphological census only", "Morphological survey only", "One Sample or Low Abundance"),
                     values = c('#46295d', '#82bac4', '#e37c78','lightpink', '#966c8b')) +
  guides(color=guide_legend(title="Found in Morphological\nor Molecular Data?")) 

scatter_plot3

lm <- lm(mol_ra_rank ~ sem_ra_rank, data = filtered_plot_data)
summary(lm) # 0.04

scatter_plot3 <- ggplot(plot_data, aes(x = sem_ra_rank, y = mol_ra_rank, color = where_found)) +
  geom_point(data = filtered_plot_data, size = 2.5) +
  geom_smooth(data = filtered_plot_data, method = 'lm') +
  geom_point(data = plot_data[plot_data$sum.mol.ra == 0 | plot_data$sem.ra == 0 | plot_data$where_found == 'in_one' | plot_data$where_found == 'in_marks',], size = 2.5) +
  labs(x = "Morphological Rank", y = "Molecular Rank") +
  annotate("text", x = min(plot_data$sem_ra_rank) +5, y = max(plot_data$mol_ra_rank)-1,
           label = paste("p-value: 0.04", "\ny = 0.64x + 5.58"),
           hjust = 1, vjust = 1, size = 4) +
  scale_color_manual(labels = c("Both", "Molecular only", "Morphological census only", "Morphological survey only", "One Sample or Low Abundance"),
                     values = c('#46295d', '#82bac4', '#e37c78','lightpink', '#966c8b')) +
  guides(color=guide_legend(title="Found in Morphological\nor Molecular Data?")) +
  geom_text(aes(label = Genus), vjust = -0.9, hjust = -0.1, size = 2, check_overlap = FALSE) 

scatter_plot3

#### separately ####

mol_RBCL$ra[is.na(mol_RBCL$ra)] <- 0
mol_RBCL <- mol_RBCL |>
  filter(!(Genus %in% c('Aulacodiscus', 'Bolidomonas', 'Biddulphiopsis', 'Caloneis', 'Cerataulina', 'Cymatosira',
                        'Cyclotella', 'Detonula', 'Diatoma', 'Discostella', 'Gomphonema', 'Guinardia', 'Lauderia',
                        'Pauliella', 'Placoneis', 'Porosira', 'Psammogramma', 'Rhoikoneis', 'Shionodiscus',
                        'Staurosira', 'Stephanodiscus', 'Hyalosira', 'Nanofrustulum', 'Pteridomonas')))


mol_RBCL = ddply(mol_RBCL, c("Genus"),
                 summarise,
                 sum.mol.ra = sum(ra),
                 mean.mol.ra = mean(ra),
                 n.mol.samples =length(unique(Row.names)))

mol_RBCL$inRBCL = 'yes'

in_one <- c('Astermophalus', 'Coscinodiscus', 'Dimeregramma', 'Ditylum', 'Ellerbeckia', 'Entomoneis', 'Fragilariopsis', 
            'Gomphonemopsis', 'Grammatophora', 'Leptocylindrus', 'Lyrella', 'Paralia', 'Plagiogramma',
            'Psammodictyon', 'Pseudo-nitzschia', 'Rhizosolenia', 'Trachyneis', 'Opephora', 'Fragilaria', 'Fogedia')

mol_RBCL$one_sample <- NA
mol_RBCL <- mol_RBCL %>%
  mutate(one_sample = ifelse(Genus %in% in_one, "yes", "no"))


###### combine data #####
all = full_join(mol_RBCL, diatom_data)

plot_data <- all %>%
  filter(!(Genus %in% c('Undatella', 'Actinoptychus', 'Rhopalodia',
                        'Hobaniella','Cyclotella', 'Donkinia')))

# add proskinia
new_values <- c(
  sum.sem.counts = 0,
  n.sem.samples = 17,
  mean.sem.count = 0,
  in_counts = "excluded"
)

# Find the row where Genus is 'Proschkinia'
row_index <- which(plot_data$Genus == "Proschkinia")

# Modify the values in that row
plot_data[row_index, c('sum.sem.counts', 'n.sem.samples', 'mean.sem.count', 'in_counts')] <- new_values

# collapse this
plot_data <- plot_data %>%
  mutate(Genus = ifelse(Genus == "Paribellus", "Parlibellus", Genus))

new_values <- tibble(
  Genus = "Parlibellus",
  sum.mol.ra = 9.181825e-02,
  mean.mol.ra =1.334568e-04,
  n.mol.samples = 43,
  inRBCL = 'yes',
  one_sample = 'no',
  sum.sem.counts = 0,
  n.sem.samples = 17,
  mean.sem.count = 0,
  in_counts = 'excluded'
  
)

# Filter out the original Parlibellus rows
filtered_data <- plot_data %>%
  filter(Genus != "Parlibellus")

# Combine the filtered data with the new values for Parlibellus
filtered_data$mean.sem.count = as.numeric(filtered_data$mean.sem.count)
filtered_data$sum.sem.counts = as.numeric(filtered_data$sum.sem.counts)
filtered_data$n.sem.samples = as.numeric(filtered_data$n.sem.samples)
plot_data <- bind_rows(filtered_data, new_values)


# fix Asteromphalus
plot_data <- plot_data |>
  mutate(in_counts = ifelse(Genus == 'Asteromphalus', 'excluded', in_counts))

plot_data <- plot_data |>
  filter(Genus != 'Paribellus')

# add new genera
new_values <- tibble(
  sum.sem.counts = 0,
  n.sem.samples = 17,
  mean.sem.count = 0,
  in_counts = 'excluded'
)


# Find the row where Genus is 'Proschkinia'
row_index <- which(plot_data$Genus == "Asteromphalus")

# Modify the values in that row
plot_data[row_index, c('sum.sem.counts', 'n.sem.samples', 'mean.sem.count', 'in_counts')] <- new_values

plot_data <- bind_rows(plot_data, new_values)

new_values <- tibble(
  Genus = "Eunotia",
  sum.mol.ra = 0,
  mean.mol.ra =0,
  n.mol.samples = 0,
  inRBCL = 'no',
  sum.sem.counts = 0,
  n.sem.samples = 17,
  mean.sem.count = 0,
  in_counts = 'excluded'
  
)

plot_data <- bind_rows(plot_data, new_values)

new_values <- tibble(
  Genus = "Fragilaria",
  sum.mol.ra = 0,
  mean.mol.ra =0,
  n.mol.samples = 0,
  inRBCL = 'no',
  sum.sem.counts = 0,
  n.sem.samples = 17,
  mean.sem.count = 0,
  in_counts = 'excluded'
  
)

plot_data <- bind_rows(plot_data, new_values)

plot_data <- plot_data %>%
  filter(!(Genus %in% c('Achnanthales', 'Bacillariales', 'Bacillariophyceae', 'Bacillariophytina', 'Bacillariophyta', 'Coscinodiscophytina', 
                        'Coscinodiscophyceae', 'Cymbellales', 'Diatomea', 'Fragilariales', 'Fragilariophyceae', 'Mediophyceae',
                        'Melosirids', 'Naviculales', 'Thalassiophysales', 'Thalassiosirales', 'Ochrophyta', 'Stauroneidaceae')))

plot_data <- plot_data |>
  filter(Genus != 0)

# make new cols for color classification
plot_data$inRBCL <- ifelse(is.na(plot_data$inRBCL), "no", plot_data$inRBCL)
plot_data$in_molecular <- ifelse(plot_data$inRBCL == "yes", "yes", "no")
plot_data$in_counts <- ifelse(is.na(plot_data$in_counts), "no", plot_data$in_counts)
#plot_data$in_counts <- ifelse(plot_data$in_counts == "excluded" | plot_data$in_counts == "yes", "yes", "no")

plot_data$where_found <- NA
plot_data$where_found <- ifelse(
  (plot_data$in_counts == "yes" | plot_data$in_counts == 'excluded') & plot_data$in_molecular == "yes" & plot_data$one_sample != 'yes', "in_both",
  ifelse(
    plot_data$in_counts == "yes" & plot_data$in_molecular == "no" , "in_morpho",
    ifelse(
      plot_data$in_counts == "no" & plot_data$in_molecular == "yes", "in_molecular",
      ifelse(
        plot_data$in_counts == "excluded" & plot_data$in_molecular == "no", "in_marks",
        ifelse(
          plot_data$one_sample == 'yes', 'in_one', NA
        )
      )
    )
  ))


# correct the abundances
total_morphological_sum <- sum(plot_data$mean.sem.count, na.rm = TRUE)


# calculate per-sample relative abundance
plot_data$sem.ra <- plot_data$mean.sem.count / total_morphological_sum

ggplot(plot_data, aes(x = sem.ra, y = mean.mol.ra, color = where_found)) +
  geom_point()

lm <- lm(mean.mol.ra ~ sem.ra, data = filtered_plot_data)
summary(lm) # 0.04

# replace NA with zero
plot_data$sem.ra[is.na(plot_data$sem.ra)] <- 0
plot_data$sum.mol.ra[is.na(plot_data$sum.mol.ra)] <- 0


# # try with the ranks and diatoms not found as zero
# plot_data <- plot_data %>%
#   dplyr::mutate(across(everything(), ~replace(., is.na(.), 0)))


plot_data$sem_ra_rank <- ifelse(plot_data$sum.sem.counts == 0, 0, rank(-plot_data$sem.ra, ties.method = "min"))
plot_data$mol_ra_rank <- ifelse(plot_data$n.mol.samples == 0, 0, rank(-plot_data$sum.mol.ra, ties.method = "min"))

filtered_plot_data <- plot_data |>
  filter(where_found == "in_both" & sum.sem.counts > 0)

# redefine the basline
plot_data$mol_ra_rank[is.na(plot_data$mol_ra_rank) | plot_data$mol_ra_rank == 0] <- 49
plot_data$sem_ra_rank[is.na(plot_data$sem_ra_rank) | plot_data$sem_ra_rank == 0] <- 32

# Create the scatter plot with the filtered data and add the trendline
scatter_plot3 <- ggplot(plot_data, aes(x = sem_ra_rank, y = mol_ra_rank, color = where_found)) +
  geom_point(data = filtered_plot_data, size = 2.5) +
  geom_smooth(data = filtered_plot_data, method = 'lm') +
  geom_point(data = plot_data[plot_data$sum.mol.ra == 0 | plot_data$sem.ra == 0 | plot_data$where_found == 'in_one' | plot_data$where_found == 'in_marks',], size = 2.5, position = position_jitter(width = 0.5, height = 0.5)) +
  labs(x = "Morphological Rank", y = "Molecular Rank") +
  annotate("text", x = min(plot_data$sem_ra_rank) +5, y = max(plot_data$mol_ra_rank)-10,
           label = paste("p-value: 0.16", "\nR-squared: 0.07"),
           hjust = 1, vjust = 1, size = 4) +
  scale_color_manual(labels = c("Both",  "Morphological survey only", "Molecular only", "Morphological census only", "One Sample or Low Abundance"),
                     values = c('#46295d', '#edc6c5','#82bac4', '#e37c78','#966c8b')) +
  guides(color=guide_legend(title="Found in Morphological\nor Molecular Data?")) 

scatter_plot3

lm <- lm(mol_ra_rank ~ sem_ra_rank, data = filtered_plot_data)
summary(lm) # 0.04

scatter_plot4 <- ggplot(plot_data, aes(x = sem_ra_rank, y = mol_ra_rank, color = where_found)) +
  geom_point(data = filtered_plot_data, size = 2.5) +
  geom_smooth(data = filtered_plot_data, method = 'lm') +
  geom_point(data = plot_data[plot_data$sum.mol.ra == 0 | plot_data$sem.ra == 0 | plot_data$where_found == 'in_one' | plot_data$where_found == 'in_marks',], size = 2.5) +
  labs(x = "Morphological Rank", y = "Molecular Rank") +
  annotate("text", x = min(plot_data$sem_ra_rank) +5, y = max(plot_data$mol_ra_rank)-1,
           label = paste("p-value: 0.04", "\ny = 0.64x + 5.58"),
           hjust = 1, vjust = 1, size = 4) +
  scale_color_manual(labels = c("Both", "Molecular only", "Morphological census only", "Morphological survey only", "One Sample or Low Abundance"),
                     values = c('#46295d', '#82bac4', '#e37c78','lightpink', '#966c8b')) +
  guides(color=guide_legend(title="Found in Morphological\nor Molecular Data?")) +
  geom_text(aes(label = Genus), vjust = -0.9, hjust = -0.1, size = 2, check_overlap = FALSE) 

scatter_plot4


#### 18S #####
mol_18S <- mol_18S |>
  filter(!(Genus %in% c('Aulacodiscus', 'Bolidomonas', 'Biddulphiopsis', 'Caloneis', 'Cerataulina', 'Cymatosira',
                        'Cyclotella', 'Detonula', 'Diatoma', 'Discostella', 'Gomphonema', 'Guinardia', 'Lauderia',
                        'Pauliella', 'Placoneis', 'Porosira', 'Psammogramma', 'Rhoikoneis', 'Shionodiscus',
                        'Staurosira', 'Stephanodiscus', 'Hyalosira', 'Nanofrustulum', 'Pteridomonas')))

mol_18S = ddply(mol_18S, c("Genus"),
                summarise,
                sum.mol.ra = sum(ra),
                mean.mol.ra = mean(ra),
                n.mol.samples =length(unique(Row.names)))

mol_18S$in18S = 'yes'

in_one <- c('Asteromphalus', 'Coscinodiscus', 'Dimeregramma', 'Ditylum', 'Ellerbeckia', 'Entomoneis', 'Fragilariopsis', 
            'Gomphonemopsis', 'Grammatophora', 'Leptocylindrus', 'Lyrella', 'Paralia', 'Plagiogramma',
            'Psammodictyon', 'Pseudo-nitzschia', 'Rhizosolenia', 'Trachyneis', 'Opephora', 'Fogedia', 'Fragilaria')

mol_18S$one_sample <- NA
mol_18S <- mol_18S %>%
  mutate(one_sample = ifelse(Genus %in% in_one, "yes", "no"))

###### combine data #####
all = full_join(mol_18S, diatom_data)

plot_data <- all %>%
  filter(!(Genus %in% c('Undatella', 'Actinoptychus', 'Rhopalodia',
                        'Hobaniella','Cyclotella', 'Donkinia')))

new_values <- tibble(
  Genus = "Proschkinia",
  sum.mol.ra = 0,
  mean.mol.ra =0,
  n.mol.samples = 0,
  in18S = 'no',
  sum.sem.counts = 0,
  n.sem.samples = 17,
  mean.sem.count = 0,
  in_counts = 'excluded'
  
)

plot_data <- bind_rows(plot_data, new_values)

# collapse this
plot_data <- plot_data %>%
  mutate(Genus = ifelse(Genus == "Paribellus", "Parlibellus", Genus))

new_values <- tibble(
  Genus = "Parlibellus",
  sum.mol.ra = 0,
  mean.mol.ra =0,
  n.mol.samples = NA,
  in18S = 'no',
  sum.sem.counts = 0,
  n.sem.samples = 17,
  mean.sem.count = 0,
  in_counts = 'excluded'
  
)

# Filter out the original Parlibellus rows
filtered_data <- plot_data %>%
  filter(Genus != "Parlibellus")

# Combine the filtered data with the new values for Parlibellus
filtered_data$mean.sem.count = as.numeric(filtered_data$mean.sem.count)
filtered_data$sum.sem.counts = as.numeric(filtered_data$sum.sem.counts)
filtered_data$n.sem.samples = as.numeric(filtered_data$n.sem.samples)
plot_data <- bind_rows(filtered_data, new_values)


# fix Asteromphalus
plot_data <- plot_data |>
  mutate(in_counts = ifelse(Genus == 'Asteromphalus', 'excluded', in_counts))

plot_data <- plot_data |>
  filter(Genus != 'Paribellus')

# add new genera
new_values <- tibble(
  Genus = "Asteromphalus",
  sum.mol.ra = 0,
  mean.mol.ra =0,
  n.mol.samples = 0,
  in18S = 'no',
  sum.sem.counts = 0,
  n.sem.samples = 17,
  mean.sem.count = 0,
  in_counts = 'excluded'
  
)

plot_data <- bind_rows(plot_data, new_values)

new_values <- tibble(
  Genus = "Eunotia",
  sum.mol.ra = 0,
  mean.mol.ra =0,
  n.mol.samples = 0,
  in18S = 'no',
  sum.sem.counts = 0,
  n.sem.samples = 17,
  mean.sem.count = 0,
  in_counts = 'excluded'
  
)

plot_data <- bind_rows(plot_data, new_values)

plot_data <- plot_data |>
  filter(Genus != 'Fragilaria')

new_values <- tibble(
  Genus = "Fragilaria",
  sum.mol.ra = 0,
  mean.mol.ra =0,
  n.mol.samples = 0,
  in18S = 'no',
  sum.sem.counts = 0,
  n.sem.samples = 17,
  mean.sem.count = 0,
  in_counts = 'excluded'
  
)

plot_data <- bind_rows(plot_data, new_values)

plot_data <- plot_data %>%
  filter(!(Genus %in% c('Achnanthales', 'Bacillariales', 'Bacillariophyceae', 'Bacillariophytina', 'Bacillariophyta', 'Coscinodiscophytina', 
                        'Coscinodiscophyceae', 'Cymbellales', 'Diatomea', 'Fragilariales', 'Fragilariophyceae', 'Mediophyceae',
                        'Melosirids', 'Naviculales', 'Thalassiophysales', 'Thalassiosirales', 'Ochrophyta', 'Stauroneidaceae')))

# make new cols for color classification
plot_data$in18S <- ifelse(is.na(plot_data$in18S), "no", plot_data$in18S)
plot_data$in_molecular <- ifelse(plot_data$in18S == "yes", "yes", "no")
plot_data$in_counts <- ifelse(is.na(plot_data$in_counts), "no", plot_data$in_counts)
#plot_data$in_counts <- ifelse(plot_data$in_counts == "excluded" | plot_data$in_counts == "yes", "yes", "no")


plot_data$where_found <- NA
plot_data$where_found <- ifelse(
  (plot_data$in_counts == "yes" | plot_data$in_counts == 'excluded') & plot_data$in_molecular == "yes" & plot_data$one_sample != 'yes', "in_both",
  ifelse(
    plot_data$in_counts == "yes" & plot_data$in_molecular == "no" , "in_morpho",
    ifelse(
      plot_data$in_counts == "no" & plot_data$in_molecular == "yes", "in_molecular",
      ifelse(
        plot_data$in_counts == "excluded" & plot_data$in_molecular == "no", "in_marks",
        ifelse(
          plot_data$one_sample == 'yes', 'in_one', NA
        )
      )
    )
  ))


# correct the abundances
total_morphological_sum <- sum(plot_data$mean.sem.count, na.rm = TRUE)


# calculate per-sample relative abundance
plot_data$sem.ra <- plot_data$mean.sem.count / total_morphological_sum

ggplot(plot_data, aes(x = sem.ra, y = mean.mol.ra, color = where_found)) +
  geom_point()

# replace NA with zero
plot_data$sem.ra[is.na(plot_data$sem.ra)] <- 0
plot_data$sum.mol.ra[is.na(plot_data$sum.mol.ra)] <- 0


# # try with the ranks and diatoms not found as zero
# plot_data <- plot_data %>%
#   dplyr::mutate(across(everything(), ~replace(., is.na(.), 0)))


plot_data$sem_ra_rank <- ifelse(plot_data$sum.sem.counts == 0, 0, rank(-plot_data$sem.ra, ties.method = "min"))
plot_data$mol_ra_rank <- ifelse(plot_data$n.mol.samples == 0, 0, rank(-plot_data$sum.mol.ra, ties.method = "min"))

filtered_plot_data <- plot_data |>
  filter(where_found == "in_both" & sum.sem.counts > 0)

# redefine the basline
plot_data$mol_ra_rank[is.na(plot_data$mol_ra_rank) | plot_data$mol_ra_rank == 0] <- 40
plot_data$sem_ra_rank[is.na(plot_data$sem_ra_rank) | plot_data$sem_ra_rank == 0] <- 32

# Create the scatter plot with the filtered data and add the trendline
scatter_plot5 <- ggplot(plot_data, aes(x = sem_ra_rank, y = mol_ra_rank, color = where_found)) +
  geom_point(data = filtered_plot_data, size = 2.5) +
  geom_smooth(data = filtered_plot_data, method = 'lm') +
  geom_point(data = plot_data[plot_data$sum.mol.ra == 0 | plot_data$sem.ra == 0 | plot_data$where_found == 'in_one' | plot_data$where_found == 'in_marks',], size = 2.5, position = position_jitter(width = 0.5, height = 0.5)) +
  labs(x = "Morphological Rank", y = "Molecular Rank") +
  annotate("text", x = min(plot_data$sem_ra_rank) +5, y = max(plot_data$mol_ra_rank)-10,
           label = paste("p-value: 0.002", "\nR-squared: 0.43"),
           hjust = 1, vjust = 1, size = 4) +
  scale_color_manual(labels = c("Both",  "Morphological survey only", "Molecular only", "Morphological census only", "One Sample or Low Abundance"),
                     values = c('#46295d', '#edc6c5','#82bac4', '#e37c78','#966c8b')) +
  guides(color=guide_legend(title="Found in Morphological\nor Molecular Data?")) 

scatter_plot5

lm <- lm(mol_ra_rank ~ sem_ra_rank, data = filtered_plot_data)
summary(lm) # 0.04

scatter_plot6 <- ggplot(plot_data, aes(x = sem_ra_rank, y = mol_ra_rank, color = where_found)) +
  geom_point(data = filtered_plot_data, size = 2.5) +
  geom_smooth(data = filtered_plot_data, method = 'lm') +
  geom_point(data = plot_data[plot_data$sum.mol.ra == 0 | plot_data$sem.ra == 0 | plot_data$where_found == 'in_one' | plot_data$where_found == 'in_marks',], size = 2.5) +
  labs(x = "Morphological Rank", y = "Molecular Rank") +
  annotate("text", x = min(plot_data$sem_ra_rank) +5, y = max(plot_data$mol_ra_rank)-1,
           label = paste("p-value: 0.04", "\ny = 0.64x + 5.58"),
           hjust = 1, vjust = 1, size = 4) +
  scale_color_manual(labels = c("Both", "Molecular only", "Morphological census only", "Morphological survey only", "One Sample or Low Abundance"),
                     values = c('#46295d', '#82bac4', '#e37c78','lightpink', '#966c8b')) +
  guides(color=guide_legend(title="Found in Morphological\nor Molecular Data?")) +
  geom_text(aes(label = Genus), vjust = -0.9, hjust = -0.1, size = 2, check_overlap = FALSE) 

scatter_plot6


ggarrange(scatter_plot5, scatter_plot3, labels = c("A", "B"), ncol = 2)

ggarrange(scatter_plot6, scatter_plot4, labels = c("A", "B"), ncol = 2)
