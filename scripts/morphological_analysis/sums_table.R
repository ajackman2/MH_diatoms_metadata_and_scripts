# table of most common diatoms
# A. Jackman
# Nov 28 2023

# libraries
library(dplyr)
library(tidyverse)
library(plyr)
library(tidyr)
library(ggplot2); theme_set(theme_bw()+
                              theme(panel.grid = element_blank(),
                                    strip.background = element_rect(fill="white"),
                                    axis.text.y = element_text(size = 12, colour = "black"),
                                    axis.title = element_text(size=15, face="bold"),
                                    strip.text = element_text(color="black", size=10),
                                    legend.text=element_text(size=10),
                                    axis.line = element_line(colour = "black"),
                                    axis.text.x = element_blank(),))


## set working directory
setwd("/Users/andreajackman/R_Stuff/diatoms_new/counting_data")


# read in the data
diatom_data <- read.csv("master_diatom_data.csv")

diatom_sum = diatom_data |>
  select(file_name, box, blade_position, Navicula:Undatella)

genus_columns <- diatom_sum[, c( "Navicula", "Tabularia", "Cylindrotheca", "Pseudonitzschia",
                                 "Gomphonemopsis", "Hyalodiscus", "Coscinodiscus", "Plagiotropis",
                                 "Nitzschia", "Rhoicosphenia", "Pseudogomphonema", "Thalassiosira",
                                 "Eucampia", "Chaetoceros", "Diploneis", "Entomoneis", "Pariraphis",
                                 "Planothidium", "Halamphora", "Achnanthes", "Amphora", "Petroneis",
                                 "Bacillaria", "Ellerbeckia", "Minidiscus", "Fragilariopsis", "Rhizosolenia",
                                 "Fogedia", "Skeletonema", "Pleurosigma", "Gomphoseptatum", "Licmophora",
                                 "Cocconeis", "Actinoptychus", "Attheya", "Cyclotella", "Dimeregramma",
                                 "Ditylum", "Donkinia", "Epithemia", "Eupyxidicula", "Encyonema", "Fallacia",
                                 "Grammatophora", "Gyrosigma", "Hanzschia", "Haslea", "Hobaniella",
                                 "Leptocylindrus", "Lyrella", "Melosira", "Odontella", "Opephora", "Paralia",
                                 "Paribellus", "Plagiogramma", "Podosira", "Psammodictyon", "Rhabdonema",
                                 "Rhopalodia", "Trigonium", "Trachyneis", "Tryblionella", "Undatella")]

sum_per_genus <- colSums(genus_columns, na.rm = TRUE)

# make a dataframe
total_counts_table <- data.frame(
  genus = names(sum_per_genus),
  sum_value = sum_per_genus
)

#View(total_counts_table)

# now separate it by the blade position
counts_position <- diatom_data |>
  select(file_name, box, blade_position, Navicula:Undatella) |>
  group_by(blade_position) |>
  summarise_at(vars(Navicula:Undatella), ~sum(., na.rm = TRUE))
#View(counts_position)

# make a new table
counts_pos_table <- counts_position |>
  pivot_longer(cols = Navicula:Undatella, names_to = "genus", values_to = "sum")

counts_pos_table <- counts_pos_table |>
  pivot_wider(names_from = blade_position, values_from = sum)
#View(counts_pos_table)

# pivot longer
counts_position_long <- diatom_data |>
  select(file_name, box, blade_position, cocconeis_small:licmophora) |>
  group_by(box, blade_position) |>
  summarise_at(vars(cocconeis_small:licmophora), ~sum(., na.rm = TRUE)) |>
  pivot_longer(cols = cocconeis_small:licmophora, names_to = "genus", values_to = "sum")
#View(counts_position_long)
