# A. Jackman
# Jan 21 2024
# check whether the patchy diatoms (colonies)
# 1) are present in all of the leaves vs some of the leaves?
# 2) are present in all of the sections from the same leaves?
# Plan:
# a. identify the patchy diatoms
# b. find the proportion of boxes the diatoms are present in
# c. find the proportion of sections from the same leaves they are present in
# d. find the proportion of 

# libraries
library(dplyr)
library(tidyverse)

# set working directory
setwd("/Users/andreajackman/R_Stuff/diatoms_new/counting_data")

# read in the data
diatom_sums <- read.csv("diatom_sum_per_box.csv")

# select cols to keep
cols_to_keep <- c('box', 'leaf', 'blade_position', 'gomphonemopsis', 'rhoicosphenia', 'pseudogonphonema')

# subset for the selected diatoms
diatom_sums <- diatom_sums |>
  select(all_of(cols_to_keep))

# check the proportion of leaves they are present in
present <- apply(diatom_sums[, patchy_diatom], 2, function(x) any(x != 0))
present # present in all leaves

# Calculate the number of boxes with each species
gomphonemopsis_boxes <- sum(diatom_sums$gomphonemopsis > 0)
rhoicosphenia_boxes <- sum(diatom_sums$rhoicosphenia > 0)
pseudogonphonema_boxes <- sum(diatom_sums$pseudogonphonema > 0)

# Calculate the proportions
proportion_gomphonemopsis <- gomphonemopsis_boxes / 15 # 0.7333333
proportion_rhoicosphenia <- rhoicosphenia_boxes / 15 # 0.7333333
proportion_pseudogomphonema <- pseudogonphonema_boxes / 15 # 0.7333333

# the proportion from each leaf
diatom_sums %>%
  group_by(leaf) %>%
  summarize(
    gomphonemopsis_boxes = sum(gomphonemopsis > 0) / n(),
    rhoicosphenia_boxes = sum(rhoicosphenia > 0) / n(),
    pseudogonphonema_boxes = sum(pseudogonphonema > 0 )/ n())

# leaf  gomphonemopsis_boxes rhoicosphenia_boxes pseudogonphonema_boxes
# <chr>                <dbl>               <dbl>                  <dbl>
#  0                    0.714               0.714                  0.714
#  12A                  0.6                 0.6                    0.6  
#  9                    1                   1                      1  







