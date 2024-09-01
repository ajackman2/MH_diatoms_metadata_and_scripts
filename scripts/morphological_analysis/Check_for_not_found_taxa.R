# libraries
library(dplyr)
library(plyr)
library(tidyverse)
library(tidyr)
library(stringr)

diat_barcode <- read.csv("imput_data/v10_correspondence.csv")

# filter for Rhoicospehnia
only_rhoico <- diat_barcode |>
  filter(genus == "Rhoicosphenia")

unique(only_rhoico$sequence) # 3

# filter for Achnanthes
only_achnan <- diat_barcode |>
  filter(genus == "Achnanthes")

unique(only_achnan$sequence) # 3

# filter for Eucampia
only_Eucampia <- diat_barcode |>
  filter(genus == "Eucampia")

unique(only_Eucampia$sequence) # 3

# filter for Pseudo-nitzschia
only_Pseudo_nitzschia <- diat_barcode |>
  filter(genus == "Pseudo-nitzschia")

unique(only_Pseudo_nitzschia$sequence) # 3
