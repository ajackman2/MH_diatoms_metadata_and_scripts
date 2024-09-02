###### SET UP #######
## load packages
library(tidyverse)
library(dplyr)
library(qualpalr)
library(ggpubr)
library(ggh4x)
library(indicspecies)
library(data.table)
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
setwd("/Users/andreajackman/R_Stuff/Parfrey/MH_diatoms_metadata_and_scripts/morphological_data")

## read in data
diatom_sums <- read.csv("diatom_sum_per_box.csv")
species_num_w_name <- read.csv("species_num_w_names.csv")
species_num <- read.csv("species_num.csv")

## indval comparison 
substrate <- as.character(species_num_w_name$blade_position) 

## run indval
indval <- multipatt(species_num, substrate, duleg = TRUE, control = how(nperm=999))

# Get indval statistic
indval.16.str <- as.data.frame(indval$str)
indval.16.str$rn <- rownames(indval.16.str)

# get p-value
indval.stat <- as.data.frame(indval$sign) #get dataframe of indval statistic
indval.stat$rn <- rownames(indval.stat) # make column of ASVs

# Prevalence as dataframe
indval.prev <- as.data.frame(indval$A)
# extract rownames into column
setDT(indval.prev, keep.rownames = TRUE)[] 
## rename columns 
colnames(indval.prev) <- paste0("prev.", colnames(indval.prev))
names(indval.prev)[names(indval.prev) == 'prev.rn'] <- 'rn'

# Fidelity as dataframe
indval.fid <- as.data.frame(indval$B) 
# extract rownames into column
setDT(indval.fid, keep.rownames = TRUE)[] 
## rename columns 
colnames(indval.fid) <- paste0("fid.", colnames(indval.fid))
names(indval.fid)[names(indval.fid) == 'fid.rn'] <- 'rn'

# Join statistics together (you could do multi_join but this might crash your computer)
str.and.stat = full_join(indval.16.str, indval.stat,
                         by="rn")
prev.and.fid = full_join(indval.prev, indval.fid,
                         by="rn")
indval_table = full_join(str.and.stat, prev.and.fid,
                         by="rn")

# remove any rows with NAs - these diatoms are missing because I reduced the sample size
indval_table <- na.omit(indval_table)

# rename the column
names(indval_table)[names(indval_table) == 'rn'] <- 'genera'

coresig = subset(indval_table, indval_table$index != "NaN" &
                   indval_table$p.value<=0.05 &
                   indval_table$stat>=0.5)

#coresig = subset(indval_table, indval_table$index != "NaN")

###### PLOT RA OF DIFF TAXA ######

# restructure diatom_data
diatom_sums <- diatom_sums |>
  pivot_longer(Navicula:Cocconeis, names_to = 'genera', values_to = 'abundance')

## only keep sig taxa
diatom_sums = inner_join(diatom_sums, coresig)

## relative abundance
diatom_sums$ra = as.numeric(diatom_sums$abundance)/as.numeric(diatom_sums$total)

## remove 0s to haveblanks
diatom_sums$ra = replace(diatom_sums$ra, diatom_sums$ra == 0, NA)

## replace index with letters
diatom_sums$index = gsub("1", "Base", diatom_sums$index)
diatom_sums$index = gsub("2", "Mid", diatom_sums$index)
diatom_sums$index = gsub("3", "Tip", diatom_sums$index)

## make bubble plot 
ggplot(diatom_sums, aes(x=as.character(box), 
                        y=paste0(genera), 
                        size=ra, color=index))+
  geom_point()+
  facet_grid(genera~ blade_position, space="free", scale="free")+
  theme(strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x="Leaf", y = "Genus", 
       size="Relative Abundance", color="Enriched On:")

#ggsave("indval_RBCL.pdf", width = 16, height=7, units = "in")
