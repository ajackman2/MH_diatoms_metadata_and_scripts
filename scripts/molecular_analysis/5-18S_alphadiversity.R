###### SET UP #######
## load packages
library(phyloseq)
library(tidyverse)
library(ggpubr)
library(ggh4x)
library(vegan)
library(car)
library(ggplot2); theme_set(theme_bw()+
                              theme(panel.grid = element_blank(),
                                    strip.background = element_rect(fill="white"),
                                    axis.text = element_text(size = 10, colour = "black"),
                                    axis.title = element_text(size=10, face="bold"),
                                    strip.text = element_text(color="black", size=10),
                                    legend.text=element_text(size=10),
                                    axis.line = element_line(colour = "black")))

## tell R whereto get data
setwd("C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Seagrass_Montague_2021-03-07/galliano_seagrass_2021and2023/imput_data")
path = "C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Seagrass_Montague_2021-03-07/galliano_seagrass_2021and2023/output/betadiv/"


## read in data
rare = readRDS("MH_phylo18_coverage_rarefied_2021to2023.RDS")

rare = subset_samples(rare, leaf_section !='NA')
rare
##### FORMAT FOR AND GET ALPHA-STATS######
meta = as.data.frame(rare@sam_data)
otu = as.data.frame(as.matrix(rare@otu_table))

meta$rich = specnumber(otu)
meta$shan = diversity(otu, index="shannon")
meta$simp = diversity(otu, index="simpson")

##### RUN ANOVA #####
arich = aov(meta$rich~meta$leaf_section)
summary(arich)
leveneTest(arich)
#plot(arich)

ashan = aov(meta$shan~meta$leaf_section)
summary(ashan)
leveneTest(ashan)
#plot(ashan)

asimp = aov(meta$simp~meta$leaf_section)
summary(asimp)
leveneTest(asimp)
#plot(asimp)


###### PLOT ALPHA DIV #####
r=ggplot(meta, aes(x=leaf_section, y=rich))+
  geom_boxplot()+
  annotate("text", x = 2, y = 50, label = "ns", size=3, color="black")

sh=ggplot(meta, aes(x=leaf_section, y=shan))+
  geom_boxplot()+
  annotate("text", x = 2, y = 0.02, label = "ns", size=3, color="black")

si=ggplot(meta, aes(x=leaf_section, y=simp))+
  geom_boxplot()+
  annotate("text", x = 2, y = 0.02, label = "ns", size=3, color="black")

ggarrange(r, sh, si, ncol=3, nrow=1)

