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
rare = readRDS("MH_16S_2021_and_2022_coveraged_rarified001percent.RDS")

rare = subset_samples(rare, LeafSection !='NA')

##### FORMAT FOR AND GET ALPHA-STATS######
meta = as.data.frame(rare@sam_data)
otu = as.data.frame(t(as.matrix(rare@otu_table)))

meta$rich = specnumber(otu)
meta$shan = diversity(otu, index="shannon")
meta$simp = diversity(otu, index="simpson")

##### RUN ANOVA #####
arich = aov(meta$rich~meta$LeafSection)
summary(arich)
leveneTest(arich)
#plot(arich)

ashan = aov(meta$shan~meta$LeafSection)
summary(ashan)
#leveneTest(ashan)
#plot(ashan)
#TukeyHSD(ashan)
ashan = oneway.test(meta$shan~meta$LeafSection, var.equal = F)
ashan

asimp = aov(meta$simp~meta$LeafSection)
summary(asimp)
leveneTest(asimp)
#plot(asimp)
## welsh anova
asimp = oneway.test(meta$simp~meta$LeafSection, var.equal = F)
asimp

###### PLOT ALPHA DIV #####
r=ggplot(meta, aes(x=LeafSection, y=rich))+
  geom_boxplot()+
  annotate("text", x = 2, y = 50, label = "F=1.4586, p=0.24", size=3, color="black")

sh=ggplot(meta, aes(x=LeafSection, y=shan))+
  geom_boxplot()+
  annotate("text", x = 2, y = 0.02, label = "F=5.579, p=0.02 (base sd)", size=3, color="black")

si=ggplot(meta, aes(x=LeafSection, y=simp))+
  geom_boxplot()+
  annotate("text", x = 2, y = 0.02, label = "Welsh F=2.0206, p= 0.172", size=3, color="black")

ggarrange(r, sh, si, ncol=3, nrow=1)
