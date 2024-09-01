###### SET UP #######
## load packages
library(phyloseq)
library(tidyverse)
library(ggpubr)
library(ggh4x)
library(vegan)
library(pairwiseAdonis)
library(ggplot2); theme_set(theme_bw()+
                              theme(panel.grid = element_blank(),
                                    strip.background = element_rect(fill="white"),
                                    axis.text = element_text(size = 10, colour = "black"),
                                    axis.title = element_text(size=10, face="bold"),
                                    strip.text = element_text(color="black", size=10),
                                    legend.text=element_text(size=10),
                                    axis.line = element_line(colour = "black")))

## tell R whereto get data
path = "C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Seagrass_Montague_2021-03-07/galliano_seagrass_2021and2023/output/betadiv/"

## read in data
rare = readRDS("C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Seagrass_Montague_2021-03-07/galliano_seagrass_2021and2023/imput_data/MH_phylo16_coverage_rarefied_2021to2023.RDS")
rare
rare = subset_samples(rare, LeafSection !='NA')
rare

##### PERMANOVA #####
## PERMANOVA ##
# get data frames
motu = as.data.frame(as.matrix(rare@otu_table))
mmeta = as.data.frame(as.matrix(rare@sam_data))
## merge the data
mm = merge(mmeta,motu, by=0)
## count colum numbers
metacols = ncol(mmeta)+1
## run permanova 
m.perm= adonis2(mm[,-c(1:metacols)] ~ LeafSection, 
                  data=mm, method = "bray")
m.perm

## POST HOC TEST ##
## post hoc test
m.pairwise = pairwise.adonis(mm[,-c(1:metacols)], mm$LeafSection)
m.pairwise

## BETADISPERSON TEST ##
## calculate the distance within your data using the same ordination method as you do for your PERMANOVA
dist <- phyloseq::distance(rare, method = "bray")
## get your metadata out from phyloseq
sample_df <- data.frame(sample_data(rare))
## calculate the betadispersion within each region
dispers <- betadisper(dist, sample_df$LeafSection) 
## betadispersion test to see if all regions have the same betadispersion
beta=permutest(dispers) 
beta


##### NMDS PLOT #####

## ordinate your data to make the NMDS plot 
ord = ordinate(rare, "NMDS","bray")
ord ## stress = 0.1218895 

## make the NMDS plot
plot_ordination(rare, ord, color="LeafSection")+
  geom_point(size=3)+
  scale_color_manual(values=c("#55D3A6", "#5582D3", "#A655D3"))+
  annotate("text", x=-0.4, y=-0.5, label="16S", size=3, color="black")+
  annotate("text", x = -0.4, y = -0.6, label = "stress =0.1236143  (Bray-Curtis)", size=3, color="black")+
  annotate("text", x = -0.4, y = -0.7, label = "PERMANOVA:Section(F=1.7292, p =0.037)", size=3, color="black")+
  annotate("text", x = -0.4, y = -0.8, label = "No sd at p<=0.05", size=3, color="black")+
  annotate("text", x = -0.4, y = -0.9, label = "BETADISPERSON F=1.0259, p=0.367", size=3, color="black")



ggsave(filename = ("nmds_plots_16S.pdf"), 
       width=7.9, height=6, units="in", 
       path=path)





