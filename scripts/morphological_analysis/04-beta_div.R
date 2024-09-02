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
setwd("/Users/andreajackman/R_Stuff/Parfrey/MH_diatoms_metadata_and_scripts/morphological_data")

## read in data
diatom_sums <- read.csv("diatom_sum_per_box.csv")

# move the last three columns to positions 3, 4, and 5
diatom_sums <- diatom_sums[, c(1:2, 27:29, 3:26, 29:ncol(diatom_sums))]
diatom_sums <- diatom_sums[, -c(30)]

##### PERMANOVA #####
## PERMANOVA ##
set.seed(2020)
## run permanova with marginal effects 
m.perm= adonis2(diatom_sums[,-c(1:5)] ~ blade_position, 
                data=diatom_sums, method = "bray")
m.perm
# Df SumOfSqs      R2      F Pr(>F)   
# blade_position  2  0.95968 0.34461 3.1548  0.004 **
#   Residual       12  1.82516 0.65539                 
# Total          14  2.78484 1.00000        

## POST HOC TEST ##
## post hoc test
m.pairwise = pairwise.adonis(diatom_sums[,-c(1:5)], diatom_sums$blade_position)
m.pairwise
# pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
# 1 B vs M  1 0.4901218 2.551507 0.2671314   0.080      0.240    
# 2 B vs T  1 0.6976320 3.748409 0.2940296   0.009      0.027   .
# 3 M vs T  1 0.2275071 2.901069 0.2661270   0.060      0.180    

## BETADISPERSON TEST ##
# calculate the distance
dist <- vegdist(diatom_sums[, 6:ncol(diatom_sums)], method = "bray")
# groups
groups <- diatom_sums$blade_position
# calculate the betadispersion within each region
dispers <- betadisp_result <- betadisper(dist, groups)
# betadispersion test to see if all regions have the same betadispersion
beta=permutest(dispers) 
beta
# Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)  
# Groups     2 0.24983 0.124915 6.1208    999  0.015 *
# Residuals 12 0.24490 0.020408  

# post-hoc test
# pairwise_beta <- pairwise.perm.manova(dist, groups)
# pairwise_beta

##### NMDS PLOT #####
# stress plot
stressplot(metaMDS(dist)) # there is an elbow around k = 2

nmds_result <- metaMDS(dist, k = 2)

# Add NMDS results to your data frame
diatom_sums$NMDS1 <- nmds_result$points[, 1]
diatom_sums$NMDS2 <- nmds_result$points[, 2]

bray <- ggplot(diatom_sums, aes(x = NMDS1, y = NMDS2, color = blade_position, shape = leaf)) +
  geom_point(size = 4) +
  scale_color_manual(values=c("#55D3A6", "#5582D3", "#A655D3"))+
  annotate("text", x = -1.9, y = 1.3, label = "stress = 0.02499072 (Bray-Curtis)", size=4, color="black")+
  annotate("text", x = -1.9, y = 1.2, label = "PERMANOVA: Section F= 3.15, p= 0.004", size=4, color="black")+
  annotate("text", x = -1.9, y = 1.1, label = "No sd", size=4, color="black")+
  annotate("text", x = -1.9, y = 1, label = "BETADISPERSON F= 6.12, p= 0.015", size=4, color="black") +
  #theme_minimal() +
  labs(color = "Blade Position", shape = "Leaf") 

bray