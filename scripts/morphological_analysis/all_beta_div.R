# combination of all beta div graphs

##### 16S #####
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

set.seed(15)
## read in data
rare = readRDS("imput_data/MH_phylo16_coverage_rarefied_2021to2023.RDS")
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
plot_16S <- plot_ordination(rare, ord, color="LeafSection")+
  geom_point(size=3)+
  scale_color_manual(values=c("#55D3A6", "#5582D3", "#A655D3")) +
  labs(color = "Blade Section") +
  annotate("text", x=-0.5, y=1.1, label="16S", size=3, color="black")+
  annotate("text", x = -0.5, y = 1, label = "stress =0.1236143  (Bray-Curtis)", size=3, color="black")+
  annotate("text", x = -0.5, y = 0.9, label = "PERMANOVA:Section(F=1.7292, p =0.037)", size=3, color="black")+
  annotate("text", x = -0.5, y = 0.8, label = "No sd at p<=0.05", size=3, color="black")+
  annotate("text", x = -0.5, y = 0.7, label = "BETADISPERSON F=1.0259, p=0.367", size=3, color="black")

#### 18S ####

## read in data
rare = readRDS("imput_data/MH_phylo18_coverage_rarefied_2021to2023.RDS")

rare = subset_samples(rare, leaf_section !='NA')

##### PERMANOVA #####
## PERMANOVA ##
# get data frames
motu = as.data.frame(as.matrix(rare@otu_table))
mmeta = as.data.frame(as.matrix(rare@sam_data))
## merge the data
mm = merge(mmeta,motu, by=0)
## count colum numbers
metacols = ncol(mmeta)+1
## run permanova with marginal effects 
m.perm= adonis2(mm[,-c(1:metacols)] ~ leaf_section, 
                data=mm, method = "bray")
m.perm

## POST HOC TEST ##
## post hoc test
m.pairwise = pairwise.adonis(mm[,-c(1:metacols)], mm$leaf_section)
m.pairwise

## BETADISPERSON TEST ##
## calculate the distance within your data using the same ordination method as you do for your PERMANOVA
dist <- phyloseq::distance(rare, method = "bray")
## get your metadata out from phyloseq
sample_df <- data.frame(sample_data(rare))
## calculate the betadispersion within each region
dispers <- betadisper(dist, sample_df$leaf_section) 
## betadispersion test to see if all regions have the same betadispersion
beta=permutest(dispers) 
beta


##### NMDS PLOT #####

## ordinate your data to make the NMDS plot 
ord = ordinate(rare, "NMDS","bray")
ord ## stress = 0.1457114 

rare@sam_data$leaf_number = as.character(rare@sam_data$leaf_number)

## make the NMDS plot
plot_18S <- plot_ordination(rare, ord, color="leaf_section")+
  geom_point(size=3)+
  scale_color_manual(values=c("#55D3A6", "#5582D3", "#A655D3")) +
  labs(color = "Blade Section") +
  annotate("text", x=-0.8, y=1.1, label="18S", size=3, color="black")+
  annotate("text", x = -0.8, y = 1, label = "stress = 0.1457114 (Bray-Curtis)", size=3, color="black")+
  annotate("text", x = -0.8, y = 0.9, label = "PERMANOVA: Section F=1.0728,  p=0.323", size=3, color="black")+
  annotate("text", x = -0.8, y = 0.8, label = "No sd", size=3, color="black")+
  annotate("text", x = -0.8, y = 0.7, label = "BETADISPERSON F=0.0897, p=0.91", size=3, color="black")

#### rbcl ####
## read in data
rare = readRDS("imput_data/MH_phyloRBCL_coverage_rarefied_2021to2023.RDS")
rare = subset_samples(rare, leafsection !='NA')

##### PERMANOVA #####
## PERMANOVA ##
# get data frames
motu = as.data.frame(as.matrix(rare@otu_table))
mmeta = as.data.frame(as.matrix(rare@sam_data))
## merge the data
mm = merge(mmeta,motu, by=0)
## count colum numbers
metacols = ncol(mmeta)+1
## run permanova with marginal effects 
m.perm= adonis2(mm[,-c(1:metacols)] ~ leafsection, 
                data=mm, method = "bray")
m.perm

## POST HOC TEST ##
## post hoc test
m.pairwise = pairwise.adonis(mm[,-c(1:metacols)], mm$leafsection)
m.pairwise

## BETADISPERSON TEST ##
## calculate the distance within your data using the same ordination method as you do for your PERMANOVA
dist <- phyloseq::distance(rare, method = "bray")
## get your metadata out from phyloseq
sample_df <- data.frame(sample_data(rare))
## calculate the betadispersion within each region
dispers <- betadisper(dist, sample_df$leafsection) 
## betadispersion test to see if all regions have the same betadispersion
beta=permutest(dispers) 
beta


##### NMDS PLOT #####

## ordinate your data to make the NMDS plot 
ord = ordinate(rare, "NMDS","bray")
ord ## stress = 0.1815948 

rare@sam_data$leafnumber = as.character(rare@sam_data$leafnumber)


## make the NMDS plot
plot_rbcl <- plot_ordination(rare, ord, color="leafsection")+
  geom_point(size=3)+
  labs(color = 'Blade Section') +
  scale_color_manual(values=c("#55D3A6", "#5582D3", "#A655D3")) +
  annotate("text", x=-0.95, y=1, label="rbcL", size=3, color="black")+
  annotate("text", x = -0.95, y = 0.9, label = "stress = 0.1865233 (Bray_Curtis)", size=3, color="black")+
  annotate("text", x = -0.95, y = 0.8, label = "PERMANOVA: Section F=1.1261,  p=0.32", size=3, color="black")+
  annotate("text", x = -0.95, y = 0.7, label = "No sd", size=3, color="black")+
  annotate("text", x = -0.95, y = 0.6, label = "BETADISPERSON F=0.0762, p=0.921", size=3, color="black")
  

#### morphological ####
## tell R whereto get data
#setwd("/Users/andreajackman/R_Stuff/Parfrey/diatoms_new/counting_data")

## read in data
diatom_sums <- read.csv("counting_data/diatom_sum_per_box.csv")

# add the plant numer
# diatom_sums$plant <- NA
# diatom_sums$plant <- ifelse(diatom_sum$box == "0_19_7", "7",
#                            ifelse(diatom_sum$box == "0_1_1Pa", "1",
#                                   ifelse(diatom_sum$box == "0_20_7", "7",
#                                          ifelse(diatom_sum$box == "0_21_7", "7",
#                                                 ifelse(diatom_sum$box == "0_3_1Da", "1",
#                                                        ifelse(diatom_sum$box == '0_8_3Mc', '3',
#                                                               ifelse(diatom_sum$box == '0_9_3', '3',
#                                                                      ifelse(diatom_sum$box == '12A_T11', "4",
#                                                                             ifelse(diatom_sum$box == '12A_T12', "4",
#                                                                                    ifelse(diatom_sum$box == '12A_T28', "10",
#                                                                                           ifelse(diatom_sum$box == '12A_T29', "10",
#                                                                                                  ifelse(diatom_sum$box == '12A_T30', "10",
#                                                                                                         ifelse(diatom_sum$box == '9', "6",
#                                                                                                                ifelse(diatom_sum$box == '9_17_6', "6",
#                                                                                                                       ifelse(diatom_sum$box == '9_18D_6', "6",
#                                                                                                                              NA)))))))))))))))
# 
# # add the section number
# diatom_sums$section <- NA
# diatom_sums$section <- ifelse(diatom_sum$box == "0_19_7", "a",
#                              ifelse(diatom_sum$box == "0_1_1Pa", "a",
#                                     ifelse(diatom_sum$box == "0_20_7", "a",
#                                            ifelse(diatom_sum$box == "0_21_7", "a",
#                                                   ifelse(diatom_sum$box == "0_3_1Da", "a",
#                                                          ifelse(diatom_sum$box == '0_8_3Mc', 'c',
#                                                                 ifelse(diatom_sum$box == '0_9_3', 'c',
#                                                                        ifelse(diatom_sum$box == '12A_T11', "a",
#                                                                               ifelse(diatom_sum$box == '12A_T12', "a",
#                                                                                      ifelse(diatom_sum$box == '12A_T28', "a",
#                                                                                             ifelse(diatom_sum$box == '12A_T29', "a",
#                                                                                                    ifelse(diatom_sum$box == '12A_T30', "a",
#                                                                                                           ifelse(diatom_sum$box == '9', "c",
#                                                                                                                  ifelse(diatom_sum$box == '9_17_6', "c",
#                                                                                                                         ifelse(diatom_sum$box == '9_18D_6', "c",
#                                                                                                                                NA)))))))))))))))

# move the last three columns to positions 3, 4, and 5
diatom_sums <- diatom_sums[, c(1:2, 27:29, 3:26, 29:ncol(diatom_sums))]
diatom_sums <- diatom_sums[, -c(30)]

##### PERMANOVA #####
## PERMANOVA ##
set.seed(2020)
## run permanova with marginal effects 
m.perm= adonis2(diatom_sums[,-c(1:7)] ~ blade_position, 
                data=diatom_sums, method = "bray")
m.perm
# Df SumOfSqs      R2      F Pr(>F)   
# blade_position  2  0.95968 0.34461 3.1548  0.004 **
#   Residual       12  1.82516 0.65539                 
# Total          14  2.78484 1.00000        

## POST HOC TEST ##
## post hoc test
m.pairwise = pairwise.adonis(diatom_sums[,-c(1:7)], diatom_sums$blade_position)
m.pairwise
# pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
# 1 B vs M  1 0.4901218 2.551507 0.2671314   0.080      0.240    
# 2 B vs T  1 0.6976320 3.748409 0.2940296   0.009      0.027   .
# 3 M vs T  1 0.2275071 2.901069 0.2661270   0.060      0.180    

## BETADISPERSON TEST ##
# calculate the distance
dist <- vegdist(diatom_sums[, 8:ncol(diatom_sums)], method = "bray")
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

bray_morph <- ggplot(diatom_sums, aes(x = NMDS1, y = NMDS2, color = blade_position)) +
  geom_point(size = 3) +
  scale_color_manual(values=c("#55D3A6", "#5582D3", "#A655D3"))+
  annotate("text", x=-1.7, y=1.25, label="morphological", size=3, color="black")+
  annotate("text", x = -1.7, y = 1.15, label = "stress = 0.02499072 (Bray-Curtis)", size=3, color="black")+
  annotate("text", x = -1.7, y = 1.05, label = "PERMANOVA: Section F= 3.15, p= 0.004", size=3, color="black")+
  annotate("text", x = -1.7, y = 0.95, label = "No sd", size=3, color="black")+
  annotate("text", x = -1.7, y = 0.85, label = "BETADISPERSON F= 6.12, p= 0.015", size=3, color="black") +
  labs(color = "Blade Section") 

ggarrange(plot_16S, plot_18S, plot_rbcl, bray_morph, labels = c("A", "B","C", "D"), ncol = 2, nrow = 2)

