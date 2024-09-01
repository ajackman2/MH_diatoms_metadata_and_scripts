# alpha diversity for the morphological counts
# A Jackman
# Jan 14 2024

# libraries
library(dplyr)
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
setwd("/Users/andreajackman/R_Stuff/Parfrey/diatoms_new/counting_data")
#path = "/Users/andreajackman/Desktop/R Stuff/diatoms_new/output/andrea/"

## read in data
diatom_data <- read.csv("corrected_master_diatom_data.csv")

# remove box_0_2_1 and box_12A_T10
diatom_data <- diatom_data |>
  dplyr::filter(box != '0_2_1') |>
  dplyr::filter(box != '12A_T10') 

# group by the box and sum columns and image_num
diatom_sum <- diatom_data |>
  dplyr::group_by(image_num, box, blade_position) |>
  dplyr::summarise(across(c(Navicula:Undatella), sum))

diatom_sum <- diatom_sum |>
  mutate(total = rowSums(across(where(is.numeric))))


# change the blade_position names
diatom_sum <- diatom_sum |>
  mutate(blade_position = case_when(
    blade_position == "proximal" ~ "B",
    blade_position == "distal" ~ "T",
    blade_position == "medial" ~ "M",
    TRUE ~ as.character(blade_position))) |>
  mutate(blade_position = as.factor(blade_position))

# change the blade_section of two sections
# change 12A_T12 from B to T
diatom_sum <- diatom_sum |>
  mutate(blade_position = if_else(box == "12A_T12", "T", blade_position))
# change 12A_T29
diatom_sum <- diatom_sum |>
  mutate(blade_position = if_else(box == "12A_T29", "M", blade_position))

# count the number of species
diatom_sum <- diatom_sum[, colSums(diatom_sum != 0) > 0]
diatom_sum <- diatom_sum |>
  mutate(species_num = rowSums(across(c(Navicula:Cocconeis), ~ . > 0)))

# average the number of species
diatom_sum <- diatom_sum %>%
  group_by(box, blade_position) %>%
  summarise(across(Navicula:species_num, mean, na.rm = TRUE))

# add a new column called leaf
diatom_sum <- diatom_sum |>
  mutate(leaf = str_extract(box, "^[^_]+"))

write.csv(diatom_sum, "diatom_sum_per_box_new.csv", row.names = FALSE)

# only the species
species_num <- diatom_sum |>
  group_by(box, blade_position) |>
  select(Navicula:Cocconeis) |>
  ungroup()

write.csv(species_num, "species_num_w_names_per_box.csv", row.names = FALSE)
species_num <- read.csv("species_num_w_names_per_box.csv")

species_num <- species_num |>
  select(Navicula:Cocconeis)

write.csv(species_num, "species_num_per_box.csv", row.names = FALSE)
diatom_sum <- read.csv("diatom_sum_per_box_new.csv")

set.seed(45)
##### FORMAT FOR AND GET ALPHA-STATS######
diatom_sum$shan = diversity(species_num, index="shannon")
diatom_sum$simp = diversity(species_num, index="simpson")

##### RUN ANOVA #####
arich = aov(diatom_sum$total~diatom_sum$blade_position)
summary(arich)
# Df  Sum Sq Mean Sq F value Pr(>F)  
# diatom_sum$blade_position  2 5335860 2667930   4.732 0.0305 *
# Residuals                 12 6765163  563764         
TukeyHSD(arich)
# diff        lwr      upr     p adj
# M-B 1120.1500  -223.5998 2463.900 0.1072110
# T-B 1339.0667   126.1028 2552.031 0.0305738
# T-M  218.9167 -1074.1072 1511.940 0.8945932

leveneTest(arich)
# Levene's Test for Homogeneity of Variance (center = median)
#       Df F value Pr(>F)
# group  2  0.1756 0.8411
#       12       
#plot(arich)

arich2 = aov(diatom_sum$species_num~diatom_sum$blade_position)

summary(arich2)
# Df Sum Sq Mean Sq F value Pr(>F)  
# diatom_sum$blade_position  2  45.49  22.745   6.167 0.0144 *
#   Residuals                 12  44.26   3.688      

TukeyHSD(arich2)
# diff        lwr      upr     p adj
# M-B 3.161111 -0.2758472 6.598069 0.0725931
# T-B 3.948148  0.8457058 7.050591 0.0136387
# T-M 0.787037 -2.5201777 4.094252 0.8041176

leveneTest(arich2)
# Levene's Test for Homogeneity of Variance (center = median)
#       Df F value Pr(>F)
# group  2  2.3399 0.1387
#       12                        
#plot(arich2)

ashan = aov(diatom_sum$shan~diatom_sum$blade_position)
summary(ashan)
# Df Sum Sq Mean Sq F value Pr(>F)
# diatom_sum$blade_position  2 0.1404 0.07018   0.323  0.729
# Residuals                 14 3.0399 0.21714     

leveneTest(ashan)
# Levene's Test for Homogeneity of Variance (center = median)
#       Df F value Pr(>F)
# group  2    0.31 0.7392
#       12                
#plot(ashan)

asimp = aov(diatom_sum$simp~diatom_sum$blade_position)
summary(asimp)
# Df Sum Sq Mean Sq F value Pr(>F)
# diatom_sum$blade_position  2 0.0674 0.03371   1.071  0.373
# Residuals                 12 0.3776 0.03146    
leveneTest(asimp)
# Levene's Test for Homogeneity of Variance (center = median)
#       Df F value Pr(>F)
# group  2   0.727 0.5035
#       12             

#plot(asimp)

###### PLOT ALPHA DIV #####
r=ggplot(diatom_sum, aes(x=blade_position, y=total))+
  #geom_point(size = 2, position=position_jitter(width=0.1))+
  geom_boxplot() +
  geom_point(alpha = 0.5) +
  annotate("text", x = c(1, 2, 3), y = c(160, 275, 250), label = c("a", "ab", 'b'), size = 5, color = 'black') +
  ylab("Number of Diatoms") +
  xlab("Blade Position") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14)) 

s= ggplot(diatom_sum, aes(x=blade_position, y=species_num))+
  #geom_point(size = 2, position=position_jitter(width=0.1))+
  geom_boxplot() +
  geom_point(alpha = 0.5) +
  annotate("text", x = c(1, 2, 3), y = c(6.4, 8.4, 9.6), label = c("a", "ab", 'b'), size = 5, color = 'black') +
  ylab("Number of Genera") +
  xlab("Blade Position") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14)) 

sh=ggplot(diatom_sum, aes(x=blade_position, y=shan))+
  #geom_point(size = 2, position=position_jitter(width=0.1))+
  geom_boxplot() +
  geom_point(alpha = 0.5) +
  annotate("text", x = c(1, 2, 3), y = c(1.3, 1.25, 1.75), label = c("a", "a", 'a'), size = 5, color = 'black') +
  ylab("Shannon's Index") +
  xlab("Blade Position") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

# si=ggplot(diatom_sum, aes(x=blade_position, y=simp, color = leaf))+
#   geom_point(size = 2, position=position_jitter(width=0.1))+
#   annotate("text", x = c(1, 2), y = c(10.75, 14), label = c("a", "a", 'a'), size = 4, color = 'darkred') +
#   ylab("Simpson's") +
#   xlab("Blade Position") 
#   # theme(axis.title = element_text(size = 16),  
#   #       axis.text = element_text(size = 14))

library(cowplot)

legend <- get_legend(r)

# Plot everything together
plot_grid(r, s, sh, legend, nrow = 1, ncol = 4, rel_widths = c(3, 3, 4, 1))

ggarrange(r, s, sh, ncol=3, nrow=1, labels = c("A", "B", "C"))

