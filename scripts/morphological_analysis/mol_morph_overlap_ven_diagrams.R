# investigating the overlap of molecular and morphological data

# load libraries
library(dplyr)
library(plyr)
library(tidyverse)
library(tidyr)
library(stringr)
library(phyloseq)
library(ggplot2); theme_set(theme_bw()+
                              theme(panel.grid = element_blank(),
                                    strip.background = element_rect(fill="white"),
                                    axis.text.y = element_text(size = 12, colour = "black"),
                                    axis.title = element_text(size=15, face="bold"),
                                    strip.text = element_text(color="black", size=10),
                                    legend.text=element_text(size=10),
                                    axis.line = element_line(colour = "black")))

## load functions
dephyloseq = function(phylo_obj){
  
  ## get the metadata
  meta = as.data.frame(as.matrix(phylo_obj@sam_data))
  
  ## how many metadata columns you have 
  metacols = ncol(meta)+1
  
  ## get out the otu table 
  ## if your metadta is empty after running this, you need to use 
  otu = as.data.frame(t(as.matrix(phylo_obj@otu_table)))
  #otu = as.data.frame(as.matrix(phylo_obj@otu_table))
  
  ## merge the metadata and otu table by the rownames (sample ids from the Illumina sequencing data)
  mo = merge(meta, otu, by=0)
  
  ## get out the taxonomy file 
  tax = as.data.frame(phylo_obj@tax_table)
  
  ## get the ASV ID out. This the matches the placeholder ASV ID in the OTU table
  tax = tax %>% rownames_to_column(var="ASVid")
  
  ## pivot longer to be able to match the ASVs in the OTU table to the taxonomy table 
  mo = mo %>% pivot_longer(cols = -c(1:metacols), names_to = "ASVid", values_to="asv_abundance")
  
  ## Join the metadata and otu table with the taoxnomy table 
  mot = full_join(mo, tax)
  
  ## Specify the output for the dephyloseq funciton 
  output = mot
}


# read in the data
diatom_data <- read.csv("counting_data/master_diatom_data.csv")
mol_18S = readRDS("imput_data/MH_18S_filtered001percent_notrarefied_phyloseq.rds")
mol_RBCL = readRDS("imput_data/MH_RBCL_filtered001percent_notrarefied_phyloseq.rds")

##### FORMAT SEM COUNTS ######
diatom_data = diatom_data |>  pivot_longer(cols=c(9:72),
                                           names_to = "Genus",
                                           values_to = "diatom_count")

diatom_data = ddply(diatom_data, c("Genus"),
                    summarise,
                    sum.sem.counts = sum(diatom_count),
                    n.sem.samples = length(unique(box)))

diatom_data$mean.sem.count = diatom_data$sum.sem.counts/diatom_data$n.sem.samples

diatom_data$in_counts = ifelse(diatom_data$mean.sem.count, "yes", "excluded")

# fixing one entry
diatom_data <- diatom_data %>%
  mutate(Genus = ifelse(Genus == "Paribellus", "Parlibellus", Genus))

##### FORMAT MOLECULAR DATA #####
## per sample read depth
mol_RBCL@sam_data$rd_filt = sample_sums(mol_RBCL)
mol_18S@sam_data$rd_filt = sample_sums(mol_18S)

## remove non-diatoms
mol_RBCL = subset_taxa(mol_RBCL, Phylum !="Eukaryota" &
                         Kingdom!="Sar") #Heterosigma akashiwo chloroplast DNA by BLAST

mol_18S = subset_taxa(mol_18S, phylum %in% c("Ochrophyta") &
                        !(class %in% c("Chrysophyceae", "Ochrophyta", "Phaeophyceae", "Pelagophyceae", "Raphidophyceae")))

## get out of phyloseq
mol_RBCL = dephyloseq(mol_RBCL)
mol_18S = dephyloseq(mol_18S) 

## get per sample relative abundance
mol_RBCL$ra = as.numeric(mol_RBCL$asv_abundance)/as.numeric(mol_RBCL$rd_filt)
mol_18S$ra = as.numeric(mol_18S$asv_abundance)/as.numeric(mol_18S$rd_filt)

## rename
names(mol_18S)[names(mol_18S)=="genus"]<-"Genus"
names(mol_18S)[names(mol_18S)=="order"]<-"Order"

inRBCL = c(unique(mol_RBCL$Genus))
in18S = c(unique(mol_18S$Genus))

mol = full_join(mol_RBCL, mol_18S)

mol = ddply(mol, c("Order", "Genus"),
            summarise,
            sum.mol.ra = sum(ra),
            mean.mol.ra = mean(ra),
            n.mol.samples =length(unique(Row.names)))

mol$in18S = ifelse(mol$Genus %in% c(in18S), "yes", "no")
mol$inRBCL = ifelse(mol$Genus %in% c(inRBCL), "yes", "no")

# remove not found in march from morphological data
diatom_data <- diatom_data %>%
  filter(!(Genus %in% c('Undatella', 'Actinoptychus', 'Trachyneis', 'Rhopalodia',
                        'Pseudonitzschia', 'Parlibellus', 'Hobaniella', 'Hanzschia',
                        'Gyrosigma', 'Cyclotella')))

###### combine data #####
all = full_join(mol, diatom_data)

# fill in order for some the genera
# Replace 'Achnanthales' for 'Achnanthes' in the Genus column
all$Order[all$Genus == "Achnanthes"] <- "Achnanthales"
all$Order[all$Genus == "Dimeregramma"] <- "Plagiogrammales"
all$Order[all$Genus == "Ditylum"] <- "Lithodesmiales"
all$Order[all$Genus == "Donkinia"] <- "Naviculales"
all$Order[all$Genus == "Encyonema"] <- "Cymbellales"
all$Order[all$Genus == "Entomoneis"] <- "Surirellales"
all$Order[all$Genus == "Epithemia"] <- "Rhopalodiales"
all$Order[all$Genus == "Eucampia"] <- "Hemiaulales"
all$Order[all$Genus == "Eupyxidicula"] <- "Bacillariophytina"
all$Order[all$Genus == "Fogedia"] <- "Naviculales"
all$Order[all$Genus == "Fragilariopsis"] <- "Bacillariales"
all$Order[all$Genus == "Gomphoseptatum"] <- "Cymbellales"
all$Order[all$Genus == "Grammatophora"] <- "Rhabdonematales"
all$Order[all$Genus == "Opephora"] <- "Fragilariales"
all$Order[all$Genus == "Paralia"] <- "Paraliales"
all$Order[all$Genus == "Pariraphis"] <- "Thalassiophysales"
all$Order[all$Genus == "Petroneis"] <- "Lyrellales"
all$Order[all$Genus == "Plagiogramma"] <- "Plagiogrammales"
all$Order[all$Genus == "Plagiotropis"] <- "Naviculales"
all$Order[all$Genus == "Podosira"] <- "Melosirales"
all$Order[all$Genus == "Psammodictyon"] <- "Bacillariales"
all$Order[all$Genus == "Rhoicosphenia"] <- "Cymbellales"
all$Order[all$Genus == "Trigonium"] <- "Biddulphiaceae"

#write.csv(all, "counting_data/combined_sem_illumina_counts_all_data.csv")

# read the data in
overlap_data <- all

# edit the data to add in no and yes
ven_data <- overlap_data
ven_data$in18S <- ifelse(is.na(ven_data$in18S), "no", ven_data$in18S)
ven_data$inRBCL <- ifelse(is.na(ven_data$inRBCL), "no", ven_data$inRBCL)
ven_data$in_counts <- ifelse(is.na(ven_data$in_counts), "no", ven_data$in_counts)
ven_data$in_counts <- gsub("excluded", "no", ven_data$in_counts)

# making a vendiagram
library(VennDiagram)

# make a new col called in_molecular
ven_data$in_molecular <- ifelse(ven_data$in18S == "yes" | ven_data$inRBCL == "yes", "yes", "no")

# make a new vec with the in_molecular and in_counts
molecular_genera_vec <- c()
morpho_genera_vec <- c()

for (i in 1:nrow(ven_data)) {
  if (ven_data$in_molecular[i] == "yes") {
    molecular_genera_vec <- c(molecular_genera_vec, ven_data$Genus[i])
  }
}

for (i in 1:nrow(ven_data)) {
  if (ven_data$in_counts[i] == "yes") {
    morpho_genera_vec <- c( morpho_genera_vec, ven_data$Genus[i])
  }
}

venn.plot <- venn.diagram(
  x = list(Molecular = molecular_genera_vec, Morphological = morpho_genera_vec),
  category.names = c("Molecular", "Morphological"),
  filename = NULL
)


grid.draw(venn.plot)

# include Mark's data

ven_data2 <- overlap_data
ven_data2$in18S <- ifelse(is.na(ven_data2$in18S), "no", ven_data2$in18S)
ven_data2$inRBCL <- ifelse(is.na(ven_data2$inRBCL), "no", ven_data2$inRBCL)
ven_data2$in_counts <- ifelse(is.na(ven_data2$in_counts), "no", ven_data2$in_counts)
ven_data2$in_marks <- ifelse(ven_data2$in_counts %in% c("excluded", "yes"), "yes", "no")
ven_data2$in_counts <- gsub("excluded", "no", ven_data2$in_counts)

# make a new col called in_molecular
ven_data2$in_molecular <- ifelse(ven_data2$in18S == "yes" | ven_data2$inRBCL == "yes", "yes", "no")

# make a new vec with the in_molecular and in_counts
molecular_genera_vec2 <- c()
morpho_genera_vec2 <- c()
marks_genera_vec <- c()

for (i in 1:nrow(ven_data2)) {
  if (ven_data2$in_molecular[i] == "yes") {
    molecular_genera_vec2 <- c(molecular_genera_vec2, ven_data2$Genus[i])
  }
}

for (i in 1:nrow(ven_data2)) {
  if (ven_data2$in_counts[i] == "yes") {
    morpho_genera_vec2 <- c( morpho_genera_vec2, ven_data2$Genus[i])
  }
}

for (i in 1:nrow(ven_data2)) {
  if (ven_data2$in_marks[i] == "yes") {
    marks_genera_vec <- c(marks_genera_vec, ven_data2$Genus[i])
  }
}

venn.plot2 <- venn.diagram(
  x = list(Molecular = molecular_genera_vec2, Morphological = morpho_genera_vec2, Mark = marks_genera_vec),
  category.names = c("Molecular", "Morphological Count", "Mark's Count"),
  filename = NULL
)


grid.draw(venn.plot2)

# ven3
# include Mark's data

ven_data3 <- overlap_data
ven_data3$in18S <- ifelse(is.na(ven_data3$in18S), "no", ven_data3$in18S)
ven_data3$inRBCL <- ifelse(is.na(ven_data3$inRBCL), "no", ven_data3$inRBCL)
ven_data3$in_counts <- ifelse(is.na(ven_data3$in_counts), "no", ven_data3$in_counts)
ven_data3$in_counts <- ifelse(ven_data3$in_counts == "excluded" | ven_data3$in_counts == "yes", "yes", "no")


# make a new col called in_molecular
ven_data3$in_molecular <- ifelse(ven_data3$in18S == "yes" | ven_data3$inRBCL == "yes", "yes", "no")

# make a new vec with the in_molecular and in_counts
molecular_genera_vec3 <- c()
morpho_genera_vec3 <- c()

for (i in 1:nrow(ven_data3)) {
  if (ven_data3$in_molecular[i] == "yes") {
    molecular_genera_vec3 <- c(molecular_genera_vec3, ven_data3$Genus[i])
  }
}

for (i in 1:nrow(ven_data3)) {
  if (ven_data3$in_counts[i] == "yes") {
    morpho_genera_vec3 <- c( morpho_genera_vec3, ven_data3$Genus[i])
  }
}


venn.plot3 <- venn.diagram(
  x = list(Molecular = molecular_genera_vec3, Morphological = morpho_genera_vec3),
  #category.names = c("Molecular", "Morphological"),
  filename = NULL,
  category.cex = 2, 
  category.pos = c(20, 20)
)


grid.newpage()
grid.draw(venn.plot3)

grob_tree <- grid.grabExpr(expr = venn.plot3)

molecular_grob <- grob_tree$children[[which(grob_tree$name == "Molecular")]]
morphological_grob <- grob_tree$children[[which(grob_tree$name == "Morphological")]]

molecular_grob$x <- unit(0.2, "npc")
molecular_grob$y <- unit(0.5, "npc")

morphological_grob$x <- unit(0.8, "npc")
morphological_grob$y <- unit(0.5, "npc")


grid.draw(grob_tree)

# make a dataframe where mark's counts and inside the morphological counts
ven_data4 <- overlap_data
ven_data4$in18S <- ifelse(is.na(ven_data4$in18S), "no", ven_data4$in18S)
ven_data4$inRBCL <- ifelse(is.na(ven_data4$inRBCL), "no", ven_data4$inRBCL)
ven_data4$in_counts <- ifelse(is.na(ven_data4$in_counts), "no", ven_data4$in_counts)
ven_data4$in_marks <- ifelse(ven_data4$in_counts == "excluded", "yes", "no")
ven_data4$in_counts <- ifelse(ven_data4$in_counts == "excluded" | ven_data4$in_counts == "yes", "yes", "no")


# make a new col called in_molecular
ven_data4$in_molecular <- ifelse(ven_data4$in18S == "yes" | ven_data4$inRBCL == "yes", "yes", "no")

# make a new vec with the in_molecular and in_counts
molecular_genera_vec4 <- c()
morpho_genera_vec4 <- c()
marks_genera_vec <- c()

for (i in 1:nrow(ven_data4)) {
  if (ven_data4$in_molecular[i] == "yes") {
    molecular_genera_vec4 <- c(molecular_genera_vec4, ven_data4$Genus[i])
  }
}

for (i in 1:nrow(ven_data4)) {
  if (ven_data4$in_counts[i] == "yes") {
    morpho_genera_vec4 <- c( morpho_genera_vec4, ven_data4$Genus[i])
  }
}

for (i in 1:nrow(ven_data4)) {
  if (ven_data4$in_marks[i] == "yes") {
    marks_genera_vec <- c(marks_genera_vec, ven_data4$Genus[i])
  }
}


venn.plot4 <- venn.diagram(
  x = list(Molecular = molecular_genera_vec4, Morphological_Counts = morpho_genera_vec4, Morphological = marks_genera_vec),
  #category.names = c("Molecular", "Morphological"),
  filename = NULL,
  category.cex = 2, 
  category.pos = c(20, 20)
)


grid.newpage()
grid.draw(venn.plot4)

# Filter out rows with 'NA' or zero ranks
# overlap_data <- overlap_data[!(is.na(overlap_data$sum.mol.ra) | overlap_data$sum.mol.ra == 0),]
# overlap_data <- overlap_data[!(is.na(overlap_data$sum.sem.counts) | overlap_data$sum.sem.counts == 0),]

# Rank genera based on molecular and morphological data separately
overlap_data <- overlap_data[order(overlap_data$sum.mol.ra, decreasing = TRUE),]
overlap_data$molecular_rank <- rank(-overlap_data$sum.mol.ra, ties.method = "min")

overlap_data <- overlap_data[order(overlap_data$sum.sem.count, decreasing = TRUE),]
overlap_data$morphological_rank <- rank(-overlap_data$sum.sem.count, ties.method = "min")

# Define colors based on presence in molecular or morphological data


# Create a scatter plot to compare ranks
plot <- ggplot(overlap_data, aes(x = molecular_rank, y = morphological_rank)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(title = "Comparison of Molecular and Morphological Abundance Ranks",
       x = "Molecular Rank", y = "Morphological Rank") 

# Print the plot
print(plot)

correlation <- cor(overlap_data$molecular_rank, overlap_data$morphological_rank)
correlation # 0.42

correlation_test <- cor.test(overlap_data$molecular_rank, overlap_data$morphological_rank)
correlation_test
#data:  overlap_data$molecular_rank and overlap_data$morphological_rank
# t = 1.9254, df = 17, p-value = 0.07107
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.0384875  0.7358994
# sample estimates:
#   cor 
# 0.4231185 

color_scale <- scale_color_gradient(low = "blue", high = "red")

ven_data3$where_found <- ifelse(ven_data3$in_counts == "yes" & ven_data3$in_molecular == "yes", "in_both",
                       ifelse(ven_data3$in_counts == "yes" & ven_data3$in_molecular == "no", "in_morpho",
                              ifelse(ven_data3$in_counts == "no" & ven_data3$in_molecular == "yes", "in_molecular", NA)))

ven_data3 <- ven_data3[order(ven_data3$sum.mol.ra, decreasing = TRUE),]
ven_data3$molecular_rank <- rank(-ven_data3$sum.mol.ra, ties.method = "min")

ven_data3 <- ven_data3[order(ven_data3$sum.sem.count, decreasing = TRUE),]
ven_data3$morphological_rank <- rank(-ven_data3$sum.sem.count, ties.method = "min")

# Create the scatter plot
scatter_plot <- ggplot(ven_data3, aes(x = molecular_rank, y = morphological_rank, label = Genus, color = where_found)) +
  geom_point() +
  labs(title = "Comparison of Molecular and Morphological Abundance Ranks",
       x = "Molecular Rank", y = "Morphological Rank") +
  theme_minimal()

scatter_plot
venn.plot <- venn.diagram(
  x = list(
    Molecular = molecular_genera_vec,
    Morphological = morpho_genera_vec
  ),
  category.names = c("Molecular", "Morphological"),
  filename = NULL,
)

# your data


# have a look at the default plot
grid.newpage()
grid.draw(venn.plot)

# have a look at the names in the plot object v
lapply(venn.plot,  names)
# We are interested in the labels
lapply(venn.plot, function(i) i$label)

# overwrite labels
venn.plot[[5]]$label  <- paste(setdiff(molecular_genera_vec, morpho_genera_vec), collapse="\n")  

venn.plot[[6]]$label <- paste(setdiff(morpho_genera_vec, molecular_genera_vec)  , collapse="\n")  

venn.plot[[7]]$label <- paste(intersect(molecular_genera_vec, morpho_genera_vec), collapse="\n")  


# plot  
grid.newpage()
grid.draw(venn.plot)


# Overwrite labels with two columns
venn.plot[[5]]$label <- paste(setdiff(molecular_genera_vec, morpho_genera_vec), collapse="\n")
venn.plot[[6]]$label <- paste(setdiff(morpho_genera_vec, molecular_genera_vec), collapse="\n")
venn.plot[[7]]$label <- paste(intersect(molecular_genera_vec, morpho_genera_vec), collapse="\n")

# Plot the Venn diagram
grid.newpage()
grid.draw(venn.plot)


library(ggvenn)


x <- list(Molecular=molecular_genera_vec , Morphological_count = morpho_genera_vec)

ggvenn(x, show_elements = T, label_sep = "\n", count_column = 4, auto_scale = TRUE, text_size = 2)

ggvenn(x, show_elements = TRUE, label_sep = "\n", text_size = 1, fill_alpha = 0)

ggvenn(x, show_elements = TRUE, label_sep = "\n", text_size = 2, fill_alpha = 0) +



