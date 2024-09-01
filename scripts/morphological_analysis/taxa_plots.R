# taxa plots
# A. Jackman
# Nov 28 2023

# libraries
library(dplyr)
library(tidyverse)
library(plyr)
library(tidyr)
library(ggplot2); theme_set(theme_bw()+
                              theme(panel.grid = element_blank(),
                                    strip.background = element_rect(fill="white"),
                                    axis.text.y = element_text(size = 12, colour = "black"),
                                    axis.title = element_text(size=15, face="bold"),
                                    strip.text = element_text(color="black", size=10),
                                    legend.text=element_text(size=10),
                                    axis.line = element_line(colour = "black"),
                                    axis.text.x = element_blank(),))


## set working directory
setwd("/Users/andreajackman/R_Stuff/Parfrey/diatoms_new")


# read in the data
diatom_data <- read.csv("counting_data/master_diatom_data.csv")

diatom_sum = diatom_data |>
  select(file_name, box, blade_position, Navicula:Undatella)

genus_columns <- diatom_sum[, c("Navicula", "Tabularia", "Cylindrotheca", "Pseudonitzschia", "Gomphonemopsis", "Hyalodiscus", "Coscinodiscus",   
                                "Plagiotropis", "Nitzschia", "Rhoicosphenia", "Pseudogomphonema", "Thalassiosira",   
                                "Eucampia", "Chaetoceros", "Diploneis", "Entomoneis", "Pariraphis",    
                                "Planothidium", "Halamphora", "Achnanthes", "Amphora", "Petroneis",     
                                "Bacillaria", "Ellerbeckia", "Minidiscus", "Fragilariopsis", "Rhizosolenia",   
                                "Fogedia", "Skeletonema", "Pleurosigma", "Gomphoseptatum", "Licmophora",     
                                "Cocconeis", "Actinoptychus", "Attheya", "Cyclotella", "Dimeregramma",   
                                "Ditylum", "Donkinia", "Epithemia", "Eupyxidicula", "Encyonema",       
                                "Fallacia", "Grammatophora", "Gyrosigma", "Hanzschia", "Haslea",          
                                "Hobaniella", "Leptocylindrus", "Lyrella", "Melosira", "Odontella",      
                                "Opephora", "Paralia", "Paribellus", "Plagiogramma", "Podosira",        
                                "Psammodictyon", "Rhabdonema", "Rhopalodia", "Trigonium", "Trachyneis",      
                                "Tryblionella", "Undatella")]

sum_per_genus <- colSums(genus_columns, na.rm = TRUE)

# make a dataframe
total_counts_table <- data.frame(
  genus = names(sum_per_genus),
  sum_value = sum_per_genus
)

#View(total_counts_table)

# now separate it by the blade position
counts_position <- diatom_data |>
  select(file_name, box, blade_position, Navicula:Undatella) |>
  group_by(blade_position) |>
  summarise_at(vars(Navicula:Undatella), ~sum(., na.rm = TRUE))
#View(counts_position)

# make a new table
counts_pos_table <- counts_position |>
  pivot_longer(cols = Navicula:Undatella, names_to = "genus", values_to = "sum")

counts_pos_table <- counts_pos_table |>
  pivot_wider(names_from = blade_position, values_from = sum)
#View(counts_pos_table)

# pivot longer
counts_position_long <- diatom_data |>
  select(file_name, box, blade_position, Navicula:Undatella) |>
  group_by(box, blade_position) |>
  summarise_at(vars(Navicula:Undatella), ~sum(., na.rm = TRUE)) |>
  pivot_longer(cols = Navicula:Undatella, names_to = "genus", values_to = "sum")
#View(counts_position_long)

##### LOOP TO SELECT TOP 15 GENERA IN SAMPLES ######
## make variable to track the section - proximal, medial, distal

# replace the NAs with zeros
counts_position_long <- counts_position_long |>
  replace_na(replace = list(relative_proportion = 0))

# manually put the box totals in
# box_9 = 3890
# box_0_19_7 = 313
# box_9_17_6 = 1640
# box_12A_T10 = 0
# box_12A_T11 = 1639
# box_12A_T12 = 1898
# box_12A_T28 = 25
# box_12A_T30 = 4379
# box_0_1_1Pa = 1478
# box_0_3_1Da = 1172
# box_0_8_3Mc = 279
# box_9_18D_6 = 2221
# box_0_21_7 = 2076
# box_0_20_7 = 2725
# box_0_9_3 = 1734
# box_0_2_1 = 1767
# box_12A_T29 = 1821

counts_position_long <- counts_position_long %>%
  mutate(box = paste0('box_', box))

box_totals <- data.frame(
  box = c("box_9", "box_0_19_7", "box_9_17_6", "box_12A_T10", "box_12A_T11", "box_12A_T12",
          "box_12A_T28", "box_12A_T30", "box_0_1_1Pa", "box_0_3_1Da", "box_0_8_3Mc",
          "box_9_18D_6", "box_0_21_7", "box_0_20_7", "box_0_9_3", "box_0_2_1", "box_12A_T29"),
  total = c(3890, 312, 1652, 0, 1639, 1911, 25, 4680, 1478, 1205, 279, 2465, 2080, 2752, 1746, 1767, 1823)
)

data <- counts_position_long %>%
  left_join(box_totals, by = "box")

## get relative abundance of diatoms in the box
data <- data |>
  mutate(relative_proportion = sum / total)

data$blade_position <- ifelse(data$box == "box_12A_T12", "distal",
                              ifelse(data$box == "box_12A_T29", "medial", data$blade_position))

# rough plot
ggplot(data, aes(x = box, y = relative_proportion, fill = genus)) +
       geom_bar(stat = "identity", position = "stack") +
       labs( x = "box", 
             y = "Sum") +
       facet_wrap(~blade_position, scales = "free_x") +
       theme_minimal() +
       theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
       coord_cartesian(ylim = c(0, 1))

# remove box_0_2_1 and box_12A_T10
data <- data |>
  filter(box != 'box_0_2_1') |>
  filter(box != 'box_12A_T10')


## summarize mto
# data.sum = ddply(data, c("blade_position","genus"),
#                 summarise,
#                 N=length(relative_abundance),
#                 sum = sum(relative_abundance))


samplegroups = unique(data$blade_position)

## sort data by relative abundance
sorted = data[order(-data$relative_proportion),]

## make empty dataframe to store taxa
top.df = NULL

genus_values <- c("Cocconeis", "Tabularia", "Navicula", "Gomphonemopsis", "Pariraphis", "Rhoicosphenia", 
                  "Nitzschia", "Cylindrotheca", "Pseudogomphonema", "Planothidium", "Plagiotropis", 
                  "Halamphora", "Achnanthes", "Thalassiosira", "Skeletonema", "Minidiscus", "Amphora", 
                  "Coscinodiscus", "Pseudonitzschia")


# Add a new column 'place' if genus is in the top 15
data <- data |>
  mutate(place = ifelse(genus %in% genus_values, "top_19", NA)) |>
  mutate(taxaplotgroup = blade_position)


##### SET UP TAXA PLOT #####

## join the top taxa and existing mto dataframe
diatoms.top = data

## make the empty "place" cells say bottom
diatoms.top$place = replace(diatoms.top$place, is.na(diatoms.top$place), "bottom")

## replace plot_names for bottom taxa with Other
diatoms.top[diatoms.top$place == "bottom",]$genus <- "Others"



##### MAKE COLOR LIST #####
## how many do we need?
taxa = unique(diatoms.top$genus)

colors = read.csv("scripts/scripts_Andrea/taxaplot_colors.csv")
colors

##### MAKE TAXAPLOT #####
# Assign taxa names to colors
scolors <- colors$plot_colors
names(scolors) <- colors$genus

## order mot.top by decreasing relative abundance
diatoms.top = diatoms.top[order(-diatoms.top$relative_proportion),]

## get list of factors in order
natural.genus.order = as.list(c(unique(diatoms.top$genus)))

## remove others from list #!#
no.others=natural.genus.order[!natural.genus.order == 'Others']

## add Others to end of list
plot.order = append(no.others, "Others")

## set plot_names levels
plot.order = unlist(plot.order)

## order dataframe by relative abundance
diatoms.top$genus = factor(diatoms.top$genus, levels=c(plot.order))

## get list to cycle through
diatoms.top$taxaplotgroup = diatoms.top$blade_position
taxaplot.groups = unique(diatoms.top$taxaplotgroup)

diatoms.top <- diatoms.top %>%
  mutate(blade_position = dplyr::recode(blade_position, "proximal" = "B", "distal" = "T", "medial" = "M"))

#filepath = "visualizations"
## field taxaplot loop
for (i in taxaplot.groups){
  for (j in i){
    #sub.df = subset(diatoms.top, diatoms.top$taxaplotgroup== c(j))
    
    myplot=ggplot(diatoms.top, aes(x = as.character(box), y=as.numeric(relative_proportion), 
                              fill=as.factor(genus)))+
      geom_bar(stat = "identity")+
      scale_fill_manual(values=scolors)+
      facet_grid(.~blade_position,
      scales="free", space="free") +
      labs(y="Genus relative abundance in samples", x = "Sample", fill="Genus")+
      guides(fill=guide_legend(ncol=1)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 0))

    
    myplot
    
  }
}
myplot


#### molecular taxa plot ####
###### SET UP #######
## load packages
library(phyloseq)
library(tidyverse)
library(plyr)
library(qualpalr)
library(ggpubr)
library(ggh4x)
library(ggplot2); theme_set(theme_bw()+
                              theme(panel.grid = element_blank(),
                                    strip.background = element_rect(fill="white"),
                                    axis.text.y = element_text(size = 10, colour = "black"),
                                    axis.title = element_text(size=10, face="bold"),
                                    strip.text = element_text(color="black", size=10),
                                    legend.text=element_text(size=10),
                                    axis.line = element_line(colour = "black"),
                                    axis.text.x = element_blank()))

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


# read in RBCL and 18S data
mol_RBCL = readRDS("imput_data/MH_RBCL_filtered001percent_notrarefied_phyloseq.rds")
mol_18S = readRDS("imput_data/MH_18S_filtered001percent_notrarefied_phyloseq.rds")

## remove non-diatoms
mol_RBCL = subset_taxa(mol_RBCL, Phylum !="Eukaryota" &
                         Kingdom!="Sar") #Heterosigma akashiwo chloroplast DNA by BLAST

mol_18S = subset_taxa(mol_18S, phylum %in% c("Ochrophyta") &
                        !(class %in% c("Chrysophyceae", "Ochrophyta", "Phaeophyceae", "Pelagophyceae", "Raphidophyceae")))

###### PREP FOR TAXAPLOT #####
## summarize at rank 6
mol_RBCL@sam_data$rd_filt = sample_sums(mol_RBCL)
mol_18S@sam_data$rd_filt = sample_sums(mol_18S)

mol_RBCL = tax_glom(mol_RBCL, taxrank = "Genus")
mol_18S = tax_glom(mol_18S, taxrank = "genus")

## get data out of phyloseq and inta a dataframe
mol_RBCL = dephyloseq(mol_RBCL)
mol_18S = dephyloseq(mol_18S)

# filter mol_18S for the diatoms
mol_18S <- mol_18S |>
  filter(class == "Diatomea")

# filter out the higher order classifications
mol_18S <- mol_18S |>
  filter(!(genus %in% c('Achnanthales', 'Bacillariales', 'Bacillariophyceae', 'Bacillariophytina',
                        'Bacillariophyta', 'Coscinodiscophyceae', 'Cymbellales', 'Diatomea',
                        'Fragilariales', 'Fragilariophyceae', 'Mediophyceae', 'Melosirids',
                        'Naviculales', 'Thalassiophysales', 'Thalassiosirales')))
  
# fix the columns
mol_18S <- mol_18S |>
  select(-domain, -phylum, -class, -order, -family, -species) |>
  mutate(plant_id = as.character(plant_id))

mol_18S$database <- '18S'

# filter out the higher order classifications
mol_RBCL<- mol_RBCL |>
  filter(!(Genus %in% c('Achnanthales', 'Bacillariales', 'Bacillariophyceae', 'Bacillariophytina',
                        'Bacillariophyta', 'Coscinodiscophyceae', 'Cymbellales', 'Diatomea',
                        'Fragilariales', 'Fragilariophyceae', 'Mediophyceae', 'Melosirids',
                        'Naviculales', 'Thalassiophysales', 'Thalassiosirales')))

mol_RBCL <- mol_RBCL |>
  mutate(plant_id = as.character(plant_id))
mol_RBCL$database <- 'RBCL'

names(mol_18S)[names(mol_18S)=="genus"]<-"Genus"
inRBCL = c(unique(mol_RBCL$Genus))
in18S = c(unique(mol_18S$Genus))

seagrassdf = mol_18S #full_join(mol_RBCL, mol_18S)

seagrassdf <- seagrassdf |>
  select(-asv_id)
  #select(-Species, -subsp, -asv_id, -Domain, -Kingdom, -Subkingdom, -Phylum, -Class, -Order, -Empty)

seagrassdf <- na.omit(seagrassdf)

## get per sample relative abundance
sum_per_sample <- aggregate(asv_abundance ~ Row.names, data = seagrassdf, FUN = sum)

# Merge the sum back to the original dataframe
seagrassdf <- merge(seagrassdf, sum_per_sample, by.x = "Row.names", by.y = "Row.names", suffixes = c("", "_sum"))

# Rename the new column
colnames(seagrassdf)[which(colnames(seagrassdf) == "asv_abundance_sum")] <- "rd_filt_new"

seagrassdf$ra = as.numeric(seagrassdf$asv_abundance) / as.numeric(seagrassdf$rd_filt_new)

## make plotnames
seagrassdf$plotnames = paste0(seagrassdf$Genus)

## make group for taxaplot
# make list of groups
grouplist = c(unique(seagrassdf$leaf_section))

## summarize data by taxaplot geoup type. 
seagrass.sum = ddply(seagrassdf, c("leaf_section", "plotnames"),
                     summarise,
                     sum = sum(ra))

## sort data by relative abundance. This is how the loop will pick the mos tabundant taxa
sorted = seagrass.sum[order(-seagrass.sum$sum),]

##### CALULCUATE MORE ABUNDANT TAXA ######
## make empty dataframe to store output from the loop
top.df = NULL

## start loop
for(i in grouplist) {
  for(j in i) {
    
    ## subset dataframe by samples
    #!# Remeber to change te substrate to your group! 
    sample = subset(sorted, sorted$leaf_section %in% c(j))
    
    ## get top 15 genera
    top = sample[c(1:15),]
    
    ## save list of top  abundance taxa
    t.tmp <- top
    top.df <- rbind.fill(top.df, t.tmp)
    
    ## close loop 
  }
}

top.df$place = "top_19"
alldata = full_join(seagrassdf, top.df)

alldata$place = replace(alldata$place, is.na(alldata$place), "bottom")

## replace plot_names for bottom taxa with Other
alldata[alldata$place == "bottom",]$plotnames <- "Others"

##### MAKE COLOR LIST #####
## how many do we need?
taxa = unique(alldata$plotnames)

## write a csv with the names
colors = read.csv("scripts/scripts_Andrea/taxaplot_colors_18S.csv")
colors

##### MAKE TAXAPLOT #####
# Assign taxa names to colors
rcolors <- colors$plot_colors
names(rcolors) <- colors$plot_names

##### ORDER THE TAXA SO OTHERS ARE AT THE BOTTOM #####
## order by decreasing relative abundance
alldata = alldata[order(-alldata$ra),]

## get list of factors in order
natural.genus.order = as.list(c(unique(alldata$plotnames)))

## remove others from list #!#
no.others=natural.genus.order[!natural.genus.order == 'Others']

## add Others to end of list
plot.order = append(no.others, "Others")

## set plot_names levels
plot.order = unlist(plot.order)

## order dataframe by relative abundance
alldata$plotnames = factor(alldata$plotnames, levels=c(plot.order))


##### MAKE THE TAXAPLOTS ######

## make the plot
mol_18S_plot <- ggplot(alldata, aes(x=as.character(Row.names), y=as.numeric(ra), 
                                fill=as.factor(plotnames)))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values=rcolors)+
  guides(fill=guide_legend(ncol=1))+
  facet_nested(.~leaf_section + database, scales="free")+
  labs(y="Relative Abundance", x="Sample", fill="Genus ") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 0))

# read in RBCL and 18S data
mol_RBCL = readRDS("imput_data/MH_RBCL_filtered001percent_notrarefied_phyloseq.rds")
mol_18S = readRDS("imput_data/MH_18S_filtered001percent_notrarefied_phyloseq.rds")

## remove non-diatoms
mol_RBCL = subset_taxa(mol_RBCL, Phylum !="Eukaryota" &
                         Kingdom!="Sar") #Heterosigma akashiwo chloroplast DNA by BLAST

mol_18S = subset_taxa(mol_18S, phylum %in% c("Ochrophyta") &
                        !(class %in% c("Chrysophyceae", "Ochrophyta", "Phaeophyceae", "Pelagophyceae", "Raphidophyceae")))

###### PREP FOR TAXAPLOT #####
## summarize at rank 6
mol_RBCL@sam_data$rd_filt = sample_sums(mol_RBCL)
mol_18S@sam_data$rd_filt = sample_sums(mol_18S)

mol_RBCL = tax_glom(mol_RBCL, taxrank = "Genus")
mol_18S = tax_glom(mol_18S, taxrank = "genus")

## get data out of phyloseq and inta a dataframe
mol_RBCL = dephyloseq(mol_RBCL)
mol_18S = dephyloseq(mol_18S)

# filter mol_18S for the diatoms
mol_18S <- mol_18S |>
  filter(class == "Diatomea")

# filter out the higher order classifications
mol_18S <- mol_18S |>
  filter(!(genus %in% c('Achnanthales', 'Bacillariales', 'Bacillariophyceae', 'Bacillariophytina',
                        'Bacillariophyta', 'Coscinodiscophyceae', 'Cymbellales', 'Diatomea',
                        'Fragilariales', 'Fragilariophyceae', 'Mediophyceae', 'Melosirids',
                        'Naviculales', 'Thalassiophysales', 'Thalassiosirales')))

# fix the columns
mol_18S <- mol_18S |>
  select(-domain, -phylum, -class, -order, -family, -species) |>
  mutate(plant_id = as.character(plant_id))

mol_18S$database <- '18S'

# filter out the higher order classifications
mol_RBCL<- mol_RBCL |>
  filter(!(Genus %in% c('Achnanthales', 'Bacillariales', 'Bacillariophyceae', 'Bacillariophytina',
                        'Bacillariophyta', 'Coscinodiscophyceae', 'Cymbellales', 'Diatomea',
                        'Fragilariales', 'Fragilariophyceae', 'Mediophyceae', 'Melosirids',
                        'Naviculales', 'Thalassiophysales', 'Thalassiosirales')))

mol_RBCL <- mol_RBCL |>
  mutate(plant_id = as.character(plant_id))
mol_RBCL$database <- 'RBCL'

names(mol_18S)[names(mol_18S)=="genus"]<-"Genus"
inRBCL = c(unique(mol_RBCL$Genus))
in18S = c(unique(mol_18S$Genus))

seagrassdf = mol_RBCL #full_join(mol_RBCL, mol_18S)

seagrassdf <- seagrassdf |>
  #select(-asv_id)
  select(-Species, -subsp, -asv_id, -Domain, -Kingdom, -Subkingdom, -Phylum, -Class, -Order, -Empty)

seagrassdf <- na.omit(seagrassdf)

## get per sample relative abundance
sum_per_sample <- aggregate(asv_abundance ~ Row.names, data = seagrassdf, FUN = sum)

# Merge the sum back to the original dataframe
seagrassdf <- merge(seagrassdf, sum_per_sample, by.x = "Row.names", by.y = "Row.names", suffixes = c("", "_sum"))

# Rename the new column
colnames(seagrassdf)[which(colnames(seagrassdf) == "asv_abundance_sum")] <- "rd_filt_new"

seagrassdf$ra = as.numeric(seagrassdf$asv_abundance) / as.numeric(seagrassdf$rd_filt_new)

## make plotnames
seagrassdf$plotnames = paste0(seagrassdf$Genus)

## make group for taxaplot
# make list of groups
grouplist = c(unique(seagrassdf$leaf_section))

## summarize data by taxaplot geoup type. 
seagrass.sum = ddply(seagrassdf, c("leaf_section", "plotnames"),
                     summarise,
                     sum = sum(ra))

## sort data by relative abundance. This is how the loop will pick the mos tabundant taxa
sorted = seagrass.sum[order(-seagrass.sum$sum),]

##### CALULCUATE MORE ABUNDANT TAXA ######
## make empty dataframe to store output from the loop
top.df = NULL

## start loop
for(i in grouplist) {
  for(j in i) {
    
    ## subset dataframe by samples
    #!# Remeber to change te substrate to your group! 
    sample = subset(sorted, sorted$leaf_section %in% c(j))
    
    ## get top 15 genera
    top = sample[c(1:17),]
    
    ## save list of top  abundance taxa
    t.tmp <- top
    top.df <- rbind.fill(top.df, t.tmp)
    
    ## close loop 
  }
}

top.df$place = "top_19"
alldata = full_join(seagrassdf, top.df)

alldata$place = replace(alldata$place, is.na(alldata$place), "bottom")

## replace plot_names for bottom taxa with Other
alldata[alldata$place == "bottom",]$plotnames <- "Others"

##### MAKE COLOR LIST #####
## how many do we need?
taxa = unique(alldata$plotnames)

## write a csv with the names
colors = read.csv("scripts/scripts_Andrea/taxaplot_colors_rbcl.csv")
colors

##### MAKE TAXAPLOT #####
# Assign taxa names to colors
scolors <- colors$plot_colors
names(scolors) <- colors$plot_names

##### ORDER THE TAXA SO OTHERS ARE AT THE BOTTOM #####
## order by decreasing relative abundance
alldata = alldata[order(-alldata$ra),]

## get list of factors in order
natural.genus.order = as.list(c(unique(alldata$plotnames)))

## remove others from list #!#
no.others=natural.genus.order[!natural.genus.order == 'Others']

## add Others to end of list
plot.order = append(no.others, "Others")

## set plot_names levels
plot.order = unlist(plot.order)

## order dataframe by relative abundance
alldata$plotnames = factor(alldata$plotnames, levels=c(plot.order))


##### MAKE THE TAXAPLOTS ######

## make the plot
mol_rbcl_plot <- ggplot(alldata, aes(x=as.character(Row.names), y=as.numeric(ra), 
                               fill=as.factor(plotnames)))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values=scolors)+
  guides(fill=guide_legend(ncol=1))+
  facet_nested(.~leaf_section + database, scales="free")+
  labs(y="Relative Abundance", x="Sample", fill="Genus  ") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 0))


ggarrange(myplot, mol_18S_plot, mol_rbcl_plot, labels = c("A", "B", "C"), ncol = 3)

#### make a common color legend ####
# Set up the plots
# Set up the color legend for the taxa plots
taxa_colors <- scale_fill_manual(values = scolors, name = "Taxa")

myplot <- ggplot(data, aes(x = box, y = relative_proportion, fill = genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = scolors) +
  labs(x = "Box", y = "Sum", fill = "Taxa") +
  facet_wrap(~blade_position, scales = "free_x") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Set up the 18S molecular taxa plot
mol_18S_plot <- ggplot(alldata, aes(x = as.character(Row.names), y = as.numeric(ra), fill = as.factor(plotnames))) +
  geom_bar(stat = "identity") +
  taxa_colors +
  facet_nested(. ~ leaf_section + database, scales = "free") +
  labs(y = "Relative Abundance", x = "Sample", fill = "Taxa") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 0))

# Set up the RBCL molecular taxa plot
mol_rbcl_plot <- ggplot(alldata, aes(x = as.character(Row.names), y = as.numeric(ra), fill = as.factor(plotnames))) +
  geom_bar(stat = "identity") +
  taxa_colors +
  facet_nested(. ~ leaf_section + database, scales = "free") +
  labs(y = "Relative Abundance", x = "Sample", fill = "Taxa") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 0))

# Arrange the plots vertically with a common color legend
combined_plot <- ggarrange(myplot, mol_18S_plot, mol_rbcl_plot, labels = c("A", "B", "C"), ncol = 1)

# Print the combined plot
combined_plot

