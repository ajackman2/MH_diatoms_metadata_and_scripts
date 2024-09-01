###### SET UP #######
## load packages
library(phyloseq)
library(tidyverse)
library(plyr)
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


# read in the new data
ind_data <- read.csv("counting_data/06-05-24_blast_mol_RBCL")

# Spread the ASVid into columns
ind_data <- ind_data %>%
  spread(key = Genus, value = asv_abundance, fill = 0)

ind_data <- ind_data %>%
  select(-leafnumber, -ASVid, -Domain, -Phylum, -Class, -Order, -Empty, -Species, -asv_id, -X, -Kingdom, -Subkingdom, -subsp)

ind_data$Row.names <- as.factor(ind_data$Row.names)
ind_data$leafsection <- as.factor(ind_data$leafsection)

metaasv <- ind_data %>%
  dplyr::group_by(Row.names, leafsection) %>%
  dplyr::summarize(
    rd_uf = sum(rd_uf),
    rd_filt = sum(rd_filt),
    Achnanthales = sum(Achnanthales),
    Amphora = sum(Amphora),
    Attheya = sum(Attheya),
    Bacillariales = sum(Bacillariales),
    Bacillariophyceae = sum(Bacillariophyceae),
    Bacillariophyta = sum(Bacillariophyta),
    Chaetoceros = sum(Chaetoceros),
    Cocconeis = sum(Cocconeis),
    Cylindrotheca = sum(Cylindrotheca),
    Cymbellales = sum(Cymbellales),
    Diploneis = sum(Diploneis),
    Extubocellulus = sum(Extubocellulus),
    Fallacia = sum(Fallacia),
    Fragilariales = sum(Fragilariales),
    Fragilariophyceae = sum(Fragilariophyceae),
    Gedaniella = sum(Gedaniella),
    Gyrosigma = sum(Gyrosigma),
    Halamphora = sum(Halamphora),
    Haslea = sum(Haslea),  
    Hyalodiscus = sum(Hyalodiscus),
    Licmophora = sum(Licmophora),
    Minidiscus = sum(Minidiscus),
    Navicula = sum(Navicula),
    Naviculales = sum(Naviculales),
    Nitzschia = sum(Nitzschia),
    Odontella = sum(Odontella),
    Parlibellus = sum(Parlibellus),
    Plagiotropis = sum(Plagiotropis),
    Planothidium = sum(Planothidium),
    Pleurosigma = sum(Pleurosigma),
    Proschkinia = sum(Proschkinia),
    Pseudogomphonema = sum(Pseudogomphonema),
    Rhabdonema = sum(Rhabdonema),
    Seminavis = sum(Seminavis),
    Serratifera = sum(Serratifera),
    Skeletonema = sum(Skeletonema),
    Sternimirus = sum(Sternimirus),
    Tabularia = sum(Tabularia),
    Thalassiosira = sum(Thalassiosira),
    Thalassiosirales = sum(Thalassiosirales),
    Tryblionella = sum(Tryblionella)
  )


## Get otu table where samples are rows and ASVs are columns 
otu.tab = metaasv[,-c(1:4)]

## indval comparison 
substrate <- as.character(metaasv$leafsection) 

## run indval
indval <- multipatt(otu.tab, substrate, duleg = TRUE, control = how(nperm=999))

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

## rename columns to join with taxonomy
names(indval_table)[names(indval_table) == 'rn'] <- 'ASV'

coresig05 = subset(indval_table, indval_table$index != "NaN" & 
                     indval_table$p.value<=0.1 &
                     indval_table$stat>=0.5)

coresig = subset(indval_table, indval_table$index != "NaN" & 
                   indval_table$p.value<=0.1 &
                   indval_table$stat>=0.7)

#### make table of 0.5 vs 0.7 #####
coresig05$in05 = "yes"
coresig$in07 = "yes"

corediff = full_join(coresig05, coresig)

#write.csv(corediff, "18S_core05_vs_07.csv")

###### PLOT RA OF DIFF TAXA ######
seagrass.df = dephyloseq(seagrass)

## only keep sig taxa
seagrass.df = inner_join(seagrass.df, coresig)

## relative abundance
seagrass.df$ra = as.numeric(seagrass.df$asv_abundance)/as.numeric(seagrass.df$rd_filt)

## remove 0s to haveblanks
seagrass.df = subset(seagrass.df, seagrass.df$ra>0)

## replace index with letters
seagrass.df$index = gsub("1", "Base", seagrass.df$index)
seagrass.df$index = gsub("2", "Mid", seagrass.df$index)
seagrass.df$index = gsub("3", "Tip", seagrass.df$index)

## make bubble plot 
ggplot(seagrass.df, aes(x=as.character(Row.names), 
                        y=paste0(genus," " , asv_id), 
                        size=ra, color=index))+
  geom_point()+
  facet_grid(order~leaf_section, space="free", scale="free")+
  theme(strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x="Sample ID", y = "Genus and ASV ID", 
       size="Relative Abundance", color="Enriched On:")

ggsave("indval_18S.pdf", width = 16, height=7, units = "in")