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
ind_data <- read.csv("counting_data/06-11-24_blast_mol_18S")

# Spread the ASVid into columns
ind_data <- ind_data %>%
  spread(key = genus, value = asv_abundance, fill = 0)

ind_data <- ind_data %>%
  select(-leaf_number, -ASVid, -domain, -phylum, -class, -order, -family, -species, -asv_id, -X)

ind_data$Row.names <- as.factor(ind_data$Row.names)
ind_data$leaf_section <- as.factor(ind_data$leaf_section)

metaasv <- ind_data %>%
  dplyr::group_by(Row.names, leaf_section) %>%
  dplyr::summarize(
    rd_uf = sum(rd_uf),
    rd_filt = sum(rd_filt),
    Achnanthales = sum(Achnanthales),
    Achnanthes = sum(Achnanthes),
    Achnanthidium = sum(Achnanthidium),
    Amphora = sum(Amphora),
    Arcocellulus = sum(Arcocellulus),
    Bacillaria = sum(Bacillaria),
    Bacillariophyceae = sum(Bacillariophyceae),
    Bacillariophytina = sum(Bacillariophytina),
    Berkeleya = sum(Berkeleya),
    Cocconeis = sum(Cocconeis),
    Cylindrotheca = sum(Cylindrotheca),
    Diatomea = sum(Diatomea),
    Diploneis = sum(Diploneis),
    Fistulifera = sum(Fistulifera),
    Grammonema = sum(Grammonema),
    Gyrosigma = sum(Gyrosigma),
    Halamphora = sum(Halamphora),
    Haslea = sum(Haslea),
    Hyalodiscus = sum(Hyalodiscus),
    Licmophora = sum(Licmophora),
    Melosira = sum(Melosira),  
    Navicula = sum(Navicula),
    Nitzschia = sum(Nitzschia),
    Odontella = sum(Odontella),
    Phaeodactylum = sum(Phaeodactylum),
    Pinnularia = sum(Pinnularia),
    Placoneis = sum(Placoneis),
    Planothidium = sum(Planothidium),
    Pleurosigma = sum(Pleurosigma),
    Skeletonema = sum(Skeletonema),
    Stauroneidaceae = sum(Stauroneidaceae),
    Stauroneis = sum(Stauroneis),
    Stephanopyxis = sum(Stephanopyxis),
    Tabularia = sum(Tabularia),
    Thalassiosira = sum(Thalassiosira),
    Tryblionella = sum(Tryblionella)
  )



## Get otu table where samples are rows and ASVs are columns 
otu.tab = metaasv[,-c(1:4)]

## indval comparison 
substrate <- as.character(metaasv$leaf_section) 

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

#### everything else to the genus level ####
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


seagrass = readRDS("imput_data/MH_18S_filtered001percent_notrarefied_phyloseq.rds")

seagrass@sam_data$rd_filt = sample_sums(seagrass)


seagrass = dephyloseq(seagrass) 

seagrass = subset_samples(seagrass, leaf_section !="NA")

seagrass <- seagrass |>
  filter(class != 'Diatomea')

ind_data <- seagrass %>%
  spread(key = genus, value = asv_abundance, fill = 0)

ind_data <- ind_data %>%
  select(-leaf_number, -ASVid, -domain, -phylum, -class, -order, -family, -species, -asv_id)

ind_data$Row.names <- as.factor(ind_data$Row.names)
ind_data$leaf_section <- as.factor(ind_data$leaf_section)
ind_data$rd_uf <- as.numeric(ind_data$rd_uf)
ind_data$rd_filt <- as.numeric(ind_data$rd_filt)

metaasv <- ind_data %>%
  dplyr::group_by(Row.names, leaf_section) %>%
  dplyr::summarize(
    rd_uf = sum(rd_uf),
    rd_filt = sum(rd_filt),
    Allium = sum(Allium),
    Arthropoda = sum(Arthropoda),
    Aureobasidium = sum(Aureobasidium),
    Biecheleria = sum(Biecheleria),
    Brooklynella = sum(Brooklynella),
    Caenogastropoda = sum(Caenogastropoda),
    Cercozoa = sum(Cercozoa),
    Chytridium = sum(Chytridium),
    Cryomyces = sum(Cryomyces),
    Cyclopoida = sum(Cyclopoida),
    Dinophyceae = sum(Dinophyceae),
    Dysteria = sum(Dysteria),
    Ectocarpales = sum(Ectocarpales),
    Ectocarpales_fa = sum(Ectocarpales_fa),
    Eukaryota = sum(Eukaryota),
    Flamella = sum(Flamella),
    Halothrix = sum(Halothrix),  
    Heterosigma = sum(Heterosigma),
    Intramacronucleata = sum(Intramacronucleata),
    Leucocryptos = sum(Leucocryptos),
    Lobulomycetaceae = sum(Lobulomycetaceae),
    Mesophyllum = sum(Mesophyllum),
    Mollusca = sum(Mollusca),
    Myoida = sum(Myoida),
    Nematoda = sum(Nematoda),
    Nicotiana = sum(Nicotiana),
    Ochrophyta = sum(Ochrophyta),
    Omegastrombidium = sum(Omegastrombidium),
    Pelagodinium = sum(Pelagodinium),
    Penicillium = sum(Penicillium),
    Peridiniphycidae = sum(Peridiniphycidae),
    Peronosporomycetes = sum(Peronosporomycetes),
    Peronosporomycetes_fa = sum(Peronosporomycetes_fa),
    Peziza = sum(Peziza),
    Phyllodocida = sum(Phyllodocida),
    Phyllopharyngea = sum(Phyllopharyngea),
    Poales_fa = sum(Poales_fa),
    Pseudopirsonia = sum(Pseudopirsonia),
    Pylaiella = sum(Pylaiella),
    Rhogostoma = sum(Rhogostoma),
    Saccharomyces = sum(Saccharomyces),
    Scytosiphon = sum(Scytosiphon),
    Spirotrichea = sum(Spirotrichea),
    Stichococcus = sum(Stichococcus),
    Syndiniales_Group_I = sum(Syndiniales_Group_I),
    Taiwania = sum(Taiwania),
    Telonema = sum(Telonema),
    Terebellida = sum(Terebellida),
    Trebouxiophyceae = sum(Trebouxiophyceae),
    Ulotrichales_fa = sum(Ulotrichales_fa),
    Ulva = sum(Ulva),
    Ulvales_fa = sum(Ulvales_fa),
    Ulvophyceae = sum(Ulvophyceae),
    Veneroida = sum(Veneroida)
  )

## Get otu table where samples are rows and ASVs are columns 
otu.tab = metaasv[,-c(1:4)]

## indval comparison 
substrate <- as.character(metaasv$leaf_section) 

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

## merge with taxonomy
#indval_table= inner_join(indval_table, tax)

coresig05 = subset(indval_table, indval_table$index != "NaN" & 
                     indval_table$p.value<=0.05 &
                     indval_table$stat>=0.5)

coresig = subset(indval_table, indval_table$index != "NaN" & 
                   indval_table$p.value<=0.05 &
                   indval_table$stat>=0.7)

#### make table of 0.5 vs 0.7 #####
coresig05$in05 = "yes"
coresig$in07 = "yes"

corediff = full_join(coresig05, coresig)

write.csv(corediff, "18S_core05_vs_07.csv")

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

