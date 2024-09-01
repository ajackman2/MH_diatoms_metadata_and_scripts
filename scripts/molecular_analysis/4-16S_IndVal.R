###### SET UP #######
## load packages
library(phyloseq)
library(tidyverse)
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


seagrass = readRDS("imput_data/MH_16S_filtered001percent_notrarefied_phyloseq.rds")

seagrass@sam_data$rd_filt = sample_sums(seagrass)

seagrass = subset_samples(seagrass, LeafSection !="NA")

seagrass = dephyloseq(seagrass) 

ind_data <- seagrass %>%
  spread(key = Genus, value = asv_abundance, fill = 0)

ind_data <- ind_data %>%
  select(-sample_id, -illumina_number, -library, -target, -PlantNumber, -LeafNumber, -SwabID, -PowerSoilDate, -PowerSoilID,  -FieldNotes,                       
         -LabNotes, -Aliquot_plate_well, -Normalized_aliquot_plate_well, -primer_18S_combo, -primer_RBCL_combo, -primer_16S_FW,                     
         -primer_16S_RV, -PCR_plate_position_aug31_and_sept1, -picogreen_well_id, -absorbance, -DNA_ng_C1_DNA_normalization,        
         -V2_DNA_normalization, -C2_DNA_normalization, -V1_DNA_normalization_uL, -vol_water_for_normalization_uL,    
         -DateCollected, -SampleTimeStart, -SwabTimeStart, -SwabTimeEnd, -unfiltered_read_depth, -rd_uf, -rd_filt,                           
         -ASVid, -Kingdom, -Phylum, -Class, -Order, -Family, -Species, -asv_id)

ind_data$Row.names <- as.factor(ind_data$Row.names)
ind_data$leaf_section <- as.factor(ind_data$LeafSection)

metaasv <- ind_data %>%
  dplyr::group_by(Row.names, LeafSection) %>%
  dplyr::summarize(across('37-13':'Yoonia-Loktanella')
  )
## Get otu table where samples are rows and ASVs are columns 
otu.tab = metaasv[,-c(1:2)]

## indval comparison 
substrate <- as.character(metaasv$LeafSection) 

## run indval
indval <- multipatt(otu.tab, substrate, duleg = TRUE, control = how(nperm=999))

# Get indval statistic
indval.16.str <- as.data.frame(indval$str)
indval.16.str$rn <- rownames(indval.16.str)

# get p-value
indval.stat <- as.data.frame(indval$sign) #get dataframe of indval statistic
indval.stat$rn <- rownames(indval.stat) # make column of ASVs

# specif as dataframe
indval.spe <- as.data.frame(indval$A)
# extract rownames into column
setDT(indval.spe, keep.rownames = TRUE)[] 
## rename columns 
colnames(indval.spe) <- paste0("spe.", colnames(indval.spe))
names(indval.spe)[names(indval.spe) == 'spe.rn'] <- 'rn'

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
spe.and.fid = full_join(indval.spe, indval.fid,
                        by="rn")
indval_table = full_join(str.and.stat, spe.and.fid,
                         by="rn")

## rename columns to join with taxonomy
names(indval_table)[names(indval_table) == 'rn'] <- 'ASV'

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

write.csv(corediff, "16S_core05_vs_07.csv")

###### PLOT RA OF DIFF TAXA ######
seagrass.df = dephyloseq(seagrass)

## only keep sig taxa
seagrass.df = inner_join(seagrass.df, coresig)

## relative abundance
seagrass.df$ra = as.numeric(seagrass.df$asv_abundance)/as.numeric(seagrass.df$rd_filt)

## remove 0s to have blanks
seagrass.df$ra = replace(seagrass.df$ra, seagrass.df$ra == 0, NA)

## replace index with letters
seagrass.df$index = gsub("1", "Base", seagrass.df$index)
seagrass.df$index = gsub("2", "Mid", seagrass.df$index)
seagrass.df$index = gsub("3", "Tip", seagrass.df$index)

## make bubble plot 
ggplot(seagrass.df, aes(x=as.character(Row.names), 
                        y=paste0(Genus," " , asv_id), 
                        size=ra, color=index))+
  geom_point()+
  facet_grid(Order~LeafSection, space="free", scale="free")+
  theme(strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x="Sample ID", y = "Genus and ASV ID", 
       size="Relative Abundance", color="Enriched On:")

ggsave("indval_16S.pdf", width = 16, height=7, units = "in")
