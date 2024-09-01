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

setwd("C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Seagrass_Montague_2021-03-07/galliano_seagrass_2021and2023/output/taxaplots")

seagrass = readRDS("C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Seagrass_Montague_2021-03-07/galliano_seagrass_2021and2023/imput_data/MH_RBCL_filtered001percent_notrarefied_phyloseq.RDS")
seagrass@sam_data$rd_filt = sample_sums(seagrass)
View(seagrass@sam_data)

blast_mol_RBCL <- read.csv("C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Seagrass_Montague_2021-03-07/galliano_seagrass_2021and2023/counting_data/05-15-24_blast_mol_RBCL")

##### REPLACE WITH UPDATED TAXONOMY ######
tax = blast_mol_RBCL[,c("ASVid", "Domain", "Kingdom","Phylum","Order","Genus","Species","asv_id")]

tax = ddply(tax, c("ASVid", "Domain", "Kingdom","Phylum","Order","Genus","Species","asv_id"),
            summarise,
            scount = length(asv_id))

tax = tax |> column_to_rownames(var="ASVid")

seagrass = phyloseq(sample_data(seagrass@sam_data),
                    tax_table(as.matrix(tax)),
                    otu_table(seagrass@otu_table, taxa_are_rows = T))




###### PREP FOR TAXAPLOT #####
## summarize at rank 6
seagrass = tax_glom(seagrass, taxrank = "Genus")

## get data out of phyloseq and inta a dataframe
seagrassdf = dephyloseq(seagrass)

## caluclate relative abundance of each rank 6 within each sample
seagrassdf$relativeabundance = as.numeric(seagrassdf$asv_abundance)/as.numeric(seagrassdf$rd_filt)

## make plotnames
seagrassdf$plotnames = paste0(seagrassdf$Order, ";", seagrassdf$Genus)

## make group for taxaplot
# make list of groups
grouplist = c(unique(seagrassdf$leafsection))

## summarize data by taxaplot geoup type. 
seagrass.sum = ddply(seagrassdf, c("leafsection", "plotnames"),
                     summarise,
                     sum = sum(relativeabundance))

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
    sample = subset(sorted, sorted$leafsection %in% c(j))
    
    ## get top 15 genera
    top = sample[c(1:15),]
    
    ## save list of top  abundance taxa
    t.tmp <- top
    top.df <- rbind.fill(top.df, t.tmp)
    
    ## close loop 
  }
}

## add identifier for top 15 taxa
top.df$place = "top_15"

###### FORMAT TOP TAXA OUTPUT #####
## join the top taxa and existing dataframe
alldata = full_join(seagrassdf, top.df)

## make the empty "place" cells say bottom. This workes because we used full_join
alldata$place = replace(alldata$place, is.na(alldata$place), "bottom")

## replace plot_names that have bottom taxa as their "place" with Other
alldata[alldata$place == "bottom",]$plotnames <- "Others"

###### GET COLORS FOR TAXAPLOT #####

# 1. find out how many colors you need
numcol <- length(unique(alldata$plotnames))

# 2. use a number seed to determine how qualpar samples your colors from its palette
set.seed(15)

# 3. use qualpalr colour palettes for easily distinguishing taxa
newpal <- qualpal(n = numcol, colorspace = "pretty")


# 4. Extract hex colors
hex = as.data.frame(newpal$hex)
colnames(hex) <- c("taxa_color")

# 5. Get list of taxa
tops = as.data.frame(c(unique(alldata$plotnames)))
colnames(tops) <- c("plotnames")

# 6. Join color list and taxa names
topcolors = cbind(tops, hex)

# 7. for the "others" plot name, replace that with grey 90 (this is just an astetic thing)
topcolors[topcolors$plotnames == "Others",]$taxa_color <- "grey90"

# 8. Make an object R can pull form for the colors
plotcolors <- topcolors$taxa_color
names(plotcolors) <- topcolors$plotnames


##### ORDER THE TAXA SO OTHERS ARE AT THE BOTTOM #####
## order by decreasing relative abundance
alldata = alldata[order(-alldata$relativeabundance),]

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
ggplot(alldata, aes(x=as.character(Row.names), y=as.numeric(relativeabundance), 
                    fill=as.factor(plotnames)))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values=plotcolors)+
  guides(fill=guide_legend(ncol=2))+
  facet_nested(.~leafsection+leafnumber, scales="free")+
  labs(y="Relative Abundance", x="Sample", fill="Taxa")

ggsave(filename = ("taxaplot_RBCL.pdf"), 
       width=17, height=8, units="in")


