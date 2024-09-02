# make a master file for the diatom data
# A. JAckman
# Nov 28 2023

# import libraries
library(dplyr)
library(stringr)

# set the working directory
setwd("/Users/andreajackman/R_Stuff/Parfrey/MH_diatoms_metadata_and_scripts/morphological_data/raw_morpho_counts")

# make a list of the csv files that are in this specific folder
csv_files <- list.files(pattern = "\\.csv$")
csv_files <- csv_files[1:15]

# make an empty dataframe to put the data into
combined_diatom_data <- data.frame()

# loop through the list, read the data, and append it to the dataframe
for (file in csv_files) {
  data <- read.csv(file) # read the file
  
  combined_diatom_data <- rbind(combined_diatom_data, data) # append the data using rbind
}

# make a new csv file with all of the data

# fix two rows that have errors for the total
combined_diatom_data[1466, 43] <- 20
combined_diatom_data[1467, 43] <- 32

# replaces the NAs with zeros
combined_diatom_data <- mutate_all(combined_diatom_data, ~ifelse(is.na(.), 0, .))

setwd("/Users/andreajackman/R_Stuff/Parfrey/MH_diatoms_metadata_and_scripts/morphological_data")
#write.csv(data, "master_diatom_data.csv", row.names = FALSE)

# reassign Fragilariopsis to Nitzschia because the ones counted originally as fragilariopsis are actually Nitzschia
data <- read.csv("master_diatom_data.csv")
data$Nitzschia <- data$Nitzschia + data$Fragilariopsis  # Add values in nitzschia and fragilariopsis columns
data$Fragilariopsis <- 0

# read the new csv
diatom_data <- read.csv("master_diatom_data.csv")
View(diatom_data)

# combine cocconeis_small and cocconeis_big because this is one genus
diatom_data$cocconeis <- diatom_data$cocconeis_small + diatom_data$cocconeis_big
diatom_data <- diatom_data[, -c(9, 10)]
diatom_data <- diatom_data[, c(1:40, 42, 41, (42:ncol(diatom_data)))]
diatom_data <- diatom_data[, -c(43)]

# select the columns to convert into upper case
columns_to_convert <- c('navicula', 'tabularia', 'cylindrotheca', 'pseudo_nitzschia', 'gomphonemopsis', 
                        'hyalodiscus', 'coscinodiscus', 'plagiotropis', 'nitzschia', 'rhoicosphenia', 'pseudogonphonema', 'thalassiosira', 
                        'eucampia', 'chaetoceros', 'diploneis', 'entomoneis', 'pariraphis_like', 'planothidium', 'halamphora',
                        'achnanthes', 'amphora', 'petroneis', 'bacillaria', 'ellerbeckia','minidiscus', 'fragilariopsis', 'rhizosolenia', 
                        'fogedia', 'skeletonema', 'pleurosigma', 'gomphoseptatum', 'licmophora', 'cocconeis')


# change the first letter to uppercase
colnames(diatom_data)[which(names(diatom_data) %in% columns_to_convert)] <- 
  sapply(columns_to_convert, function(col) paste(toupper(substr(col, 1, 1)), substr(col, 2, nchar(col)), sep = ""))

# rename columns so they are formatted correctly
diatom_data <- diatom_data %>%
  rename(Pseudonitzschia = Pseudo_nitzschia)

diatom_data <- diatom_data %>%
  rename(Pariraphis = Pariraphis_like)

diatom_data <- diatom_data %>%
  rename(Pseudogomphonema = Pseudogonphonema)

# add cols for all of the diatoms from the diversity survey
diatom_data$Actinoptychus <- NA
diatom_data$Attheya <- NA
diatom_data$Cyclotella <- NA
diatom_data$Dimeregramma <- NA
diatom_data$Ditylum <- NA
diatom_data$Donkinia <- NA
diatom_data$Entomoneis <- NA
diatom_data$Epithemia <- NA
diatom_data$Eupyxidicula <- NA
diatom_data$Encyonema <- NA
diatom_data$Fallacia <- NA
diatom_data$Grammatophora <- NA
diatom_data$Gyrosigma <- NA
diatom_data$Hanzschia <- NA
diatom_data$Haslea <- NA
diatom_data$Hobaniella <- NA
diatom_data$Leptocylindrus <- NA
diatom_data$Lyrella <- NA
diatom_data$Melosira <- NA
diatom_data$Odontella <- NA
diatom_data$Opephora <- NA
diatom_data$Paralia <- NA
diatom_data$Paribellus <- NA
diatom_data$Petroneis <- NA
diatom_data$Plagiogramma <- NA
diatom_data$Podosira <- NA
diatom_data$Psammodictyon <- NA
diatom_data$Rhabdonema <- NA
diatom_data$Rhopalodia <- NA
diatom_data$Trigonium <- NA
diatom_data$Trachyneis <- NA
diatom_data$Tryblionella <- NA
diatom_data$Undatella <- NA

# convert all of the NAs to zeros
diatom_data <- mutate_all(diatom_data, ~ifelse(is.na(.), 0, .))

# move total to the last column
diatom_data <- diatom_data[, c(setdiff(names(diatom_data), 'total'), 'total')]

# make the columns numeric values and the blade_position a factor
diatom_data <- diatom_data |>
  mutate_at(vars(Navicula:total), ~as.numeric(.)) 

#write.csv(diatom_data, "master_diatom_data.csv", row.names = FALSE)

#### fix the sample size for analysis (make sure all samples have equal sample size) ####
# need to correct the sample size of each file to n = 9 because some of them have more than n=9
# except 12AT29 which is 9
# most are 11 images but some are 10 or around 30
# will randomly pick which images to omit

# box 9 - has 30 images
set.seed(123)  # set seed
sample(1:30, 9, replace = FALSE)
# [1] 15 19 14  3 10 18 22 11  5

box_9 <- diatom_data[diatom_data$image_num %in% c(15, 19, 14, 3, 10, 18, 22, 11, 5) & diatom_data$box == 9, ]

# box 12AT10 - has 10 images
sample(1:10, 9, replace = FALSE)
# [1] 8  7  2  1  6  3  4 10  9
box_12A_T10 <- diatom_data[diatom_data$image_num %in% c(8, 7, 2, 1, 6, 3, 4, 10, 9) & diatom_data$box == '12A_T10' , ]

# box 12AT11 - has 10 images
sample(1:10, 9, replace = FALSE)
# [1] 7  5  4 10  2  9  3  1  8
box_12A_T11 <- diatom_data[diatom_data$image_num %in% c(7, 5, 4, 10, 2, 9, 3, 1, 8) & diatom_data$box == '12A_T11' , ]

# box 12AT12
sample(1:11, 9, replace = FALSE)
# [1] 10  7  5 11  9  3  1  2  8
box_12A_T12 <- diatom_data[diatom_data$image_num %in% c(10, 7,5, 11,  9,  3,  1,  2,  8) & diatom_data$box == '12A_T12' , ]

# box 12AT28
sample(1:10, 9, replace = FALSE)
# [1] 5  8  4  7 10  2  1  9  3
box_12A_T28 <- diatom_data[diatom_data$image_num %in% c(5,  8,  4,  7, 10,  2,  1,  9,  3) & diatom_data$box == '12A_T28' , ]

# 12A_T29
box_12A_T29 <- diatom_data[diatom_data$box == '12A_T29' , ]

# box 12AT30
sample(1:11, 9, replace = FALSE)
# [1] 9  6  5  7  1  2  4  3 10
box_12A_T30 <- diatom_data[diatom_data$image_num %in% c(9,  6,  5,  7,  1,  2,  4,  3, 10) & diatom_data$box == '12A_T30' , ]

# box 0_19_7
sample(1:30, 9, replace = FALSE)
# [1] 6 11  8 22 27  7 16 17 26
box_0_19_7 <- diatom_data[diatom_data$image_num %in% c(6, 11,  8, 22, 27,  7, 16, 17, 26) & diatom_data$box == '0_19_7' , ]

# box 0_1-1
sample(1:11, 9, replace = FALSE)
# [1] 2  1 11  4  5  7  3 10  9
box_0_1_1 <- diatom_data[diatom_data$image_num %in% c(2,  1, 11,  4,  5,  7,  3, 10,  9) & diatom_data$box == '0_1_1Pa' , ]

# box 0_8_3
sample(1:11, 9, replace = FALSE)
# [1] 9 11  7  3  6  4  1  2 10
box_0_8_3 <- diatom_data[diatom_data$image_num %in% c(9, 11,  7,  3,  6,  4,  1,  2, 10) & diatom_data$box == '0_8_3Mc' , ]

# box 0_9_17_6
sample(1:11, 9, replace = FALSE)
# [1] 3  7 11 10  6  2  5  1  9
box_9_17_6 <- diatom_data[diatom_data$image_num %in% c(3,  7, 11, 10,  6,  2,  5,  1,  9) & diatom_data$box == '9_17_6' , ]

# box 0_2_1
sample(1:11, 9, replace = FALSE)
# [1] 10  2 11  9  6  4  7  8  1
box_0_2_1 <- diatom_data[diatom_data$image_num %in% c(10,2, 11,  9,  6,  4,  7,  8,  1) & diatom_data$box == '0_2_1' , ]

# box 9_18D
sample(1:11, 9, replace = FALSE)
# [1] 6  3  8 11 10  9  1  7  4
box_9_18D <- diatom_data[diatom_data$image_num %in% c(6,  3,  8, 11, 10,  9,  1,  7,  4) & diatom_data$box == '9_18D_6' , ]

# box 0_3_1
sample(1:11, 9, replace = FALSE)
# [1] 7 11  6  9 10  3  2  1  5
box_0_3_1 <- diatom_data[diatom_data$image_num %in% c(7, 11,  6,  9, 10,  3,  2,  1,  5) & diatom_data$box == '0_3_1Da' , ]

# box 0_9_3
sample(1:11, 9, replace = FALSE)
# [1] 8  5  7  3  4  9  6  1 11
box_0_9_3 <- diatom_data[diatom_data$image_num %in% c(8,  5,  7,  3,  4,  9,  6,  1, 11) & diatom_data$box == '0_9_3' , ]

# box 0_20_7
sample(1:11, 9, replace = FALSE)
# [1] 9 7 6 2 1 8 5 3 4
box_0_20_7 <- diatom_data[diatom_data$image_num %in% c(9, 7, 6, 2, 1, 8, 5, 3, 4) & diatom_data$box == '0_20_7' , ]

# box 0_21_7
sample(1:11, 9, replace = FALSE)
# [1] 7  4  1  8  6 11 10  9  2
box_0_21_7 <- diatom_data[diatom_data$image_num %in% c(7,  4,  1,  8,  6, 11, 10,  9,  2) & diatom_data$box == '0_21_7' , ]

# bind these all together into a dataframe

# make an empty dataframe to put the data into
diatom_sample <- data.frame()
diatom_sample <- rbind(diatom_sample, box_9, box_12A_T10, box_12A_T11, box_12A_T12, box_12A_T28, box_12A_T29, box_12A_T30, box_0_19_7,
      box_0_1_1, box_0_8_3, box_9_17_6, box_0_2_1, box_9_18D, box_0_3_1, box_0_9_3, box_0_20_7, box_0_21_7)


# write_csv
#write.csv(diatom_sample, "corrected_master_diatom_data.csv", row.names = FALSE)

# read the new csv
diatom_data <- read.csv("corrected_master_diatom_data.csv")
View(diatom_data)

# combine cocconeis_small and cocconeis_big
# combine cocconeis_small and cocconeis_bug
diatom_data$cocconeis <- diatom_data$cocconeis_small + diatom_data$cocconeis_big
diatom_data <- diatom_data[, -c(9, 10)]
diatom_data <- diatom_data[, c(1:40, 42, 41, (42:ncol(diatom_data)))]
diatom_data <- diatom_data[, -c(43)]

# cols to change the names of
columns_to_convert <- c('navicula', 'tabularia', 'cylindrotheca', 'pseudo_nitzschia', 'gomphonemopsis', 
                        'hyalodiscus', 'coscinodiscus', 'plagiotropis', 'nitzschia', 'rhoicosphenia', 'pseudogonphonema', 'thalassiosira', 
                        'eucampia', 'chaetoceros', 'diploneis', 'entomoneis', 'pariraphis_like', 'planothidium', 'halamphora',
                        'achnanthes', 'amphora', 'petroneis', 'bacillaria', 'ellerbeckia','minidiscus', 'fragilariopsis', 'rhizosolenia', 
                        'fogedia', 'skeletonema', 'pleurosigma', 'gomphoseptatum', 'licmophora', 'cocconeis')


# change the first letter to uppercase
colnames(diatom_data)[which(names(diatom_data) %in% columns_to_convert)] <- 
  sapply(columns_to_convert, function(col) paste(toupper(substr(col, 1, 1)), substr(col, 2, nchar(col)), sep = ""))

# rename pariraphis_like and Pseudo_nitzschia
diatom_data <- diatom_data %>%
  rename(Pseudonitzschia = Pseudo_nitzschia)

diatom_data <- diatom_data %>%
  rename(Pariraphis = Pariraphis_like)

diatom_data <- diatom_data %>%
  rename(Pseudogomphonema = Pseudogonphonema)

# convert all of the NAs to zeros
diatom_data <- mutate_all(diatom_data, ~ifelse(is.na(.), 0, .))

# make the columns numeric values and the blade_position a factor
diatom_data <- diatom_data |>
  mutate_at(vars(Navicula:total), ~as.numeric(.)) 

#write.csv(diatom_data, "corrected_master_diatom_data.csv", row.names = FALSE)


