


# STEP 1: Extraction of open-source data

#set your working directory
#identify the species you wish to include, here, 20 randomly selected species from the original study have been selected as example data.
#save these species' names as a vector

example_species <- c("Acacia paraguariensis", 
                     "Acosmium cardenasii",
                     "Aeschynomene denticulata",
                     "Arachis chiquitana",
                     "Bauhinia tuichiensis",
                     "Blepharodon philibertioides",
                     "Calea dalyi",
                     "Cascaronia astragalina",
                     "Eugenia boliviana",
                     "Fosterella vasquezii",
                     "Handroanthus selachidentatus",
                     "Inga approximata",
                     "Manihot fabianae",
                     "Mitracarpus bicrucis",
                     "Myrcia myrcioides",
                     "Nymphoides herzogii",
                     "Oxypetalum fuscum",
                     "Pectis harryi",
                     "Pitcairnia chiquitana",
                     "Senna coimbrae")

#Extract species occurrence records from BIEN using the R package 'BIEN' Maitner et al. (2017)
install.packages("BIEN")
library(BIEN)

occ_bien <- BIEN_occurrence_species(example_species, 
                                    cultivated = TRUE,
                                    only.new.world = FALSE,
                                    all.taxonomy = TRUE,
                                    native.status = TRUE,
                                    observation.type = TRUE,
                                    political.boundaries = TRUE)

#The code above retrieved 322 occurrence records for 35 variables (on 20th July 2020), however, this may change due to the continous updating of the BIEN database.

#save as a dataframe
occ_bien_df <- as.data.frame(occ_bien)

#now save to your working directory as Tab Separated Values (.tsv)
write.table(x = occ_bien_df, "occurrences_bien.tsv", sep='\t',row.names = F)

#Now we will extract species occurrence records from the GBIF database using the 'rgbif' package (Chamberlain & Boettiger, 2017)

install.packages("rgbif")
library(rgbif)
library(ape)
library(spocc)
library(geiger)
library(dismo)
library(geiger)
library(maptools)
library(XML)
library(ggmap)
library(sp)
library(rworldmap)
library(maps)
library(rgeos)
library(mapr)
library(ggplot2)
library(RColorBrewer)

#it is important to use the same species' names vector as we did for the BIEN extraction.
#we know that there will be approximately the same number of records from GBIF, so a limit of 1000 records per species will be sufficient. 
occ_gbif <- occ(example_species, from="gbif", limit=1000, has_coords=TRUE) 

#to download all associated metadata, use the following code. This is important for assessing which records are errorneous or not.
occ_gbif_meta <- occ2df(occ_gbif$gbif, what="all")
occ_gbif_matrix <- as.matrix(occ_gbif_meta)

#on 20th July 2020, this downloaded 512 records with 132 variables, again this may change due to the continous updating of the GBIF database.

#save GBIF data to your working directory
write.table(x = occ_gbif_matrix, "occurrences_gbif.tsv", sep='\t',row.names = F)




# STEP 2: Taxonomic standardization and Basic Coordinate Check

#This step is completed in Microsoft Excel, using the example data, zero BIEN records were removed and 10 removed from GBIF (3 during taxonomic standardization and 7 during Basic Coordinate Check; see methods section in manuscript for more information).
#create a new column headed 'standardized_names' where synonyms included in the open-source data sets will be standardized to their accepted names.
#save your new csv files as standardized_gbif.csv and standardized_bien.csv in your working directory.



# STEP 3: Automated Clean

#For this stage, we will use the automated coordinate cleaning R package 'CoordinateCleaner' (Zizka et al. 2019)
install.packages("tibble")
install.packages("CoordinateCleaner")
install.packages("tidyverse")
library(CoordinateCleaner)
library(tidyverse)
library(tibble)

stn_bien <- read.csv("standardized_names_bien.csv", header = TRUE) #322 records of 37 variables
stn_gbif <- read.csv("standardized_names_gbif.csv", header = TRUE) #499 records of 134 variables

flags_bien <- clean_coordinates(x = stn_bien, 
                                lon = "longitude", 
                                lat = "latitude", 
                                species = "standardized_names", 
                                countries = NULL,
                                tests = c("capitals",
                                          "centroids",
                                          "equal",
                                          "gbif",
                                          "institutions",
                                          "outliers",
                                          "seas",
                                          "zeros"))

summary(flags_bien) #produces an overview of failed tests and how many records failed each one.

flags_gbif <- clean_coordinates(x = stn_gbif, 
                                lon = "longitude", 
                                lat = "latitude", 
                                species = "standardized_names", 
                                countries = NULL,
                                tests = c("capitals",
                                          "centroids",
                                          "equal",
                                          "gbif",
                                          "institutions",
                                          "outliers",
                                          "seas",
                                          "zeros"))

summary(flags_gbif) #and the same for GBIF

plot(flags_gbif, lon = "longitude", lat = "latitude") #we can visualize the flagged and clean records from each data set.
plot(flags_bien, lon = "longitude", lat = "latitude") #as you can see there is a large geographic spread of erroneous records for BIEN.

bien.data.clean <- occ_bien[flags_bien$.summary,] #create a new dataframe excluding the flagged records, these records will be used during the manual clean.
gbif.data.clean <- occ_gbif_matrix[flags_gbif$.summary,] #do the same for GBIF.

#From this point onwards, I will only use the BIEN data as an example as the process is exactly the same for the GBIF data.

bien_out <- cbind(occ_bien, flags_bien)
bien.flagged <- bien_out[!bien_out$.summary,] #now attach the flag reasons to each flagged record from the BIEN data set.

write.csv(x = bien.flagged, "bien_flags.csv", sep='\t',row.names = F) #save to your working directory, there should be 19 flagged records with their flag reasons attached on the end.

write.csv(x = bien.data.clean, "bien_auto_clean.csv", sep='\t',row.names = F) #for cleaned BIEN records, there should be 303.


#STEP 4: Manual Clean

#This stage was not performed in R, instead a GIS was used to plot the automated cleaned records and outliers were removed by eye. See the methods section in the accompanying manuscript for further details.


#STEP 5: Computation of threat categories using rCAT (Moat )

install.packages("sp")
install.packages("rCAT")
install.packages("tidyverse")
library(sp)
library(rCAT)
library(tidyverse)

#For the purposes of this example code, we will just compute EOO metrics for the automated clean BIEN data using the ConBatch() function.

bien_auto_EOOs <- ConBatch(bien.data.clean$scrubbed_species_binomial,
                           bien.data.clean$latitude,
                           bien.data.clean$longitude, 
                           2000) #make sure to include the cell size (2000m)

write.csv(x = bien_auto_EOOs, "bien_EOOS.csv", sep='\t', row.names = F)

#Now we have the computed EOO, AOO (not needed for this study) and associated threat categories for our example species that have been processed through the automated clean.
#This final process was performed for 1) uncleaned GBIF and BIEN data, 2) automated GBIF and BIEN data and 3) manually cleaned GBIF and BIEN data.
