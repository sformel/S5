#Script for S5
#Purpose: Import ITS phyloseq object and replace/clean sample data

#Last updated 2022-12-16

#load libraries----
library(readxl)
library(plyr)
library(tidyverse)
library(phyloseq)

#Import phyloseq objects-----
#made in script called "ITS_seqs_BI.R" and moved into S5 data folder

S5 <- readRDS("../Relevant Data/feature_tables_as_phyloseq_objects/S5_ITS_ps.rds")

#Coordinate and replace sample data-----

#root oil data and soil CN data are incomplete

data.FP <- "../Relevant Data/data_available_from_GRIIDC/GOMRI_R5x2860000002_data.xlsx"

sample_key <- as.data.frame(read_excel(data.FP,
                                       sheet = "sample_key", na = "NA"))

microbial<- as.data.frame(read_excel(data.FP,
                                     sheet = "microbial", na = "NA"))

pmorph <- as.data.frame(read_excel(data.FP,
                        sheet = "plant_morphology", na = "NA"))

infl <- as.data.frame(read_excel(data.FP,
                                 sheet = "inflorescences", na = "NA"))

CN <- as.data.frame(read_excel(data.FP,
                               sheet = "CN", na = "NA"))

ph_cond <- as.data.frame(read_excel(data.FP,
                                    sheet = "pH_cond", na = "NA"))

oil <- as.data.frame(read_excel(data.FP,
                                sheet = "oil", na = "NA"))

#put it together
samdf <- right_join(x = sample_key, y = microbial) %>%
  right_join(x = ph_cond, y = .) %>%
  right_join(x = pmorph, y = ., by = 'sampleID_stem') %>%
  right_join(x = oil, y = .) %>%
  right_join(x = infl, y = .) %>%
  right_join(x = CN, y = .)

#reorder to match phyloseq object
samdf <- samdf[match(sample_data(S5)$sampleID, samdf$sampleID_ITS),]

#match?
any(sample_data(S5)$sampleID!=samdf$sampleID_ITS)


#cleaning and fixing----

#reorder factors
samdf$sampling_period <- factor(samdf$sampling_period, levels = c("NA", 
                                                                  "050216",
                                                                  "N16",
                                                                  "20317",
                                                                  "40117",
                                                                  "51817",
                                                                  "J17",
                                                                  "91817",
                                                                  "N17",
                                                                  "121317",
                                                                  "30318",
                                                                  "51418",
                                                                  "J18"))

#check structure
str(samdf)

#correct all character columns to factors
samdf <- samdf %>% 
  mutate_if(is.character,as.factor)

#correct table and block
samdf$table <- factor(samdf$table)
samdf$block <- factor(samdf$block)


#rename and reorder levels

samdf$special_char <- factor(samdf$special_char, levels = c("Pre-Experiment",
                                                            "Experimental Sample",
                                                            "Prev-Oiled Inoc.", 
                                                            "Not Prev-Oiled Inoc"))

samdf$oil_added <- plyr::revalue(samdf$oil_added, c("Y" = "Oil Added", "N" = "No Oil Added"))

samdf$plant_trt <- plyr::revalue(samdf$plant_trt, c("Y" = "Plant", "N" = "No Plant"))


samdf$orig_soil <- plyr::revalue(samdf$orig_soil, c("Y" = "Prev-Oiled Inoc.", "N" = "Not Prev-Oiled Inoc"))

samdf$Treatment <- interaction(samdf$oil_added, samdf$orig_soil,sep = " ; ")

#rename categories
colnames(samdf) <- samdf %>% 
  dplyr::rename("Hopanes" = hopanes, "Total Chrys." = total_methy_chry, "Total Naph." = total_methy_naph,"Total Dibenz." = total_methy_dibenz, "Total Phenan." = total_methy_phen) %>% 
  colnames()

samdf$pH <- round(as.numeric(samdf$pH),digits = 2)
samdf$conductivity <- round(as.numeric(samdf$conductivity),digits = 2)

#make select columns integers

cols <- c("num_nodes", "num_stems_live", "num_stems_outside", "num_infl")
samdf[cols] <- lapply(samdf[cols], as.integer)

#match rownames
rownames(samdf) <- rownames(sample_data(S5))

#double check that rownames matches sampleID
any(rownames(samdf)!=samdf$sampleID_ITS)

#get rid of non-essential columns for ITS analysis
samdf <- samdf %>%
  select(-c("DNA_ext_date_16S", "PCR1_16S", "PCR2_16S", "DNA_pur_conc_ITS", "date.x", "date.y", "date", "collected_by.x", "collected_by.y", "collected_by", "sampleID_16S"))
  
#put new data into pseq object
sample_data(S5) <- samdf

#get rid of everything but the S5 phyloseq object
rm(list=setdiff(ls(), "S5"))