#Summarize how many samples remain for 16S after removing chloroplast/mitchondrial DNA
#Last updated 2022-12-16

library(readxl)
library(plyr)
library(dplyr)
library(phyloseq)
library(data.table)
library(reshape2)
library(car)
library(cowplot)
library(vegan)

## Import 16S phyloseq object and clean

#read in data

S5 <- readRDS("../Relevant Data/feature_tables_as_phyloseq_objects/S5_16S_pseq.rds")

taxa_names(S5) <- paste0("Seq", seq(ntaxa(S5)))

S5.noCL.noMT <- S5 %>%
  subset_taxa(Family != "Mitochondria" & Order != "Chloroplast")

S5.16S <- S5.noCL.noMT

sort(sample_sums(S5.16S))

S5 <- S5.16S

sample_data(S5)$seq_count <- sample_sums(S5)

#fix factors
sample_data(S5) <- sample_data(S5) %>%
  data.frame() %>%
  mutate(across(where(is.character), as.factor)) %>%
  dplyr::mutate(across(c(plantID, table, block), factor))



#Plot Shannon diversity for each tissue type by time point.
S5.rare <- rarefy_even_depth(physeq = S5,
                              sample.size = 1039)
S5 <- S5.rare

sort(sample_sums(S5))

S5SampleData <- data.frame(sample_data(S5))


Summary <- S5SampleData %>% 
  group_by(sampling_period, tissue) %>% 
  dplyr::summarise(n = n())

#write_csv(Summary, "Summary_16S_sample_count_after_filtering.csv")
