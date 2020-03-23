#example of how to utilize data from different sheets on excel file

#By Steve Formel
#Last updated March 12, 2020

library(readxl)

sample_key <- read_excel("GOMRI_R5x2860000002_data.xlsx", 
                                     sheet = "sample_key", na = "NA")


infl <- read_excel("data/S5_final_data_4Mar2020.xlsx", 
                         sheet = "inflorescences", na = "NA")

#joining data example----
library(dplyr)

df <- right_join(x = sample_key, y = infl, by = 'sampleID_stem')


#example plot

library(ggplot2)

ggplot(data = df,
       aes(x = num_infl,
           y = avg_mass_infl)) +
  geom_point(aes(color = oil_added)) +
  facet_wrap(~ sampling_period)

#Example of how to join multiple data sheets----

pmorph <- read_excel("GOMRI_R5x2860000002_data.xlsx", 
                   sheet = "plant_morphology", na = "NA")

df <- right_join(x = sample_key, y = infl) %>%
  right_join(x = pmorph, y = ., by = 'sampleID_stem')

#example plot

library(ggplot2)

ggplot(data = df,
       aes(x = num_stems_live,
           y = avg_mass_infl)) +
  geom_point(aes(color = oil_added)) +
  facet_wrap(~ sampling_period)

