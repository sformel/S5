#Plant biomass and stem count over time
#Last updated: August 30, 2022
#By Steve Formel

## Description
#Run PERMANOVA and make PCA of 16S based on Aitchison

#load libraries----

library(cowplot)
library(readxl)
library(tidyverse)
library(compositions)
library(vegan)
library(plotrix)
library(phyloseq)
library(data.table)
library(reshape2)
library(car)
library(compositions)
library(ggcorrplot)
library(ggforce)

#Read in data and clean-----

#Import ITS data and clean

#read in data
source("./import_and_clean_S5_ITS.R")

S5.exp <-
  prune_samples(sample_data(S5)$special_char == "Experimental Sample", S5)
S5.all <- S5

S5.ITS <- S5.exp

S5.rare <- rarefy_even_depth(physeq = S5.ITS,
                              sample.size = 3914)
S5.ITS <- S5.rare

## Import 16S phyloseq object and clean
#read in data

S5 <- readRDS("../Relevant Data/feature_tables_as_phyloseq_objects/S5_16S_pseq.rds")

taxa_names(S5) <- paste0("Seq", seq(ntaxa(S5)))

S5.noCL.noMT <- S5 %>%
  subset_taxa(Family != "Mitochondria" & Order != "Chloroplast")

S5.16S <- S5.noCL.noMT

S5 <- S5.16S

S5.rare <- rarefy_even_depth(physeq = S5,
                              sample.size = 10441)
S5 <- S5.rare

sort(sample_sums(S5))

#fix factors
sample_data(S5) <- sample_data(S5) %>%
  data.frame() %>%
  mutate(across(where(is.character), as.factor)) %>%
  dplyr::mutate(across(c(plantID, table, block), factor))

#Set up plotting objects

#Break out some elements for easy adjustment of plots----
annotation_size <- 5 / ggplot2::.pt
point_size <- 5 / ggplot2::.pt
stroke_size <- 0.5 / ggplot2::.pt

#color palette
cPAL <- c("#F0E442", "#0072B2")

S5_theme <- theme_bw(base_size = 7) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(7, 4, 2, 2, "mm"),
    axis.title.x = element_text(margin = margin(
      t = 2,
      r = 0,
      b = 0,
      l = 0,
      unit = "mm"
    )),
    axis.title.y = element_text(margin = margin(
      t = 0,
      r = 2,
      b = 0,
      l = 0,
      unit = "mm"
    ))
  )


SF1_aes <- aes(x = PC1,
               y = PC2,
               shape = orig_soil,
               fill = oil_added)

my.scale_aes <- list(geom_point(size = 2),
                     scale_shape_manual(values = c(21, 24)),
                     scale_fill_manual(values = c(cPAL, "black", "gray")))


#16S analysis ----


## Soil - Did the inoculation work for 16S?

S5.soil <- prune_samples(sample_data(S5)$tissue == "Soil", S5)
S5.soil <-
  prune_samples((
    sample_data(S5.soil)$plantID != "3" |
      sample_data(S5.soil)$sampling_period != "N17"
  ),
  S5.soil)


#takes about 5 min on SF macbook pro
df.soil <- otu_table(S5.soil) %>%
  as.data.frame()

df.meta.soil <- sample_data(S5.soil) %>%
  data.frame()

#make generic for downstream analysis
df <- df.soil
df.meta <- df.meta.soil

clr.df <- clr(df)
pca = prcomp(clr.df, center = FALSE, scale. = FALSE)

PCA.loadings = pca$rotation
PCA.scores = pca$x

#For below analyses
comm.df <- clr.df %>%
  as.data.frame()

#for PERMANOVA with all chemistry
comm.df.16S <- comm.df

#### Plot PCA - Aitchison

#Just SOIL

#get summary data
pca.sum <- as.data.frame(summary(pca)$importance)
df.PCA <- cbind(PCA.scores, df.meta)

#get rid of plant 3 at N17, crazy outlier

ggplot(data = df.PCA %>%
         filter(!(plantID == "3" & sampling_period == "N17")),
       aes(
         x = PC1,
         y = PC2,
         fill = orig_soil,
         shape = oil_added
       )) +
  my.scale_aes +
  ggforce::geom_mark_ellipse(aes(linetype = plant_trt),
                             alpha = 0.2) +
  scale_linetype_manual(labels = c("No Plant",
                                   "Planted"),
                        values = c("dashed", "solid")) +
  labs(
    title = NULL,
    caption = "Ellipses are visual guides and have no statistical meaning",
    shape = "Oil Addition",
    fill = "Inoculum",
    linetype = "Plant Treatment",
    x = paste("PC1 (", round(pca.sum$PC1[2] * 100, digits = 1), "%)"),
    y = paste("PC2 (", round(pca.sum$PC2[2] * 100, digits = 1), "%)")
  ) +
  facet_grid(
    cols = vars(sampling_period),
    labeller = labeller(
      plant_trt = c("Y" = "Planted", "N" = "No Plant"),
      sampling_period = c(
        "N16" = "6 Months",
        "J17" = "1 Year",
        "N17" = "1.5 Years",
        "J18" = "2 Years"
      )
    )
  ) +
  theme(legend.position = "bottom" , legend.direction = "vertical") +
  scale_color_manual(values = cPAL)


### NMDS

#plantID S3-N17 is a clear outlier

A <- comm.df[rownames(comm.df) != "S3-N17", ]
A <- metaMDS(comm = A, distance = "euclidean")

# cbind(df.meta[rownames(df.meta) != "S3-N17", ], A$points) %>%
#   ggplot(aes(
#     x = MDS1,
#     y = MDS2,
#     fill = orig_soil,
#     shape = oil_added
#   )) +
#   my.scale_aes +
#   ggforce::geom_mark_ellipse(aes(linetype = plant_trt),
#                              alpha = 0.2) +
#   scale_linetype_manual(labels = c("No Plant",
#                                    "Planted"),
#                         values = c("dashed", "solid")) +
#   labs(
#     title = NULL,
#     caption = "Ellipses are visual guides and have no statistical meaning",
#     shape = "Oil Addition",
#     fill = "Inoculum",
#     linetype = "Plant Treatment"
#   ) +
#   facet_grid(
#     cols = vars(sampling_period),
#     labeller = labeller(
#       plant_trt = c("Y" = "Planted", "N" = "No Plant"),
#       sampling_period = c(
#         "N16" = "6 Months",
#         "J17" = "1 Year",
#         "N17" = "1.5 Years",
#         "J18" = "2 Years"
#       )
#     )
#   ) +
#   theme(legend.position = "bottom" , legend.direction = "vertical") +
#   scale_color_manual(values = cPAL) +
#   geom_text(aes(label = plantID))

## NMDS Plot for publication

#The panels I made above were too much and there were some critiques of the PCA, so I've made something that is visually simpler and is an NMDS.  It's hard to see the yellow though, so I'm going to change the colors to something darker.


# p.16S.soil <- cbind(df.meta, A$points) %>%
#   group_by(plant_trt, oil_added, orig_soil, sampling_period) %>%
#   dplyr::summarise(MDS1.mean = mean(MDS1),
#                    MDS2.mean = mean(MDS2)) %>%
#   ggplot(aes(
#     x = MDS1.mean,
#     y = MDS2.mean,
#     group = interaction(plant_trt, oil_added, orig_soil)
#   )) +
#   geom_path(
#     aes(color = orig_soil,
#         linetype = fct_rev(plant_trt)),
#     arrow = arrow(type = "open",
#                   length = unit(2, "mm")),
#     size = 0.5,
#     show.legend = TRUE
#   ) +
#   geom_point(
#     data = . %>%
#       filter(sampling_period == "N16"),
#     aes(shape = oil_added,
#         fill = orig_soil),
#     size = 3
#   ) +
#   scale_shape_manual(values = c(21, 24)) +
#   scale_fill_manual(values = c(cPAL, "black", "gray")) +
#   scale_color_manual(values = c(cPAL, "black", "gray")) +
#   S5_theme +
#   theme(
#     legend.position = "bottom" ,
#     legend.direction = "vertical",
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank()
#   ) +
#   guides(
#     fill = guide_legend(override.aes = list(
#       shape = 22,
#       size = 4,
#       linetype = NULL
#     )),
#     color = guide_legend(override.aes = list(
#       shape = 22,
#       size = 4,
#       linetype = NULL
#     )),
#     shape = guide_legend(override.aes = list(
#       shape = c(21, 24),
#       size = 4,
#       linetype = NULL
#     ))
#   ) +
#   labs(
#     title = NULL,
#     caption = paste0("Stress = ", round(A$stress, 2)),
#     shape = "Oil Addition",
#     color = "Inoculum",
#     fill = "Inoculum",
#     linetype = "Plant Presence",
#     x = "MDS1",
#     y = "MDS2"
#   )

#This was extremely helpful for tweaking the plot
# nflplotR::ggpreview(plot = p.16S.soil,
#                     width = 85,
#                     height = 85,
#                     units = "mm",
#                     dpi = 300,
#                     scale = 1,
#                     device = "pdf")

# ggsave(plot = p.16S.soil,
#        filename = "S5_Figure6.pdf",
#        path = "figures/",
#        width = 85,
#        height = 85,
#        units = "mm",
#        dpi = 300,
#        scale = 1)
# 
# #For draft and co-author review
# ggsave(plot = p.16S.soil,
#        filename = "S5_Figure6.png",
#        path = "figures/",
#        width = 85,
#        height = 85,
#        units = "mm",
#        dpi = 300,
#        scale = 1)

#below I will combine this plot with the ITS analog

# Permanova - Aitchison

#Just SOIL

set.seed(1)
mod1 <- adonis(
  comm.df ~ plant_trt * orig_soil * oil_added,
  data = df.meta,
  permutations = 9999,
  method = "euclidean",
  strata = df.meta$sampling_period
)
mod1

#No change in statistical significance when rarefied to 10441

#check with dbRDA (actually RDA because it's euclidean, but the syntax was more obvious to me in dbRDA than in RDA)

mod1 <-
  dbrda(
    comm.df ~ plant_trt * orig_soil * oil_added + Condition(sampling_period),
    data = df.meta,
    distance = "euclidean"
  )

anova(mod1,
      by = "margin",
      permutations = 9999)


#Soil as points instead of centroids - this was preferred by reviewers over the vectors created above.

p.16S.soil.point <- cbind(df.meta, A$points) %>%
  ggplot(aes(
    x = MDS1,
    y = MDS2,
    group = interaction(plant_trt, oil_added, orig_soil)
  )) +
  geom_point(aes(shape = oil_added,
        fill = orig_soil),
    size = point_size,
    stroke = stroke_size) +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_manual(values = c(cPAL, "black", "gray")) +
  scale_color_manual(values = c(cPAL, "black", "gray")) +
  S5_theme +
  theme(
    legend.position = "right" ,
    legend.direction = "vertical",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.spacing.y = unit(annotation_size/10, 'pt'),
    legend.margin=margin(c(0,0,0,0))
  ) +
  coord_cartesian(xlim = c(-13, 13),
                  ylim = c(-13, 13)) +
  guides(
    fill = guide_legend(override.aes = list(
      shape = 22,
      size = point_size,
      linetype = NULL
    )),
    color = guide_legend(override.aes = list(
      shape = 22,
      size = point_size,
      linetype = NULL
    )),
    shape = guide_legend(override.aes = list(
      shape = c(21, 24),
      size = point_size,
      linetype = NULL
    ))
  ) +
  labs(
    title = "Prokaryote Soil Communities",
    caption = paste0("Stress = ", round(A$stress, 2)),
    shape = "Oil Addition",
    color = "Inoculum",
    fill = "Inoculum",
    linetype = "Plant Presence",
    x = "MDS1",
    y = "MDS2"
  ) +
  facet_grid(rows = vars(sampling_period), 
             cols = vars(plant_trt),
             labeller = labeller(sampling_period = c("N16" = "6 months",
                                                     "J17" = "1 year",
                                                     "N17" = "1.5 years",
                                                     "J18" = "2 years"),
                                 plant_trt = c("Y" = "Planted", 
                                               "N" = "No Plant")))

p.16S.soil.point

#This was extremely helpful for tweaking the plot
# nflplotR::ggpreview(plot = p.16S.soil.point,
#                     width = 85,
#                     height = 85,
#                     units = "mm",
#                     dpi = 300,
#                     scale = 1,
#                     device = "pdf")

ggsave(plot = p.16S.soil.point,
       filename = "S5_figure4.pdf",
       path = "figures/",
       width = 85,
       height = 85,
       units = "mm",
       dpi = 300,
       scale = 1)

#For draft and co-author review
ggsave(plot = p.16S.soil.point,
       filename = "S5_figure4.png",
       path = "figures/",
       width = 85,
       height = 85,
       units = "mm",
       dpi = 300,
       scale = 1)

#convert pdf to TIFF for publisher
pdftools::pdf_convert(
  pdf = "figures/S5_figure4.pdf",
  format = "tiff",
  filenames = "figures/S5_figure4.tif",
  dpi = 300
)

#ITS analysis----

S5 <- S5.ITS
### Plot of Inocula and Soil together
#takes about 5 min on SF macbook pro

df.all <- otu_table(S5) %>%
  as.data.frame()

df.meta.all <- sample_data(S5) %>%
  data.frame()

#make generic for downstream analysis
df <- df.all
df.meta <- df.meta.all

clr.df <- clr(df)
pca = prcomp(clr.df, center = FALSE, scale. = FALSE)

PCA.loadings = pca$rotation
PCA.scores = pca$x

#For below analyses
comm.df <- clr.df %>%
  as.data.frame()

#for PERMANOVA with all chemistry
comm.df.ITS <- comm.df

#### Plot PCA - Aitchison
#All tissues + Inocula

#make variable for inocula

df.meta$orig_soil[str_detect(df.meta$sampleID_stem, pattern = "O$")] <-
  "Prev-Oiled Inoc."
df.meta$orig_soil[str_detect(df.meta$sampleID_stem, pattern = "U$")] <-
  "Not Prev-Oiled Inoc"


#get summary data
pca.sum <- as.data.frame(summary(pca)$importance)
df.PCA <- cbind(PCA.scores, df.meta)

ggplot(data = df.PCA,
       aes(
         x = PC1,
         y = PC2,
         shape = orig_soil,
         fill = tissue
       )) +
  geom_point(size = 2) +
  scale_shape_manual(values = c(21, 24)) +
  theme_bw() +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 4)),
         shape = guide_legend(override.aes = list(size = 4))) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    title = "PCA on ITS Metagenome - Aitchison",
    caption = "Ellipses are visual guides and have no statistical meaning",
    x = paste("PC1 (", round(pca.sum$PC1[2] * 100, digits = 1), "%)"),
    y = paste("PC2 (", round(pca.sum$PC2[2] * 100, digits = 1), "%)"),
    shape = "Inoculum"
  ) +
  ggforce::geom_mark_ellipse(aes(color = tissue),
                             alpha = 0.1)


### Plot and PERMANOVA (adonis) of ITS - Aitchison

#By Tissue through time

#takes about 5 min on SF macbook pro
df.all <- otu_table(S5) %>%
  as.data.frame()

df.meta.all <- sample_data(S5) %>%
  data.frame()

#make generic for downstream analysis
df <- df.all
df.meta <- df.meta.all

clr.df <- clr(df)
pca = prcomp(clr.df, center = FALSE, scale. = FALSE)

PCA.loadings = pca$rotation
PCA.scores = pca$x

#For below analyses
comm.df <- clr.df %>%
  as.data.frame()

#for PERMANOVA with all chemistry
comm.df.ITS <- comm.df

#### Plot PCA - Aitchison
#All tissues


#get summary data
pca.sum <- as.data.frame(summary(pca)$importance)
df.PCA <- cbind(PCA.scores, df.meta)

ggplot(data = df.PCA,
       aes(
         x = PC1,
         y = PC2,
         shape = orig_soil,
         fill = tissue
       )) +
  geom_point(size = 2) +
  scale_shape_manual(values = c(21, 24)) +
  theme_bw() +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 4)),
         shape = guide_legend(override.aes = list(size = 4))) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    title = "PCA on ITS Metagenome - Aitchison",
    caption = "Ellipses are visual guides and have no statistical meaning",
    x = paste("PC1 (", round(pca.sum$PC1[2] * 100, digits = 1), "%)"),
    y = paste("PC2 (", round(pca.sum$PC2[2] * 100, digits = 1), "%)"),
    shape = "Inoculum"
  ) +
  facet_grid(cols = vars(sampling_period),
             rows = vars(oil_added)) +
  ggforce::geom_mark_ellipse(aes(linetype = orig_soil))


## Permanova - Aitchison

#Control for sampling period.  This means permutations will only happen within a sampling period.  It isn't exactly repeated-measures, but it's in th right d

#Tissue x Block
set.seed(1)
mod1 <- adonis(
  comm.df ~ tissue*block,
  data = df.meta,
  permutations = 9999,
  method = "euclidean",
  df.meta$sampling_period
)
mod1

B <- vegdist(comm.df, method = "euclidean")

anova(betadisper(B, group = df.meta$tissue))


#Tissue x Treatment
set.seed(1)
mod1 <- adonis(
  comm.df ~ tissue*orig_soil*oil_added,
  data = df.meta,
  permutations = 9999,
  method = "euclidean",
  strata = df.meta$sampling_period
)

mod1

B <- vegdist(comm.df, method = "euclidean")

anova(betadisper(B, group = df.meta$tissue))

## Soil - Did the inoculation work?

S5.soil <- prune_samples(sample_data(S5)$tissue == "Soil", S5)

#takes about 5 min on SF macbook pro
df.soil <- otu_table(S5.soil) %>%
  as.data.frame()

df.meta.soil <- sample_data(S5.soil) %>%
  data.frame()

#make generic for downstream analysis
df <- df.soil
df.meta <- df.meta.soil

clr.df <- clr(df)
pca = prcomp(clr.df, center = FALSE, scale. = FALSE)

PCA.loadings = pca$rotation
PCA.scores = pca$x

#For below analyses
comm.df <- clr.df %>%
  as.data.frame()

#for PERMANOVA with all chemistry
comm.df.ITS <- comm.df

#### Plot PCA - Aitchison

#Just SOIL

#get summary data
pca.sum <- as.data.frame(summary(pca)$importance)
df.PCA <- cbind(PCA.scores, df.meta)

ggplot(data = df.PCA,
       aes(
         x = PC1,
         y = PC2,
         fill = orig_soil,
         shape = oil_added
       )) +
  my.scale_aes +
  ggforce::geom_mark_ellipse(aes(linetype = plant_trt),
                             alpha = 0.2) +
  scale_linetype_manual(labels = c("No Plant",
                                   "Planted"),
                        values = c("dashed", "solid")) +
  labs(
    title = NULL,
    caption = "Ellipses are visual guides and have no statistical meaning",
    shape = "Oil Addition",
    fill = "Inoculum",
    linetype = "Plant Treatment",
    x = paste("PC1 (", round(pca.sum$PC1[2] * 100, digits = 1), "%)"),
    y = paste("PC2 (", round(pca.sum$PC2[2] * 100, digits = 1), "%)")
  ) +
  facet_grid(
    cols = vars(sampling_period),
    labeller = labeller(
      plant_trt = c("Y" = "Planted", "N" = "No Plant"),
      sampling_period = c(
        "N16" = "6 Months",
        "J17" = "1 Year",
        "N17" = "1.5 Years",
        "J18" = "2 Years"
      )
    )
  ) +
  theme(legend.position = "bottom" , legend.direction = "vertical") +
  scale_color_manual(values = cPAL)

### NMDS

A <- comm.df
A <- metaMDS(comm = A, distance = "euclidean")

cbind(df.meta, A$points) %>%
  ggplot(aes(
    x = MDS1,
    y = MDS2,
    fill = orig_soil,
    shape = oil_added
  )) +
  my.scale_aes +
  #ggforce::geom_mark_ellipse(aes(linetype = plant_trt),
  #                           alpha = 0.2) +
  scale_linetype_manual(labels = c("No Plant",
                                   "Planted"),
                        values = c("dashed", "solid")) +
  labs(
    title = NULL,
    caption = "Ellipses are visual guides and have no statistical meaning",
    shape = "Oil Addition",
    fill = "Inoculum",
    linetype = "Plant Treatment"
  ) +
  #   facet_grid(cols = vars(sampling_period),
  #              labeller = labeller(plant_trt = c("Y" = "Planted", "N" = "No Plant"),
  #                                  sampling_period = c("N16" = "6 Months", "J17" = "1 Year", "N17" = "1.5 Years", "J18" = "2 Years"))) +
  theme(legend.position = "bottom" , legend.direction = "vertical") +
  scale_color_manual(values = cPAL)

## Alternative NMDS Plot

#The panels I made above were too much and there were some critiques of the PCA, so I've made something that is visually simpler and is an NMDS.  It's hard to see the yellow though, so I'm going to change the colors to something darker.
# 
# p.ITS.soil <- cbind(df.meta, A$points) %>%
#   group_by(plant_trt, oil_added, orig_soil, sampling_period) %>%
#   dplyr::summarise(MDS1.mean = mean(MDS1),
#                    MDS2.mean = mean(MDS2)) %>%
#   ggplot(aes(
#     x = MDS1.mean,
#     y = MDS2.mean,
#     group = interaction(plant_trt, oil_added, orig_soil)
#   )) +
#   # geom_link2(aes(alpha = sampling_period),
#   #            size = 3) +
#   geom_path(
#     aes(color = orig_soil,
#         linetype = fct_rev(plant_trt)),
#     arrow = arrow(type = "open",
#                   length = unit(2, "mm")),
#     size = 0.5,
#     show.legend = TRUE
#   ) +
#   geom_point(
#     data = . %>%
#       filter(sampling_period == "N16"),
#     aes(shape = oil_added,
#         fill = orig_soil),
#     size = 3
#   ) +
#   scale_shape_manual(values = c(21, 24)) +
#   scale_fill_manual(values = c(cPAL, "black", "gray")) +
#   scale_color_manual(values = c(cPAL, "black", "gray")) +
#   S5_theme +
#   theme(
#     legend.position = "bottom" ,
#     legend.direction = "vertical",
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank()
#   ) +
#   guides(
#     fill = guide_legend(override.aes = list(
#       shape = 22,
#       size = 4,
#       linetype = NULL
#     )),
#     color = guide_legend(override.aes = list(
#       shape = 22,
#       size = 4,
#       linetype = NULL
#     )),
#     shape = guide_legend(override.aes = list(
#       shape = c(21, 24),
#       size = 4,
#       linetype = NULL
#     ))
#   ) +
#   labs(
#     title = NULL,
#     caption = paste0("Stress = ", round(A$stress, 2)),
#     shape = "Oil Addition",
#     color = "Inoculum",
#     fill = "Inoculum",
#     linetype = "Plant Presence",
#     x = "MDS1",
#     y = "MDS2"
#   )

#This was extremely helpful for tweaking the plot
# nflplotR::ggpreview(plot = p.ITS.soil,
#                     width = 85,
#                     height = 85,
#                     units = "mm",
#                     dpi = 300,
#                     scale = 1,
#                     device = "pdf")

# ggsave(filename = "S5_Figure5.pdf",
#        path = "figures/",
#        width = 85,
#        height = 85,
#        units = "mm",
#        dpi = 300,
#        scale = 1)
#
# #For draft and co-author review
# ggsave(filename = "S5_Figure5.png",
#        path = "figures/",
#        width = 85,
#        height = 85,
#        units = "mm",
#        dpi = 300,
#        scale = 1)
#

p.ITS.soil.point <- cbind(df.meta, A$points) %>%
  ggplot(aes(
    x = MDS1,
    y = MDS2,
    group = interaction(plant_trt, oil_added, orig_soil)
  )) +
  geom_point(aes(shape = oil_added,
        fill = orig_soil),
    size = point_size,
    stroke = stroke_size) +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_manual(values = c(cPAL, "black", "gray")) +
  scale_color_manual(values = c(cPAL, "black", "gray")) +
  S5_theme +
  theme(
    legend.position = "right" ,
    legend.direction = "vertical",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.spacing.y = unit(annotation_size/10, 'pt'),
    legend.margin=margin(c(0,0,0,0))
  ) +
  coord_cartesian(xlim = c(-15, 15),
                  ylim = c(-15, 15)) +
  guides(
    fill = guide_legend(override.aes = list(
      shape = 22,
      size = point_size,
      linetype = NULL
    )),
    color = guide_legend(override.aes = list(
      shape = 22,
      size = point_size,
      linetype = NULL
    )),
    shape = guide_legend(override.aes = list(
      shape = c(21, 24),
      size = point_size,
      linetype = NULL
    ))
  ) +
  labs(
    title = "Fungal Soil Communities",
    caption = paste0("Stress = ", round(A$stress, 2)),
    shape = "Oil Addition",
    color = "Inoculum",
    fill = "Inoculum",
    linetype = "Plant Presence",
    x = "MDS1",
    y = "MDS2"
  ) +
  facet_grid(rows = vars(sampling_period), 
             cols = vars(plant_trt),
             labeller = labeller(sampling_period = c("N16" = "6 months",
                                                     "J17" = "1 year",
                                                     "N17" = "1.5 years",
                                                     "J18" = "2 years"),
                                 plant_trt = c("Plant" = "Planted", 
                                               "No Plant" = "No Plant")))

p.ITS.soil.point

#This was extremely helpful for tweaking the plot
# nflplotR::ggpreview(plot = p.16S.soil.point,
#                     width = 85,
#                     height = 85,
#                     units = "mm",
#                     dpi = 300,
#                     scale = 1,
#                     device = "pdf")

ggsave(plot = p.ITS.soil.point,
       filename = "S5_figure5.pdf",
       path = "figures/",
       width = 85,
       height = 85,
       units = "mm",
       dpi = 300,
       scale = 1)

#For draft and co-author review
ggsave(plot = p.ITS.soil.point,
       filename = "S5_figure5.png",
       path = "figures/",
       width = 85,
       height = 85,
       units = "mm",
       dpi = 300,
       scale = 1)

#convert pdf to TIFF for publisher
pdftools::pdf_convert(
  pdf = "figures/S5_figure5.pdf",
  format = "tiff",
  filenames = "figures/S5_figure5.tif",
  dpi = 300
)

## Permanova - Aitchison

#Just SOIL

set.seed(1)

mod1 <- adonis(
  comm.df ~ plant_trt * orig_soil * oil_added,
  data = df.meta,
  permutations = 9999,
  method = "euclidean",
  strata = df.meta$sampling_period
)

mod1

#check with dbRDA (actually RDA because it's euclidean, but the syntax was more obvious to me in dbRDA than in RDA)

mod1 <-
  dbrda(
    comm.df ~ plant_trt + orig_soil + oil_added + Condition(sampling_period),
    data = df.meta,
    distance = "euclidean"
  )

anova(mod1,
      by = "margin",
      permutations = 9999)


#sampling period, plant_trt, orig_soil and the interaction of plant_trt:orig_soil are significant.  But it's still flimsy.  If I take sampling period out, it's still the same, results otherwise.  So I mine as well keep it in.

## dbRDA with oil components + CN

mod1 <-
  dbrda(
    comm.df ~ pH + conductivity + Hopanes + Total.Naph. + Total.Dibenz. + Total.Phenan. + Total.Chrys. + Condition(sampling_period + table + block),
    data = df.meta,
    distance = "euclidean"
  )

mod2 <- dbrda(comm.df ~ 1,
              data = df.meta,
              distance = "euclidean")

ordistep(mod2, scope = formula(mod1))

anova(mod1,
      by = "margin",
      permutations = 9999)

#Significant for Hopanes, pH, table, block.  After controlling for table and block, Hopanes are significant.  Oil addition wasn't significant in the permanova, so that's not very reassuring tha these results mean much.

## Roots - Did the inoculation work?

S5.roots <- prune_samples(sample_data(S5)$tissue == "Root", S5)

#takes about 5 min on SF macbook pro
df.roots <- otu_table(S5.roots) %>%
  as.data.frame()

df.meta.roots <- sample_data(S5.roots) %>%
  data.frame()

#make generic for downstream analysis
df <- df.roots
df.meta <- df.meta.roots

clr.df <- clr(df)
pca = prcomp(clr.df, center = FALSE, scale. = FALSE)

PCA.loadings = pca$rotation
PCA.scores = pca$x

#For below analyses
comm.df <- clr.df %>%
  as.data.frame()

#for PERMANOVA with all chemistry
comm.df.ITS <- comm.df

#### Plot PCA - Aitchison

#Just ROOTS

#get summary data
pca.sum <- as.data.frame(summary(pca)$importance)
df.PCA <- cbind(PCA.scores, df.meta)

#remove plant 60, J18 as an outlier

ggplot(data = df.PCA %>%
         filter(!(plantID == "60" & sampling_period == "J18")),
       aes(
         x = PC1,
         y = PC2,
         fill = orig_soil,
         shape = oil_added
       )) +
  my.scale_aes +
  ggforce::geom_mark_ellipse(aes(color = orig_soil),
                             alpha = 0.2) +
  labs(
    title = NULL,
    caption = "Ellipses are visual guides and have no statistical meaning",
    shape = "Oil Addition",
    fill = "Inoculum",
    color = "Inoculum",
    x = paste("PC1 (", round(pca.sum$PC1[2] * 100, digits = 1), "%)"),
    y = paste("PC2 (", round(pca.sum$PC2[2] * 100, digits = 1), "%)")
  ) +
  facet_grid(cols = vars(sampling_period),
             labeller = labeller(
               sampling_period = c(
                 "N16" = "6 Months",
                 "J17" = "1 Year",
                 "N17" = "1.5 Years",
                 "J18" = "2 Years"
               )
             )) +
  theme(legend.position = "bottom" , legend.direction = "vertical") +
  scale_color_manual(values = cPAL)

## NMDS

A <- comm.df[rownames(comm.df) != "S5_R60_J18", ]
A <- metaMDS(comm = A, distance = "euclidean")

cbind(df.meta[rownames(df.meta) != "S5_R60_J18", ], A$points) %>%
  ggplot(aes(
    x = MDS1,
    y = MDS2,
    fill = orig_soil,
    shape = oil_added
  )) +
  my.scale_aes +
  ggforce::geom_mark_ellipse(aes(linetype = plant_trt),
                             alpha = 0.2) +
  scale_linetype_manual(labels = c("No Plant",
                                   "Planted"),
                        values = c("dashed", "solid")) +
  labs(
    title = NULL,
    caption = "Ellipses are visual guides and have no statistical meaning",
    shape = "Oil Addition",
    fill = "Inoculum",
    linetype = "Plant Treatment"
  ) +
  facet_grid(cols = vars(sampling_period),
             labeller = labeller(
               sampling_period = c(
                 "N16" = "6 Months",
                 "J17" = "1 Year",
                 "N17" = "1.5 Years",
                 "J18" = "2 Years"
               )
             )) +
  theme(legend.position = "bottom" , legend.direction = "vertical") +
  scale_color_manual(values = cPAL)

## Alternative NMDS Plot

#The panels I made above were too much and there were some critiques of the PCA, so I've made something that is visually simpler and is an NMDS.  It's hard to see the yellow though, so I'm going to change the colors to something darker.


# p.ITS.roots <-
#   cbind(df.meta[rownames(df.meta) != "S5_R60_J18", ], A$points) %>%
#   group_by(oil_added, orig_soil, sampling_period) %>%
#   dplyr::summarise(MDS1.mean = mean(MDS1),
#                    MDS2.mean = mean(MDS2)) %>%
#   ggplot(aes(
#     x = MDS1.mean,
#     y = MDS2.mean,
#     group = interaction(oil_added, orig_soil)
#   )) +
#   # geom_link2(aes(alpha = sampling_period),
#   #            size = 3) +
#   geom_path(
#     aes(color = orig_soil),
#     arrow = arrow(type = "open",
#                   length = unit(2, "mm")),
#     size = 0.5,
#     show.legend = TRUE
#   ) +
#   geom_point(
#     data = . %>%
#       filter(sampling_period == "N16"),
#     aes(shape = oil_added,
#         fill = orig_soil),
#     size = 3
#   ) +
#   scale_shape_manual(values = c(21, 24)) +
#   scale_fill_manual(values = c(cPAL, "black", "gray")) +
#   scale_color_manual(values = c(cPAL, "black", "gray")) +
#   S5_theme +
#   theme(
#     legend.position = "bottom" ,
#     legend.direction = "vertical",
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank()
#   ) +
#   guides(
#     fill = guide_legend(override.aes = list(
#       shape = 22,
#       size = 4,
#       linetype = NULL
#     )),
#     color = guide_legend(override.aes = list(
#       shape = 22,
#       size = 4,
#       linetype = NULL
#     )),
#     shape = guide_legend(override.aes = list(
#       shape = c(21, 24),
#       size = 4,
#       linetype = NULL
#     ))
#   ) +
#   labs(
#     title = NULL,
#     caption = paste0("Stress = ", round(A$stress, 2)),
#     shape = "Oil Addition",
#     color = "Inoculum",
#     fill = "Inoculum",
#     x = "MDS1",
#     y = "MDS2"
#   )
# 
# 
#This was extremely helpful for tweaking the plot
# nflplotR::ggpreview(plot = p.ITS.roots,
#                     width = 85,
#                     height = 85,
#                     units = "mm",
#                     dpi = 300,
#                     scale = 1,
#                     device = "pdf")

# ggsave(filename = "S5_Figure5.pdf",
#        path = "figures/",
#        width = 85,
#        height = 85,
#        units = "mm",
#        dpi = 300,
#        scale = 1)
#
# #For draft and co-author review
# ggsave(filename = "S5_Figure5.png",
#        path = "figures/",
#        width = 85,
#        height = 85,
#        units = "mm",
#        dpi = 300,
#        scale = 1)
#


p.ITS.root.point <- cbind(df.meta[rownames(df.meta) != "S5_R60_J18", ], A$points) %>%
  ggplot(aes(
    x = MDS1,
    y = MDS2,
    group = interaction(oil_added, orig_soil)
  )) +
  geom_point(aes(shape = oil_added,
        fill = orig_soil),
    size = point_size,
    stroke = stroke_size) +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_manual(values = c(cPAL, "black", "gray")) +
  scale_color_manual(values = c(cPAL, "black", "gray")) +
  S5_theme +
  theme(
    legend.position = "right" ,
    legend.direction = "vertical",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.spacing.y = unit(annotation_size/10, 'pt'),
    legend.margin=margin(c(0,0,0,0))
  ) +
  coord_cartesian(xlim = c(-8, 8),
                  ylim = c(-8, 8)) +
  guides(
    fill = guide_legend(override.aes = list(
      shape = 22,
      size = point_size,
      linetype = NULL
    )),
    color = guide_legend(override.aes = list(
      shape = 22,
      size = point_size,
      linetype = NULL
    )),
    shape = guide_legend(override.aes = list(
      shape = c(21, 24),
      size = point_size,
      linetype = NULL
    ))
  ) +
  labs(
    title = "Fungal Root Communities",
    caption = paste0("Stress = ", round(A$stress, 2)),
    shape = "Oil Addition",
    color = "Inoculum",
    fill = "Inoculum",
    x = "MDS1",
    y = "MDS2"
  ) +
  facet_grid(rows = vars(sampling_period), 
             labeller = labeller(sampling_period = c("N16" = "6 months",
                                                     "J17" = "1 year",
                                                     "N17" = "1.5 years",
                                                     "J18" = "2 years")))

p.ITS.root.point

#This was extremely helpful for tweaking the plot
# nflplotR::ggpreview(plot = p.16S.soil.point,
#                     width = 85,
#                     height = 85,
#                     units = "mm",
#                     dpi = 300,
#                     scale = 1,
#                     device = "pdf")

ggsave(plot = p.ITS.root.point,
       filename = "S5_figure6.pdf",
       path = "figures/",
       width = 85,
       height = 85,
       units = "mm",
       dpi = 300,
       scale = 1)

#For draft and co-author review
ggsave(plot = p.ITS.root.point,
       filename = "S5_figure6.png",
       path = "figures/",
       width = 85,
       height = 85,
       units = "mm",
       dpi = 300,
       scale = 1)

#convert pdf to TIFF for publisher
pdftools::pdf_convert(
  pdf = "figures/S5_figure6.pdf",
  format = "tiff",
  filenames = "figures/S5_figure6.tif",
  dpi = 300
)

## Permanova - Aitchison

#Just ROOTS

set.seed(1)
mod1 <- adonis(
  comm.df ~ orig_soil * oil_added,
  data = df.meta,
  permutations = 9999,
  method = "euclidean",
  strata = df.meta$sampling_period
)
mod1

#check with dbRDA (actually RDA because it's euclidean, but the syntax was more obvious to me in dbRDA than in RDA)

mod1 <-
  dbrda(
    comm.df ~ orig_soil + oil_added + Condition(sampling_period),
    data = df.meta,
    distance = "euclidean"
  )

anova(mod1,
      by = "margin",
      permutations = 9999)

#Orig_soil, oil_added, orig_soil:oil_added. dbRDA backs up the main effects idea.

## dbRDA with oil components + CN

comm.df.sub <- subset(comm.df, rownames(comm.df) != "S5_R29_N17")
df.meta.sub <- df.meta %>%
  filter(plantID != 29 | sampling_period != "N17")

mod1 <-
  dbrda(
    comm.df.sub ~ percent_N + percent_C + C.N + table + block++stem_diam + num_nodes + num_stems_live + pH + conductivity + Hopanes + Total.Naph. + Total.Dibenz. + Total.Phenan. + Total.Chrys. + Condition(sampling_period),
    data = df.meta.sub,
    distance = "euclidean"
  )

mod2 <- dbrda(comm.df.sub ~ 1,
              data = df.meta.sub,
              distance = "euclidean")

ordistep(mod2, scope = formula(mod1))

anova(mod1,
      by = "margin",
      permutations = 9999)

#Block is the only significant factor.

## Leaves - Did the inoculation work?

#Just LEAVES

S5.leaf <- prune_samples(sample_data(S5)$tissue == "Leaf", S5)

#takes about 5 min on SF macbook pro
df.leaf <- otu_table(S5.leaf) %>%
  as.data.frame()

df.meta.leaf <- sample_data(S5.leaf) %>%
  data.frame()

#make generic for downstream analysis
df <- df.leaf
df.meta <- df.meta.leaf

clr.df <- clr(df)
pca = prcomp(clr.df, center = FALSE, scale. = FALSE)

PCA.loadings = pca$rotation
PCA.scores = pca$x

#For below analyses
comm.df <- clr.df %>%
  as.data.frame()

#for PERMANOVA with all chemistry
comm.df.ITS <- comm.df


#### Plot PCA - Aitchison

#Just LEAVES

#get summary data
pca.sum <- as.data.frame(summary(pca)$importance)
df.PCA <- cbind(PCA.scores, df.meta)

ggplot(data = df.PCA,
       aes(
         x = PC1,
         y = PC2,
         fill = orig_soil,
         shape = oil_added
       )) +
  my.scale_aes +
  ggforce::geom_mark_ellipse(aes(color = orig_soil),
                             alpha = 0.2) +
  labs(
    title = NULL,
    caption = "Ellipses are visual guides and have no statistical meaning",
    shape = "Oil Addition",
    fill = "Inoculum",
    color = "Inoculum",
    x = paste("PC1 (", round(pca.sum$PC1[2] * 100, digits = 1), "%)"),
    y = paste("PC2 (", round(pca.sum$PC2[2] * 100, digits = 1), "%)")
  ) +
  facet_grid(cols = vars(sampling_period),
             labeller = labeller(
               sampling_period = c(
                 "N16" = "6 Months",
                 "J17" = "1 Year",
                 "N17" = "1.5 Years",
                 "J18" = "2 Years"
               )
             )) +
  theme(legend.position = "bottom" , legend.direction = "vertical") +
  scale_color_manual(values = cPAL)

## NMDS

A <- comm.df
A <- metaMDS(comm = A, distance = "euclidean")

cbind(df.meta, A$points) %>%
  ggplot(aes(
    x = MDS1,
    y = MDS2,
    fill = orig_soil,
    shape = oil_added
  )) +
  my.scale_aes +
  ggforce::geom_mark_ellipse(aes(linetype = plant_trt),
                             alpha = 0.2) +
  scale_linetype_manual(labels = c("No Plant",
                                   "Planted"),
                        values = c("dashed", "solid")) +
  labs(
    title = NULL,
    caption = "Ellipses are visual guides and have no statistical meaning",
    shape = "Oil Addition",
    fill = "Inoculum",
    linetype = "Plant Treatment"
  ) +
  facet_grid(cols = vars(sampling_period),
             labeller = labeller(
               sampling_period = c(
                 "N16" = "6 Months",
                 "J17" = "1 Year",
                 "N17" = "1.5 Years",
                 "J18" = "2 Years"
               )
             )) +
  theme(legend.position = "bottom" , legend.direction = "vertical") +
  scale_color_manual(values = cPAL)

## Alternative NMDS Plot

#The panels I made above were too much and there were some critiques of the PCA, so I've made something that is visually simpler and is an NMDS.  It's hard to see the yellow though, so I'm going to change the colors to something darker.
# 
# p.ITS.leaves <- cbind(df.meta, A$points) %>%
#   group_by(oil_added, orig_soil, sampling_period) %>%
#   dplyr::summarise(MDS1.mean = mean(MDS1),
#                    MDS2.mean = mean(MDS2)) %>%
#   ggplot(aes(
#     x = MDS1.mean,
#     y = MDS2.mean,
#     group = interaction(oil_added, orig_soil)
#   )) +
#   geom_path(
#     aes(color = orig_soil),
#     arrow = arrow(type = "open",
#                   length = unit(2, "mm")),
#     size = 0.5,
#     show.legend = TRUE
#   ) +
#   geom_point(
#     data = . %>%
#       filter(sampling_period == "N16"),
#     aes(shape = oil_added,
#         fill = orig_soil),
#     size = 3
#   ) +
#   scale_shape_manual(values = c(21, 24)) +
#   scale_fill_manual(values = c(cPAL, "black", "gray")) +
#   scale_color_manual(values = c(cPAL, "black", "gray")) +
#   S5_theme +
#   theme(
#     legend.position = "bottom" ,
#     legend.direction = "vertical",
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank()
#   ) +
#   guides(
#     fill = guide_legend(override.aes = list(
#       shape = 22,
#       size = 4,
#       linetype = NULL
#     )),
#     color = guide_legend(override.aes = list(
#       shape = 22,
#       size = 4,
#       linetype = NULL
#     )),
#     shape = guide_legend(override.aes = list(
#       shape = c(21, 24),
#       size = 4,
#       linetype = NULL
#     ))
#   ) +
#   labs(
#     title = NULL,
#     caption = paste0("Stress = ", round(A$stress, 2)),
#     shape = "Oil Addition",
#     color = "Inoculum",
#     fill = "Inoculum",
#     x = "MDS1",
#     y = "MDS2"
#   )

#This was extremely helpful for tweaking the plot
# nflplotR::ggpreview(plot = p.ITS.leaves,
#                     width = 85,
#                     height = 85,
#                     units = "mm",
#                     dpi = 300,
#                     scale = 1,
#                     device = "pdf")

# ggsave(filename = "S5_Figure5.pdf",
#        path = "figures/",
#        width = 85,
#        height = 85,
#        units = "mm",
#        dpi = 300,
#        scale = 1)
#
# #For draft and co-author review
# ggsave(filename = "S5_Figure5.png",
#        path = "figures/",
#        width = 85,
#        height = 85,
#        units = "mm",
#        dpi = 300,
#        scale = 1)
#


p.ITS.leaf.point <- cbind(df.meta, A$points) %>%
  ggplot(aes(
    x = MDS1,
    y = MDS2,
    group = interaction(oil_added, orig_soil)
  )) +
  geom_point(aes(shape = oil_added,
        fill = orig_soil),
    size = point_size,
    stroke = stroke_size) +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_manual(values = c(cPAL, "black", "gray")) +
  scale_color_manual(values = c(cPAL, "black", "gray")) +
  S5_theme +
  theme(
    legend.position = "right" ,
    legend.direction = "vertical",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.spacing.y = unit(annotation_size/10, 'pt'),
    legend.margin=margin(c(0,0,0,0))
  ) +
  coord_cartesian(xlim = c(-6, 6),
                  ylim = c(-6, 6)) +
  guides(
    fill = guide_legend(override.aes = list(
      shape = 22,
      size = point_size,
      linetype = NULL
    )),
    color = guide_legend(override.aes = list(
      shape = 22,
      size = point_size,
      linetype = NULL
    )),
    shape = guide_legend(override.aes = list(
      shape = c(21, 24),
      size = point_size,
      linetype = NULL
    ))
  ) +
  labs(
    title = "Fungal Leaf Communities",
    caption = paste0("Stress = ", round(A$stress, 2)),
    shape = "Oil Addition",
    color = "Inoculum",
    fill = "Inoculum",
    x = "MDS1",
    y = "MDS2"
  ) +
  facet_grid(rows = vars(sampling_period), 
             labeller = labeller(sampling_period = c("N16" = "6 months",
                                                     "J17" = "1 year",
                                                     "N17" = "1.5 years",
                                                     "J18" = "2 years")))

p.ITS.leaf.point

#This was extremely helpful for tweaking the plot
# nflplotR::ggpreview(plot = p.16S.soil.point,
#                     width = 85,
#                     height = 85,
#                     units = "mm",
#                     dpi = 300,
#                     scale = 1,
#                     device = "pdf")

ggsave(plot = p.ITS.leaf.point,
       filename = "S5_figure7.pdf",
       path = "figures/",
       width = 85,
       height = 85,
       units = "mm",
       dpi = 300,
       scale = 1)

#For draft and co-author review
ggsave(plot = p.ITS.leaf.point,
       filename = "S5_figure7.png",
       path = "figures/",
       width = 85,
       height = 85,
       units = "mm",
       dpi = 300,
       scale = 1)

#convert pdf to TIFF for publisher
pdftools::pdf_convert(
  pdf = "figures/S5_figure7.pdf",
  format = "tiff",
  filenames = "figures/S5_figure7.tif",
  dpi = 300
)
## Permanova - Aitchison

#Just Leaves
set.seed(1)
mod1 <- adonis(
  comm.df ~ orig_soil * oil_added,
  data = df.meta,
  permutations = 9999,
  method = "euclidean",
  strata = df.meta$sampling_period
)
mod1

#check with dbRDA (actually RDA because it's euclidean, but the syntax was more obvious to me in dbRDA than in RDA)

mod1 <-
  dbrda(
    comm.df ~ orig_soil + oil_added + Condition(sampling_period),
    data = df.meta,
    distance = "euclidean"
  )

anova(mod1,
      by = "margin",
      permutations = 9999)

#orig_soil, and orig_soil:oil_added are the only levels that are significant. dbRDA supports main effect.

## dbRDA with variables but not treatments

comm.df.sub <- subset(comm.df, rownames(comm.df) != "S5_L08_N16")
df.meta.sub <- df.meta %>%
  filter(plantID != 8 | sampling_period != "N16")

mod1 <-
  dbrda(
    comm.df.sub ~ percent_N + percent_C + C.N + table + block++stem_ht + stem_diam + num_nodes + num_stems_live + pH + conductivity + Condition(sampling_period),
    data = df.meta.sub,
    distance = "euclidean"
  )

mod2 <- dbrda(comm.df.sub ~ 1,
              data = df.meta.sub,
              distance = "euclidean")

ordistep(mod2, scope = formula(mod1))

anova(mod1,
      by = "margin",
      permutations = 9999)

#The only significant thing here is table and block.

#Save plots -----

#two columns width = 180mm single column = 85mm

#soil

# soil.legend <- get_legend(
#   p.16S.soil +
#     theme(
#       legend.position = "bottom",
#       legend.box = "horizontal",
#       legend.margin = margin(
#         t = 1,
#         r = 1,
#         b = 1,
#         l = 1
#       ),
#       legend.title.align = 0.5,
#       legend.background = element_blank(),
#       legend.box.background = element_rect(colour = "lightgray",
#                                            size = stroke_size),
#       legend.key.width = unit(5, "mm")
#     )
# )
# 
# top <-
#   plot_grid(
#     p.16S.soil + S5_theme + theme(legend.position = "none"),
#     p.ITS.soil + S5_theme + theme(legend.position = "none"),
#     labels = c("Prokaryotes - Soil", "Fungi - Soil"),
#     hjust = c(-0.38, -0.71)
#   )
# 
# p <- plot_grid(top,
#                soil.legend,
#                ncol = 1,
#                rel_heights = c(1, 0.3))
# 
# # #This was extremely helpful for tweaking the plot
# # nflplotR::ggpreview(plot = p,
# #                     width = 180,
# #                     height = 85,
# #                     units = "mm",
# #                     dpi = 300,
# #                     scale = 1,
# #                     device = "pdf")
# 
# ggsave(
#   filename = "S5_figure4.pdf",
#   path = "figures/",
#   width = 180,
#   height = 85,
#   units = "mm",
#   dpi = 300,
#   scale = 1
# )
# 
# # #For draft and co-author review
# ggsave(
#   filename = "S5_figure4.png",
#   path = "figures/",
#   width = 180,
#   height = 85,
#   units = "mm",
#   dpi = 300,
#   scale = 1
# )
# 
# pdftools::pdf_convert(
#   pdf = "figures/S5_figure4.pdf",
#   format = "tiff",
#   filenames = "figures/S5_figure4.tif",
#   dpi = 300
# )
# 
# #roots
# 
# roots.legend <- get_legend(
#   p.ITS.roots +
#     theme(
#       legend.position = "bottom",
#       legend.box = "horizontal",
#       legend.margin = margin(
#         t = 1,
#         r = 1,
#         b = 1,
#         l = 1
#       ),
#       legend.title.align = 0.5,
#       legend.background = element_blank(),
#       legend.box.background = element_rect(colour = "lightgray",
#                                            size = stroke_size),
#       legend.key.width = unit(5, "mm")
#     )
# )
# 
# top <-
#   plot_grid(
#     p.ITS.roots + S5_theme + theme(legend.position = "none"),
#     p.ITS.leaves + S5_theme + theme(legend.position = "none"),
#     labels = c("Fungi - Roots", "Fungi - Leaves"),
#     hjust = c(-0.38, -0.71)
#   )
# 
# p <- plot_grid(top,
#                roots.legend,
#                ncol = 1,
#                rel_heights = c(1, 0.3))
# 
# # #This was extremely helpful for tweaking the plot
# nflplotR::ggpreview(plot = p,
#                     width = 180,
#                     height = 85,
#                     units = "mm",
#                     dpi = 300,
#                     scale = 1,
#                     device = "pdf")
# 
# ggsave(
#   filename = "S5_figure5.pdf",
#   path = "figures/",
#   width = 180,
#   height = 85,
#   units = "mm",
#   dpi = 300,
#   scale = 1
# )
# 
# # #For draft and co-author review
# ggsave(
#   filename = "S5_figure5.png",
#   path = "figures/",
#   width = 180,
#   height = 85,
#   units = "mm",
#   dpi = 300,
#   scale = 1
# )
# 
# pdftools::pdf_convert(
#   pdf = "figures/S5_figure5.pdf",
#   format = "tiff",
#   filenames = "figures/S5_figure5.tif",
#   dpi = 300
# )
