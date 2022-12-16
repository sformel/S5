#Plant biomass and stem count over time
#Last updated: May 14, 2022
#By Steve Formel

#load libraries----

library(cowplot)
library(readxl)
library(tidyverse)
library(compositions)
library(vegan)
library(plotrix)
library(brms)

#Break out some elements for easy adjustment of plots----
annotation_size <- 6 / ggplot2::.pt
point_size <- 5 / ggplot2::.pt
stroke_size <- 0.5 / ggplot2::.pt

S5_theme <- theme_bw(base_size = 7) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(2, 4, 2, 2, "mm"),
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


#color palette
cPAL <- c("#F0E442", "#0072B2")

#Read in data and clean

data.FP <-
  "../Relevant Data/data_available_from_GRIIDC/GOMRI_R5x2860000002_data.xlsx"

biomass <- read_excel(data.FP,
                      sheet = "biomass", na = "NA")

sample_key <- read_excel(data.FP,
                         sheet = "sample_key", na = "NA")

df <- right_join(x = sample_key, y = biomass, by = 'sampleID_stem')

#fix factors
df <- df %>%
  mutate(across(where(is.character), as.factor)) %>%
  dplyr::mutate(across(c(plantID, table, block), factor)) %>%
  filter(sampling_period == "J18", plant_trt == "Y")

#Gussy up the factor levels
df$oil_added <-
  plyr::revalue(df$oil_added, c("Y" = "Oil Added", "N" = "No Oil Added"))
df$orig_soil <-
  plyr::revalue(df$orig_soil,
                c("Y" = "Prev-Oiled Inoc.", "N" = "Not Prev-Oiled Inoc"))

## ANOVA on AG and BG biomass

m <- lm(sqrt(total_ag) ~ oil_added * orig_soil, data = df)

summary(m)
hist(resid(m))
shapiro.test(resid(m))

car::Anova(m, type ="III")

m <- lm(sqrt(total_bg) ~ oil_added * orig_soil, data = df)

summary(m)
hist(resid(m))
shapiro.test(resid(m))
car::Anova(m, type ="III")

#No difference in AG or BG biomass by oil_added*orig_soil, but must sqrt transform for normality.


#AG:BG ratio
df$AGBG_ratio <- df$total_ag/df$total_bg
m <- lm(sqrt(AGBG_ratio) ~ oil_added * orig_soil, data = df)

summary(m)
hist(resid(m))
shapiro.test(resid(m))
car::Anova(m, type ="III")

#Total Biomass
df$total_biomass <- df$total_ag + df$total_bg
m <- lm(sqrt(total_biomass) ~ oil_added * orig_soil, data = df)

summary(m)
hist(resid(m))
shapiro.test(resid(m))
car::Anova(m, type ="III")


## Plot
p <- df %>%
  group_by(oil_added, orig_soil) %>%
  summarise_at("total_live_mass",
               funs(mean,
                    sd,
                    se = sd(.) / sqrt(n()))) %>%
  ggplot(aes(
    x = orig_soil,
    y = mean,
    shape = oil_added,
    fill = orig_soil
  )) +
  geom_errorbar(
    aes(min = mean - 2 * se,
        max = mean + 2 * se),
    width = point_size / 6,
    size = stroke_size*1.5,
    position = position_dodge(width = point_size / 4)
  ) +
  geom_point(
    stroke = stroke_size,
    size = point_size,
    position = position_dodge(width = point_size / 4)
  ) +
  geom_point(data = df, 
             inherit.aes = FALSE, 
             aes(
    x = orig_soil,
    y = total_live_mass,
    shape = oil_added,
    ),
  alpha = 0.5,
  stroke = stroke_size,
    size = point_size*1.5,
    position = position_dodge(width = point_size / 4)
  ) +
  labs(
    y = "Mean Biomass (g)",
    x = NULL,
    fill = "Soil Inoculum",
    shape = "Oil Addition",
    caption = "Error bars represent \U00B1 2 SE"
  ) +
  S5_theme +
  labs(
    y = "Biomass (g)",
    x = "Soil Inoculum",
    fill = "Soil Inoculum",
    color = "Soil Inoculum",
    shape = "Oil Addition"
  ) +
  theme(legend.position = "right") +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.margin = margin(
      t = 1,
      r = 1,
      b = 1,
      l = 1,
      unit = 'mm'
    ),
    legend.title.align = 0.5,
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "lightgray",
                                         size = stroke_size),
    legend.key.width = unit(5, "mm")
  ) +
  expand_limits(x = 0,
                y = 0) +
  scale_fill_manual(values = cPAL) +
  scale_shape_manual(values = c(21, 24)) +
  guides(
    fill = guide_legend(
      title.position = "top",
      nrow = 2,
      override.aes = list(shape = 21,
                          size = point_size)
    ),
    shape = guide_legend(
      title.position = "top",
      nrow = 2,
      override.aes = list(size = point_size + 1)
    )
  )

p

#pull means to top layer
p$layers <- p$layers[c(3,1,2)]

p <- p +
  expand_limits(x = 0,
                y = 0)
#This was extremely helpful for tweaking the plot
# nflplotR::ggpreview(plot = p,
#                     width = 85,
#                     height = 85,
#                     units = "mm",
#                     dpi = 300,
#                     scale = 1,
#                     device = "pdf")

ggsave(
  filename = "S5_figure2.pdf",
  path = "figures/",
  width = 85,
  height = 85,
  units = "mm",
  dpi = 300,
  scale = 1
)

#For draft and co-author review
ggsave(
  filename = "S5_figure2.png",
  path = "figures/",
  width = 85,
  height = 85,
  units = "mm",
  dpi = 300,
  scale = 1
)

#convert pdf to TIFF for publisher
pdftools::pdf_convert(
  pdf = "figures/S5_figure2.pdf",
  format = "tiff",
  filenames = "figures/S5_figure2.tif",
  dpi = 300
)
