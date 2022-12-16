#Map of Soil Inocula Sampling Area
#Last updated: 2022-12-16

#load libraries----

library(ggmap)
library(ggsn)
library(dplyr)
library(sf)
library(cowplot)

#Break out some elements for easy adjustment of plots
annotation_size <- 6 / ggplot2::.pt

S5_theme <- theme_bw(base_size = 8) +
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
#Site plot-----

#Read in Custom shapefile made in Google Earth and ArcGIS
shapefile.site <-
  st_read("../Relevant Data/shapefiles/porkchop_island/S5_Porkchop_Island-polygon.shp")
#make data frame of annotations
annotations.site <- data.frame(
  x = c(-89.891,-89.887),
  y = c(29.444, 29.439),
  name = c("Heavily Oiled",
           "Reference")
)
#plot
site.map <- ggplot() +
  geom_sf(
    data = shapefile.site,
    color = 'gray',
    fill = '#6ed471',
    alpha = 0.7,
    size = .2
  ) +
  blank() +
  north(shapefile.site,
        scale = 0.15,
        anchor = c("x" = -89.889,
                   "y" = 29.437)) +
  ggsn::scalebar(
    data = shapefile.site,
    location = "bottomleft",
    dist = 0.25,
    dist_unit = "km",
    transform = TRUE,
    st.size = annotation_size,
    border.size = 0.2,
    st.dist = 0.05
  ) +
  labs(x = "Longitude (DD)",
       y = "Latitude (DD)")  +
  annotate(
    "text",
    x = annotations.site$x,
    y = annotations.site$y,
    label = annotations.site$name,
    hjust = 0,
    fontface = 2,
    size = annotation_size
  ) +
  S5_theme

#Inset of Southern Louisiana----

shapefile.state <-
  st_read("../Relevant Data/shapefiles/states/states.shp") %>%
  filter(STATE_FIPS == 22) %>%
  st_crop(., y = c(
    xmin = -92,
    ymin = 28,
    xmax = -89,
    ymax = 31
  ))

#make data frame of annotations
annotations.state <- data.frame(
  x = c(-90.1,-90.1),
  y = c(29.95, 29.42),
  name = c("New Orleans",
           "Bay Jimmy")
)

#plot
inset <- ggplot() +
  geom_sf(data = shapefile.state,
          fill = "white",
          size = 0.2) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "#f7f8fa"),
    panel.border = element_rect(colour = "black",
                                fill = NA)
  ) +
  geom_text(
    data = annotations.state[1, ],
    aes(x = x,
        y = y,
        label = name),
    hjust = 1,
    fontface = 2,
    size = annotation_size
  ) +
  geom_point(data = annotations.state[1, ],
             #decided to only include New Orleans
             aes(x = x + 0.1,
                 y = y),
             size = 1 / ggplot2::.pt)


#Make Final Plot and Save----

arrow_points <- data.frame(x1 = 0.84,
                           x2 = 0.68,
                           y1 = 0.67,
                           y2 = 0.48)

p <- ggdraw(site.map) +
  draw_plot(
    inset,
    x = 0.64,
    y = 0.58,
    width = 0.3,
    height = 0.3
  ) +
  geom_segment(
    data = arrow_points,
    aes(
      x = x1,
      y = y1,
      xend = x2,
      yend = y2
    ),
    arrow = arrow(length = unit(2, "mm")),
    lineend = "round"
  )

p

#This was extremely helpful for tweaking the plot
# nflplotR::ggpreview(plot = p,
#                     width = 85,
#                     height = 85,
#                     units = "mm",
#                     dpi = 300,
#                     scale = 1,
#                     device = "pdf")

# ggsave(
#   filename = "S5_supp_fig1.pdf",
#   path = "figures/",
#   width = 85,
#   height = 85,
#   units = "mm",
#   dpi = 300,
#   scale = 1
# )
# 
# #For draft and co-author review
# ggsave(
#   filename = "S5_supp_fig1.png",
#   path = "figures/",
#   width = 85,
#   height = 85,
#   units = "mm",
#   dpi = 300,
#   scale = 1
# )
# 
# #convert pdf to TIFF for publisher
# pdftools::pdf_convert(
#   pdf = "figures/S5_supp_fig1.pdf",
#   format = "tiff",
#   filenames = "figures/S5_supp_fig1.tif",
#   dpi = 300
# )
