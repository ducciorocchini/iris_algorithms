
# Spectral Scatterplot of Sentinel-2 Pixels

# Classifies a Sentinel-2B image with im.classify() (k-means)
# and plots the resulting clusters in Blue (B02) vs NIR (B08)
# space, with marginal density curves attached to the top and
# right edges.

# Data: Sentinel-2B, tile T30NXN, 26 January 2022
# Bands: B02 (Blue), B03 (Green), B04 (Red), B08 (NIR), all 10 m




library(terra)
library(ggplot2)
library(dplyr)
library(imageRy)
library(cowplot)



project_root <- "/Users/irisnanaobeng/Desktop/Scatter Plot Algorithim/Scatter_Plot_Algorithim"


blue  <- rast(file.path(project_root, "T30NXN_20220126T102209_B02_10m.jp2"))
green <- rast(file.path(project_root, "T30NXN_20220126T102209_B03_10m.jp2"))
red   <- rast(file.path(project_root, "T30NXN_20220126T102209_B04_10m.jp2"))
nir   <- rast(file.path(project_root, "T30NXN_20220126T102209_B08_10m.jp2"))



sentinel_stack <- c(blue, green, red, nir)
names(sentinel_stack) <- c("Blue", "Green", "Red", "NIR")



crop_extent <- ext(
  xmin(sentinel_stack) + 60000,
  xmin(sentinel_stack) + 64000,
  ymin(sentinel_stack) + 60000,
  ymin(sentinel_stack) + 64000
)

sentinel_crop <- crop(sentinel_stack, crop_extent)



sentinel_crop <- sentinel_crop / 10000


rgb_stack <- sentinel_crop[[c("Blue", "Green", "Red")]]

set.seed(42)

classified <- imageRy::im.classify(
  rgb_stack,
  num_clusters = 3,
  seed = 42,
  do_plot = FALSE
)


pixel_data <- data.frame(
  Blue    = as.numeric(values(sentinel_crop[["Blue"]], mat = FALSE)),
  NIR     = as.numeric(values(sentinel_crop[["NIR"]],  mat = FALSE)),
  Cluster = as.numeric(values(classified,              mat = FALSE))
)

pixel_data <- pixel_data |>
  filter(
    !is.na(Blue),
    !is.na(NIR),
    !is.na(Cluster)
  )


cluster_sizes  <- pixel_data |> count(Cluster, sort = TRUE)
cluster_labels <- cluster_sizes$Cluster

pixel_data <- pixel_data |>
  mutate(
    Cluster_Label = case_when(
      Cluster == cluster_labels[1] ~ "A",
      Cluster == cluster_labels[2] ~ "B",
      Cluster == cluster_labels[3] ~ "C"
    )
  )


set.seed(1)

plot_data <- pixel_data |>
  slice_sample(n = 300)


cluster_colours <- c(
  A = "#011959",
  B = "#D68A00",
  C = "#A3311D"
)

cluster_shapes <- c(
  A = 16,
  B = 17,
  C = 15
)


main_plot <- ggplot(plot_data, aes(x = Blue, y = NIR)) +
  
  geom_point(
    aes(
      colour = Cluster_Label,
      shape  = Cluster_Label
    ),
    size  = 1.8,
    alpha = 0.6
  ) +
  
  scale_colour_manual(values = cluster_colours) +
  scale_shape_manual(values  = cluster_shapes)  +
  
  labs(
    title    = "Spectral Scatterplot of Classified Sentinel-2 Pixels",
    subtitle = "300 sampled pixels with marginal cluster distributions",
    x        = "Blue reflectance (B02)",
    y        = "Near-infrared reflectance (B08)"
  ) +
  
  theme_minimal(base_size = 13) +
  
  theme(
    legend.position = "none",
    plot.title      = element_text(face = "bold"),
    plot.subtitle   = element_text(size = 11)
  )


blue_dens <- pixel_data |>
  group_by(Cluster_Label) |>
  group_modify(~ {
    d <- density(.x$Blue, na.rm = TRUE)
    data.frame(Blue = d$x, density = d$y)
  }) |>
  ungroup()

nir_dens <- pixel_data |>
  group_by(Cluster_Label) |>
  group_modify(~ {
    d <- density(.x$NIR, na.rm = TRUE)
    data.frame(NIR = d$x, density = d$y)
  }) |>
  ungroup()


blue_max <- blue_dens |>
  group_by(Cluster_Label) |>
  filter(density == max(density)) |>
  ungroup()

blue_labels <- blue_max |>
  mutate(
    label = Cluster_Label,
    
    hjust = case_when(
      Cluster_Label == "A" ~ 0,
      Cluster_Label == "B" ~ 0.5,
      Cluster_Label == "C" ~ 0
    ),
    
    vjust = case_when(
      Cluster_Label == "A" ~ 1,
      Cluster_Label == "B" ~ 0,
      Cluster_Label == "C" ~ 1
    ),
    
    nudge_x = case_when(
      Cluster_Label == "A" ~ 0.003,
      Cluster_Label == "B" ~ 0,
      Cluster_Label == "C" ~ 0.003
    ),
    
    nudge_y = case_when(
      Cluster_Label == "A" ~ 3,
      Cluster_Label == "B" ~ 5,
      Cluster_Label == "C" ~ 3
    )
  )


xdens <- axis_canvas(main_plot, axis = "x") +
  
  geom_ribbon(
    data = blue_dens,
    aes(
      x    = Blue,
      ymin = 0,
      ymax = density,
      fill = Cluster_Label,
      group = Cluster_Label
    ),
    alpha = 0.35
  ) +
  
  geom_line(
    data = blue_dens,
    aes(
      x      = Blue,
      y      = density,
      colour = Cluster_Label,
      group  = Cluster_Label
    ),
    linewidth = 0.5
  ) +
  
  geom_text(
    data = blue_labels,
    aes(
      x      = Blue + nudge_x,
      y      = density + nudge_y,
      label  = label,
      hjust  = hjust,
      vjust  = vjust,
      colour = Cluster_Label
    ),
    size       = 13 / .pt,
    fontface   = "bold",
    inherit.aes = FALSE
  ) +
  
  scale_colour_manual(values = cluster_colours, guide = "none") +
  scale_fill_manual(values   = cluster_colours, guide = "none") +
  
  scale_y_continuous(
    limits = c(0, max(blue_dens$density) * 1.10),
    expand = c(0, 0)
  ) +
  
  coord_cartesian(clip = "off")

ydens <- axis_canvas(main_plot, axis = "y", coord_flip = TRUE) +
  
  geom_ribbon(
    data = nir_dens,
    aes(
      x    = NIR,
      ymin = 0,
      ymax = density,
      fill = Cluster_Label,
      group = Cluster_Label
    ),
    alpha = 0.35
  ) +
  
  geom_line(
    data = nir_dens,
    aes(
      x      = NIR,
      y      = density,
      colour = Cluster_Label,
      group  = Cluster_Label
    ),
    linewidth = 0.5
  ) +
  
  scale_colour_manual(values = cluster_colours, guide = "none") +
  scale_fill_manual(values   = cluster_colours, guide = "none") +
  
  scale_y_continuous(
    limits = c(0, max(nir_dens$density) * 1.05),
    expand = c(0, 0)
  ) +
  
  coord_flip()

p1 <- insert_xaxis_grob(
  main_plot + theme(legend.position = "none"),
  xdens,
  grid::unit(3 * 14, "pt"),
  position = "top",
  clip = "off"
)

p2 <- insert_yaxis_grob(
  p1,
  ydens,
  grid::unit(3 * 14, "pt"),
  position = "right",
  clip = "off"
)

marginal_plot <- ggdraw(p2)

print(marginal_plot)


ggsave(
  file.path(project_root, "figures", "spectral_scatterplot_marginals.png"),
  plot   = marginal_plot,
  width  = 9,
  height = 7,
  dpi    = 300
)