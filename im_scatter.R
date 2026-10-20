# im_scatter_run.R
# Scatterplot with marginal densities.
# Works with any number of clusters.
# Cluster colours match imageRy::im.classify() by default.


library(imageRy)
library(terra)
library(ggplot2)
library(dplyr)
library(cowplot)


im.scatter <- function(
    band_x,
    band_y,
    classified,
    n_sample    = 300,
    seed_sample = 1,
    colours     = NULL,
    shapes      = NULL,
    labels      = NULL,
    title       = "Scatterplot of Classified NDVI Pixels",
    subtitle    = "300 sampled pixels with marginal cluster distributions",
    xlab        = "NDVI (winter)",
    ylab        = "NDVI (summer)",
    marginal_pt = 3 * 14
) {
  
  if (!terra::compareGeom(band_x, band_y,     stopOnError = FALSE) ||
      !terra::compareGeom(band_x, classified, stopOnError = FALSE)) {
    stop("band_x, band_y and classified must share the same extent, ",
         "resolution and number of pixels.")
  }
  
  pixel_data <- data.frame(
    X       = as.numeric(terra::values(band_x,     mat = FALSE)),
    Y       = as.numeric(terra::values(band_y,     mat = FALSE)),
    Cluster = as.numeric(terra::values(classified, mat = FALSE))
  )
  pixel_data <- pixel_data |>
    dplyr::filter(!is.na(X), !is.na(Y), !is.na(Cluster))
  
  if (nrow(pixel_data) == 0) stop("No valid pixels found.")
  
  # Number of clusters 
  cluster_sizes <- pixel_data |>
    dplyr::count(Cluster, sort = TRUE)
  
  n_clusters <- nrow(cluster_sizes)
  
  #  Labels: A, B, C, ... as needed 
  if (is.null(labels)) {
    labels <- LETTERS[seq_len(n_clusters)]
  }
  if (length(labels) != n_clusters) {
    stop(sprintf("Found %d clusters but %d labels supplied.",
                 n_clusters, length(labels)))
  }
  
  pixel_data <- pixel_data |>
    dplyr::mutate(
      Cluster_Label = labels[match(Cluster, cluster_sizes$Cluster)]
    )
  
  # Colours: match imageRy::im.classify() default palette 
  if (is.null(colours)) {
   
    auto_colours <- scales::hue_pal()(n_clusters)
    colours <- stats::setNames(auto_colours, labels)
  }
  
  # Shapes: cycle to fit n_clusters 
  if (is.null(shapes)) {
    base_shapes <- c(16, 17, 15, 18, 8, 3, 4, 1, 2, 0)
    shapes <- stats::setNames(
      rep(base_shapes, length.out = n_clusters),
      labels
    )
  }
  
  set.seed(seed_sample)
  plot_data <- if (nrow(pixel_data) > n_sample) {
    dplyr::slice_sample(pixel_data, n = n_sample)
  } else {
    pixel_data
  }
  
  main_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = X, y = Y)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(colour = Cluster_Label, shape = Cluster_Label),
      size = 1.8, alpha = 0.6
    ) +
    ggplot2::scale_colour_manual(values = colours) +
    ggplot2::scale_shape_manual(values  = shapes) +
    ggplot2::labs(
      title    = title,
      subtitle = subtitle,
      x        = xlab,
      y        = ylab
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      legend.position  = "none",
      plot.title       = ggplot2::element_text(face = "bold"),
      plot.subtitle    = ggplot2::element_text(size = 11),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        linewidth = 0.3, colour = "grey88"
      )
    )
  
  x_dens <- pixel_data |>
    dplyr::group_by(Cluster_Label) |>
    dplyr::group_modify(~ {
      d <- stats::density(.x$X, na.rm = TRUE)
      data.frame(value = d$x, density = d$y)
    }) |>
    dplyr::ungroup()
  
  y_dens <- pixel_data |>
    dplyr::group_by(Cluster_Label) |>
    dplyr::group_modify(~ {
      d <- stats::density(.x$Y, na.rm = TRUE)
      data.frame(value = d$x, density = d$y)
    }) |>
    dplyr::ungroup()
  
  # --- Peak positions for labels --------------------------------------
  x_peaks <- x_dens |>
    dplyr::group_by(Cluster_Label) |>
    dplyr::filter(density == max(density)) |>
    dplyr::ungroup()
  
  y_span <- diff(range(x_dens$density, na.rm = TRUE))
  
  x_peaks <- x_peaks |>
    dplyr::mutate(
      label   = Cluster_Label,
      hjust   = 0.5,
      vjust   = 0,
      nudge_x = 0,
      nudge_y = y_span * 0.05
    )
  
  xdens <- cowplot::axis_canvas(main_plot, axis = "x") +
    ggplot2::geom_ribbon(
      data = x_dens,
      ggplot2::aes(x = value, ymin = 0, ymax = density,
                   fill = Cluster_Label, group = Cluster_Label),
      alpha = 0.35
    ) +
    ggplot2::geom_line(
      data = x_dens,
      ggplot2::aes(x = value, y = density,
                   colour = Cluster_Label, group = Cluster_Label),
      linewidth = 0.5
    ) +
    ggplot2::geom_text(
      data = x_peaks,
      ggplot2::aes(x = value + nudge_x, y = density + nudge_y,
                   label = label, hjust = hjust, vjust = vjust,
                   colour = Cluster_Label),
      size = 11 / ggplot2::.pt, fontface = "bold", inherit.aes = FALSE
    ) +
    ggplot2::scale_colour_manual(values = colours, guide = "none") +
    ggplot2::scale_fill_manual(  values = colours, guide = "none") +
    ggplot2::scale_y_continuous(
      limits = c(0, max(x_dens$density) * 1.25),
      expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(clip = "off")
  
  ydens <- cowplot::axis_canvas(main_plot, axis = "y", coord_flip = TRUE) +
    ggplot2::geom_ribbon(
      data = y_dens,
      ggplot2::aes(x = value, ymin = 0, ymax = density,
                   fill = Cluster_Label, group = Cluster_Label),
      alpha = 0.35
    ) +
    ggplot2::geom_line(
      data = y_dens,
      ggplot2::aes(x = value, y = density,
                   colour = Cluster_Label, group = Cluster_Label),
      linewidth = 0.5
    ) +
    ggplot2::scale_colour_manual(values = colours, guide = "none") +
    ggplot2::scale_fill_manual(  values = colours, guide = "none") +
    ggplot2::scale_y_continuous(
      limits = c(0, max(y_dens$density) * 1.05),
      expand = c(0, 0)
    ) +
    ggplot2::coord_flip()
  
  p1 <- cowplot::insert_xaxis_grob(
    main_plot + ggplot2::theme(legend.position = "none"),
    xdens, grid::unit(marginal_pt, "pt"),
    position = "top", clip = "off"
  )
  p2 <- cowplot::insert_yaxis_grob(
    p1, ydens, grid::unit(marginal_pt, "pt"),
    position = "right", clip = "off"
  )
  
  cowplot::ggdraw(p2)
}


# Load NDVI data 
ndvi_winter <- im.import("Sentinel2_NDVI_2020-02-21.tif")
ndvi_spring <- im.import("Sentinel2_NDVI_2020-05-21.tif")
ndvi_summer <- im.import("Sentinel2_NDVI_2020-08-01.tif")
ndvi_autumn <- im.import("Sentinel2_NDVI_2020-11-27.tif")

ndvi_stack <- c(ndvi_winter, ndvi_spring, ndvi_summer, ndvi_autumn)
names(ndvi_stack) <- c("Winter", "Spring", "Summer", "Autumn")


# Cluster on three bands 
ndvi_for_clustering <- c(ndvi_winter, ndvi_spring, ndvi_summer)
names(ndvi_for_clustering) <- c("Winter", "Spring", "Summer")

set.seed(42)
classified <- im.classify(
  ndvi_for_clustering,
  num_clusters = 3,
  seed = 42,
  do_plot = FALSE
)


# Build the figure 
figure <- im.scatter(
  band_x     = ndvi_stack[["Winter"]],
  band_y     = ndvi_stack[["Summer"]],
  classified = classified
)


ggsave(
  "ndvi_scatterplot_marginals.png",
  plot   = figure,
  width  = 9,
  height = 7,
  dpi    = 300
)


print(figure)


# Demonstration: same pipeline with 6 clusters 
classified_6 <- im.classify(
  ndvi_for_clustering,
  num_clusters = 6,
  seed = 42,
  do_plot = FALSE
)

figure_6 <- im.scatter(
  band_x     = ndvi_stack[["Winter"]],
  band_y     = ndvi_stack[["Summer"]],
  classified = classified_6
)

ggsave(
  "ndvi_scatterplot_marginals_6clusters.png",
  plot   = figure_6,
  width  = 9,
  height = 7,
  dpi    = 300
)

print(figure_6)