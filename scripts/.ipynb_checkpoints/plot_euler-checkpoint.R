library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(tibble)
library(patchwork)
ancestry_colors <- c(
  "European"             = "#ffe119",
  "Eastern Bantu"        = "#f58231",
  "West Central African" = "#e6194B",
  "African-American"     = "#3cb44b",
  "Admixed American"     = "#3cb44b",
  "Indigenous American"  = "#c6dbef",
  "Afro Caribbean"       = "#3cb44c",
  "Afro South American"  = "#3cb44d"
  
)


label_map_afr_am <- c(
  "European" = "European",
  "Eastern_Bantu" = "Eastern Bantu",
  "West_Central_African" = "West Central African",
  #"African-American" = "African-American",
  "African-American" = "Admixed American",
  "Indigenous_American" = "Indigenous American"
  
)

label_map_afr_car <- c(
  "European" = "European",
  "Eastern_Bantu" = "Eastern Bantu",
  "West_Central_African" = "West Central African",
  #"Afro-Caribbean" = "Afro Caribbean",
  "Afro-Caribbean" = "Admixed American",
  "Indigenous_American" = "Indigenous American"
)


label_map_afro_south <- c(
  "European" = "European",
  "Eastern_Bantu" = "Eastern Bantu",
  "West_Central_African" = "West Central African",
  #"Afro-South_American" = "Afro South American",
  "Afro-South_American" = "Admixed American",
  "Indigenous_American" = "Indigenous American"
)

ellipse_to_polys <- function(euler_obj, label_map = NULL, npoints = 400) {
  ell <- as.data.frame(euler_obj$ellipses)
  ell$set <- rownames(ell)

  if (!is.null(label_map)) {
    ell$plot_label <- dplyr::recode(ell$set, !!!label_map)
  } else {
    ell$plot_label <- ell$set
  }

  ell$ellipse_id <- seq_len(nrow(ell))

  poly_list <- purrr::pmap(
    ell,
    function(h, k, a, b, phi, set, plot_label, ellipse_id, ...) {
      t <- seq(0, 2*pi, length.out = npoints)

      x0 <- a * cos(t)
      y0 <- b * sin(t)

      x <- h + x0 * cos(phi) - y0 * sin(phi)
      y <- k + x0 * sin(phi) + y0 * cos(phi)

      tibble::tibble(
        x = x,
        y = y,
        set = set,
        plot_label = plot_label,
        ellipse_id = ellipse_id
      )
    }
  )

  polys <- dplyr::bind_rows(poly_list)

  polys
}




# make_ancestry_plot <- function(euler_obj, label_map, panel_title, ancestry_colors) {
# 
#   # build polygon dataframe from this eulerr object
#   polys <- ellipse_to_polys(
#     euler_obj,
#     label_map = label_map,
#     npoints = 400
#   )
# 
#   # plot
#   ggplot(polys, aes(x = x, y = y)) +
#     geom_polygon(
#       aes(
#         group = ellipse_id,
#         fill  = plot_label
#       ),
#       alpha = 0.5,
#       color = "black",
#       linewidth = 0.4
#     ) +
#     scale_fill_manual(
#       values = ancestry_colors,
#       breaks = names(ancestry_colors),
#       name   = "Ancestry source"
#     ) +
#     coord_equal() +
#     theme_classic(base_size = 12) +
#     theme(
#       axis.line   = element_blank(),
#       axis.text   = element_blank(),
#       axis.ticks  = element_blank(),
#       axis.title  = element_blank(),
#       panel.border = element_rect(
#         color = "black",
#         fill = NA,
#         linewidth = 0.4
#       ),
#       legend.position = "right",
#       legend.title = element_text(size = 10, face = "bold", color = "black"),
#       legend.text  = element_text(size = 10, color = "black"),
#       plot.margin = margin(5,5,5,5)
#     ) +
#     labs(
#       title = panel_title
#     )
# }
# 
# 
# p_AA <- make_ancestry_plot(
#   euler_obj      = AfricanAmericanProxiesEuler$euler,
#   label_map      = label_map_afr_am,
#   panel_title    = "African American ancestry components",
#   ancestry_colors = ancestry_colors
# )
# p_AA
# 
# 
# p_Carib <- make_ancestry_plot(
#   euler_obj      = AfroCaribbeanProxiesEuler$euler,
#   label_map      = label_map_afr_car,
#   panel_title    = "Afro-Caribbean ancestry components",
#   ancestry_colors = ancestry_colors
# )
# 
# p_South <- make_ancestry_plot(
#   euler_obj      = AfroSouthAmericanProxies$euler,
#   label_map      = label_map_afro_south,
#   panel_title    = "Afro-South American ancestry components",
#   ancestry_colors = ancestry_colors
# )
# 
# combo <- p_AA + p_Carib + p_South + plot_layout(nrow = 1, guides = "collect") &
#   theme(legend.position = "right")
# combo
 ###start here#

#the code above works fine but its just awkward with legend stitching manually edited it below (breaks) to ignore the admixed group


make_ancestry_plot2 <- function(euler_obj, label_map, panel_title, ancestry_colors) {
  
  # build polygon dataframe from this eulerr object
  polys <- ellipse_to_polys(
    euler_obj,
    label_map = label_map,
    npoints = 400
  )
  
  # plot
  ggplot(polys, aes(x = x, y = y)) +
    geom_polygon(
      aes(
        group = ellipse_id,
        fill  = plot_label
      ),
      alpha = 0.5,
      color = "black",
      linewidth = 0.4
    ) +
    scale_fill_manual(
      values = ancestry_colors,
      breaks = c("European", "Eastern Bantu", "West Central African", "Indigenous American","Admixed American"),
      name   = "Ancestry source"
    ) +
    coord_equal() +
    theme_classic(base_size = 12) +
    theme(
      axis.line   = element_blank(),
      axis.text   = element_blank(),
      axis.ticks  = element_blank(),
      axis.title  = element_blank(),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.4
      ),
      legend.position = "right",
      legend.title = element_text(size = 10, face = "bold", color = "black"),
      legend.text  = element_text(size = 10, color = "black"),
      plot.margin = margin(5,5,5,5)
    ) +
    labs(
      title = panel_title
    )
}

p_AA2 <- make_ancestry_plot2(
  euler_obj      = AfricanAmericanProxiesEuler$euler,
  label_map      = label_map_afr_am,
  panel_title    = "African American ",
  ancestry_colors = ancestry_colors
)
p_AA2 
p_Carib2 <- make_ancestry_plot2(
  euler_obj      = AfroCaribbeanProxiesEuler$euler,
  label_map      = label_map_afr_car,
  panel_title    = "Afro-Caribbean",
  ancestry_colors = ancestry_colors
)
p_Carib2
p_South2 <- make_ancestry_plot2(
  euler_obj      = AfroSouthAmericanProxies$euler,
  label_map      = label_map_afro_south,
  panel_title    = "Afro-South American",
  ancestry_colors = ancestry_colors
)
p_South2
combo2 <- p_AA2 + p_Carib2 + p_South2 + plot_layout(nrow = 1, guides = "collect") &
  theme(legend.position = "right")
combo2

ggsave(
  "/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/geovar/results/admixed_americas.png",
  combo2,
  width  = 10,
  height = 4,
  dpi    = 600,
  bg     = "white"
)


ggsave(
  "all_ancestry_panels_oneLegend_clean.pdf",
  combo2,
  width  = 10,
  height = 4,
  bg     = "white"
)
