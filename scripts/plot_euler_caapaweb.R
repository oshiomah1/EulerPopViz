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
  "Fulani" = "black",
  "Indigenous American"  = "#c6dbef",
  "Rainforest Hunter Gatherer" = "brown",
  "Northeast African" = "blue",
  "Asian" = "pink",
  "Khoe San" ="#996"
)


ancestry_colors <- c(
  "European"             = "#ffe119",  # bright yellow
  "Eastern Bantu"        = "#f58231",  # vivid orange
  "West Central African" = "#e6194B",  # red
  "Fulani"               = "#800080",  # deep purple
  "Indigenous American"  = "#56B4E9",  # sky blue
  "Rainforest Hunter Gatherer" = "#654321",  # rich dark brown
  "Northeast African"    = "#0072B2",  # deep blue
  "Asian"                = "#FF69B4",  # bright pink
  "Khoe San"             = "#3cb44b" , # dark goldenrod
  "Back to African"  ="#180501"
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




make_ancestry_plot <- function(euler_obj, label_map, panel_title, ancestry_colors) {

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
      breaks = names(ancestry_colors),
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

label_map_KS <- c(
  "European" = "European",
  "Khoe-San" = "Khoe San",
  "Eastern_Bantu" = "Eastern Bantu",
  "Asian" = "Asian",
  "West_Central_African" = "West Central African"
)


label_map_KS2 <- c(
  "European" = "European",
  "Khoe-San" = "Khoe San",
  "Eastern_Bantu" = "Eastern Bantu",
  "Asian" = "Asian"
)
label_map_East.Back.to.Africa <- c(
  "European" = "European",
  "Eastern_Bantu" = "Eastern Bantu",
  "Back-to-African" = "Back to African",
  "Fulani" = "Fulani",
  "Northeast_African" = "Northeast African"
)

label_map_Major.Continental.Africa<- c(
  "Khoe-San" = "Khoe San",
  "Eastern_Bantu" = "Eastern Bantu",
  "West_Central_African" = "West Central African",
  "Northeast_African" = "Northeast African",
  "Rainforest_Hunter-Gatherer" = "Rainforest Hunter Gatherer"
)
p_KS <- make_ancestry_plot(
  euler_obj      =KS.Components.Euler$euler,
  label_map      = label_map_KS,
  panel_title    = "Khoe San ancestry components",
  ancestry_colors = ancestry_colors
) 
p_KS

p_KS2 <- make_ancestry_plot(
  euler_obj      =KS.Components.Euler2$euler,
  label_map      = label_map_KS2,
  panel_title    = "Khoe San ancestry components",
  ancestry_colors = ancestry_colors
) 
p_KS2
 

p_label_map_Major.Continental.Africa <- make_ancestry_plot(
  euler_obj      = Major.Continental.Africa.Euler$euler,
  label_map      = label_map_Major.Continental.Africa,
  panel_title    = "Major Continental Africa",
  ancestry_colors = ancestry_colors
)
p_label_map_Major.Continental.Africa

p_label_map_East.Back.to.Africa<- make_ancestry_plot(
  euler_obj      = East.Back.to.Africa.Euler$euler,
  label_map      = label_map_East.Back.to.Africa,
  panel_title    = "Eastern & Back to Africa",
  ancestry_colors = ancestry_colors
)
p_label_map_East.Back.to.Africa

ggsave(
  "/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/geovar/results/KS.png",
  p_KS,
  width  = 10,
  height = 4,
  dpi    = 600,
  bg     = "white"
)
ggsave(
  "/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/geovar/results/KS_all_chr_nowestafr.png",
  p_KS2,
  width  = 10,
  height = 4,
  dpi    = 600,
  bg     = "white"
)
ggsave(
  "/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/geovar/results/MajorContinentalAfricaall_chr.png",
  p_label_map_Major.Continental.Africa,
  width  = 10,
  height = 4,
  dpi    = 600,
  bg     = "white"
)

ggsave(
  "/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/geovar/results/EastBacktoAfrica.png",
  p_label_map_East.Back.to.Africa,
  width  = 10,
  height = 4,
  dpi    = 600,
  bg     = "white"
)
 