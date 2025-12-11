# List of required packages
packages <- c(
  "eulerr",
  "UpSetR",
  "jsonlite",
  "gridExtra",
  "cowplot",
  "magick",
  "grid"
)
# for (pkg in packages) {
#   if (!requireNamespace(pkg, quietly = TRUE)) {
#     install.packages(pkg)
#   }
# }

# Load packages
lapply(packages, library, character.only = TRUE)


# Define all population selections in a named list
pop_sets <- list(
  "African-American + Proxies" = c("European", "Eastern_Bantu", "West_Central_African", "African-American","Indigenous_American"),
  "Afro-South American + Proxies" = c("European", "Eastern_Bantu", "West_Central_African", "Afro-South_American","Indigenous_American"),
  "Afro-Caribbean + Proxies" = c("European", "Eastern_Bantu", "West_Central_African", "Afro-Caribbean","Indigenous_American"),
  "KS Components" = c("European", "Khoe-San", "Eastern_Bantu", "Asian"),
  "East Africa / Back to Africa" = c("European", "Eastern_Bantu", "Back-to-African", "Fulani", "Northeast_African"),
  "Major Continental Africa" = c("Khoe-San", "Eastern_Bantu", "West_Central_African", "Northeast_African", "Rainforest_Hunter-Gatherer")
)

# Define your label replacements ONCE
label_map <- c(
  "European" = "European",
  "Eastern_Bantu" = "Eastern Bantu",
  "West_Central_African" = "West Central African",
  "African-American" = "African-American",
  "Indigenous_American" = "Indigenous American",
  "Afro-South_American" = "Afro-South American",
  "Afro-Caribbean" = "Afro Caribbean",
  "Asian" = "Asian",
  "Khoe-San" = "Khoe-San",
  "Back-to-African" = "Back to African",
  "Fulani" = "Fulani",
  "Northeast_African" = "Northeast African",
  "Rainforest_Hunter-Gatherer" = "Rainforest Hunter-Gatherer"
)


sum_snps <- function(pop_combo, pops, counts, cutoff) {
  pop_tf <- rep(FALSE, length(pops))
  pop_tf[pop_combo] <- TRUE
  populations <- pops[pop_combo]
  if (length(pop_combo) > 1) {
    shared_common_snps <- which(rowSums(counts[,populations]>=cutoff)==length(populations))
  } else {
    shared_common_snps <- which(counts[,populations]>=cutoff)
  }
  snp_counts <- data.frame(t(c(pop_tf, sum(counts[shared_common_snps, ncol(counts)]))))
  colnames(snp_counts) <- c(pops, "shared_common_snps")
  return(snp_counts)
}

generate_euler_plot <- function(pairwise_file, poplist, selected_pops=c(), cutoff=3, common_pop="") {
  data <- read.table(
    file=pairwise_file,
    col.names=c("geovar_code", "counts"),
    colClasses=c("character", "numeric")
  )
  pops <- read.table(file=poplist)$V1
  if (length(selected_pops) < 1) {
    selected_pops <- pops
  }
  
  split.pops <- sapply(data[,1], function(a){strsplit(as.character(a),"")[[1]]})
  split.pops <- apply(split.pops,2,as.numeric)
  split.pops <- as.data.frame(t(split.pops))
  row.names(split.pops) <- 1:nrow(split.pops)
  colnames(split.pops) <- c(as.character(pops))
  data <- cbind(split.pops, data)
  data <- data[,c(selected_pops,"geovar_code","counts")]
  
  if (common_pop!=""){
    data <- data[which(data[common_pop]>=cutoff),]
  }
  
  combo_func <- Map(combn, list(1:length(selected_pops)), seq_along(1:length(selected_pops)), simplify=FALSE)
  combinations <- unlist(combo_func, recursive=FALSE)
  
  shared_common_snps <- do.call(rbind,lapply(combinations, FUN=sum_snps, pops=selected_pops, counts=data, cutoff=cutoff))
  shared_common_snps[nrow(shared_common_snps),"unique_snps"] <- shared_common_snps[nrow(shared_common_snps),"shared_common_snps"]
  # print("here are the common snps")
  #  print(shared_common_snps)
  
  for (i in (nrow(shared_common_snps)-1):1) {
    shared_common_snps[i,"unique_snps"] <- shared_common_snps[i,"shared_common_snps"] - sum(
      shared_common_snps[
        apply(
          shared_common_snps[,1:(ncol(shared_common_snps)-2)] - shared_common_snps[rep(i,nrow(shared_common_snps)),1:(ncol(shared_common_snps)-2)], 
          MARGIN=1, 
          FUN=function(x){!any(x<0)}
        ),"unique_snps"
      ], na.rm=TRUE
    )
  }
  
  euler_data <- shared_common_snps$unique_snps
  #print(shared_common_snps)
  #names(euler_data) <- apply(shared_common_snps, 1, function(x){paste(head(names(x),-2)[which(x==1)],collapse="&")}) #this is where code breaks, edit this so it applyies to n - ith column
  names(euler_data) <- apply(shared_common_snps[selected_pops], 1, function(x){paste(names(x)[which(x==1)],collapse="&")}) #new fix
  #print(names(euler_data))
  return(list("sets"=euler_data,"euler"=eulerr::euler(euler_data, shape = "ellipse"),"size"=sum(shared_common_snps$unique_snps)))
}

AfricanAmericanProxies = c("European", "Eastern_Bantu", "West_Central_African", "African-American","Indigenous_American")
AfroCaribbeanProxies = c("European", "Eastern_Bantu", "West_Central_African","Afro-Caribbean","Indigenous_American")
AfroSouthAmericanProxies = c("European", "Eastern_Bantu", "West_Central_African", "Afro-South_American","Indigenous_American")

KS.Components = c("European", "Khoe-San", "Eastern_Bantu", "Asian","West_Central_African")
East.Back.to.Africa = c("European", "Eastern_Bantu", "Back-to-African", "Fulani", "Northeast_African")
Major.Continental.Africa = c("Khoe-San", "Eastern_Bantu", "West_Central_African", "Northeast_African", "Rainforest_Hunter-Gatherer")
KS.Components2 = c("European", "Khoe-San", "Eastern_Bantu", "Asian")


AfricanAmericanProxiesEuler <- generate_euler_plot("/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/geovar/data/caapa_chr21.txt", 
                               "/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/geovar/data/pop_list_clean.txt", cutoff=3,selected_pops= AfricanAmericanProxies)
AfroCaribbeanProxiesEuler <- generate_euler_plot("/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/geovar/data/caapa_chr21.txt", 
                                                 "/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/geovar/data/pop_list_clean.txt", cutoff=3,selected_pops= AfroCaribbeanProxies)

AfroSouthAmericanProxies <- generate_euler_plot("/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/geovar/data/caapa_chr21.txt", 
                                                   "/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/geovar/data/pop_list_clean.txt", cutoff=3,selected_pops= AfroSouthAmericanProxies)
KS.Components.Euler =  generate_euler_plot("/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/geovar/data/caapa_chr21.txt", 
                                     "/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/geovar/data/pop_list_clean.txt", cutoff=3,selected_pops= KS.Components)
KS.Components.Euler2 =  generate_euler_plot("/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/geovar/results/allchr.txt", 
                                           "/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/geovar/data/pop_list_clean.txt", cutoff=3,selected_pops= KS.Components2)
East.Back.to.Africa.Euler =  generate_euler_plot("/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/geovar/data/caapa_chr21.txt", 
                                                  "/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/geovar/data/pop_list_clean.txt", cutoff=3,selected_pops= East.Back.to.Africa)
Major.Continental.Africa.Euler =  generate_euler_plot("/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/geovar/results/allchr.txt", 
                                                "/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/geovar/data/pop_list_clean.txt", cutoff=3,selected_pops= Major.Continental.Africa)


 
label_map_afr_am <- c(
  "European" = "European",
  "Eastern_Bantu" = "Eastern Bantu",
  "West_Central_African" = "West Central African",
  "African-American" = "African-American",
  "Indigenous_American" = "Indigenous American"
)

label_map_afr_car <- c(
  "European" = "European",
  "Eastern_Bantu" = "Eastern Bantu",
  "West_Central_African" = "West Central African",
  "Afro-Caribbean" = "Afro Caribbean",
  "Indigenous_American" = "Indigenous American"
)


label_map_afro_south <- c(
  "European" = "European",
  "Eastern_Bantu" = "Eastern Bantu",
  "West_Central_African" = "West Central African",
  "Afro-South_American" = "Afro South American",
  "Indigenous_American" = "Indigenous American"
)

KS.Components = c("European", "Khoe-San", "Eastern_Bantu", "Asian","West_Central_African")
East.Back.to.Africa = c("European", "Eastern_Bantu", "Back-to-African", "Fulani", "Northeast_African")
Major.Continental.Africa = c("Khoe-San", "Eastern_Bantu", "West_Central_African", "Northeast_African", "Rainforest_Hunter-Gatherer")


label_map_KS <- c(
  "European" = "European",
  "Khoe-San" = "Khoe San",
  "Eastern_Bantu" = "Eastern Bantu",
  "Asian" = "Asian",
  "West_Central_African" = "West Central African"
)



label_map_East.Back.to.Africa <- c(
  "European" = "European",
  "Eastern_Bantu" = "Eastern Bantu",
  "Back-to-African" = "Back to African",
  "Fulani" = "Fulani",
  "Northeast_African" = "Northeast African"
)




# 
plot(
  KS.Components.Euler$euler,
  labels = label_map_KS,
  fills = list(alpha = 0.6),
  edges = TRUE
)
plot(
  KS.Components.Euler$euler,
  fills = list(alpha = 0.6),
  edges = TRUE
)


plot(
  East.Back.to.Africa.Euler$euler,
  labels=label_map_East.Back.to.Africa, 
  fills = list(alpha = 0.6),
  edges = TRUE
)
plot(
  East.Back.to.Africa.Euler$euler,
  fills = list(alpha = 0.6),
  edges = TRUE
)
plot(
  Major.Continental.Africa.Euler$euler,
  labels=label_map_Major.Continental.Africa,
  fills = list(alpha = 0.6),
  edges = TRUE
)
 
# 
# 
# UpSetR::upset(
#   fromExpression(AfricanAmericanProxiesEuler$sets),
#   order.by = "freq",
#   mainbar.y.label = "Number of Common Variants",
#   sets.x.label = "Set Size",
#   #main = paste("UpSet Plot -", name),
#   text.scale = c(2, 1.5, 1.5, 1.2, 1.5, 2)
# )
# 
#  

##############

library(dplyr)
library(purrr)
library(stringr)
library(tibble)

sets_to_binary_df <- function(set_vector) {
  # set_vector is like AfricanAmericanProxiesEuler$sets
  # names(set_vector) are combos separated by "&"
  combos <- strsplit(names(set_vector), "&")
  
  # all unique population labels across all combinations
  all_groups <- sort(unique(unlist(combos)))
  
  # build a block of rows for each combo
  df_list <- map2(combos, as.numeric(set_vector), function(groups, n_rows) {
    # make n_rows rows initialized to FALSE for all groups
    m <- matrix(FALSE, nrow = n_rows, ncol = length(all_groups),
                dimnames = list(NULL, all_groups))
    # mark the groups in this intersection TRUE
    m[, groups] <- TRUE
    as_tibble(m)
  })
  
  bind_rows(df_list)
}


aa_upset_df <- sets_to_binary_df(AfricanAmericanProxiesEuler$sets)
carib_upset_df <- sets_to_binary_df(AfroCaribbeanProxiesEuler$sets)
south_upset_df <- sets_to_binary_df(AfroSouthAmericanProxies$sets)
Major.Continental.Africa_df <- sets_to_binary_df(Major.Continental.Africa.Euler$sets)
KS2_df <- sets_to_binary_df(KS.Components.Euler2$sets)
library(ComplexUpset)
library(ggplot2)

label_map_afr_am_upset <- c(
  "African-American" = "African-American",
  "Eastern Bantu" = "Eastern_Bantu",
  "European" = "European",
  "Indigenous American" = "Indigenous_American",
  "West Central African" = "West_Central_African"
)

label_map_afr_car_upset <-  c(
  "Afro-Caribbean" = "Afro-Caribbean",
  "Eastern Bantu" = "Eastern_Bantu",
  "European" = "European",
  "Indigenous American" = "Indigenous_American",
  "West Central African" = "West_Central_African"
)


label_map_afro_south_upset <-c(
  "Eastern Bantu" = "Eastern_Bantu",
  "Afro-South American" = "Afro-South_American",
  "European" = "European",
  "Indigenous American" = "Indigenous_American",
  "West Central African" = "West_Central_African"
)

label_map_KS2 <- c(
  "European" = "European",
  "Khoe San" = "Khoe-San",
  "Eastern Bantu" = "Eastern_Bantu",
  "Asian" = "Asian"
)


label_map_Major.Continental.Africa<- c(
  "Khoe San" = "Khoe-San",
  "Eastern Bantu" = "Eastern_Bantu",
  "West Central African" = "West_Central_African",
  "Northeast African" = "Northeast_African",
  "Rainforest Hunter Gatherer" = "Rainforest_Hunter-Gatherer"
)


aa_upset_df <- sets_to_binary_df(AfricanAmericanProxiesEuler$sets)
carib_upset_df <- sets_to_binary_df(AfroCaribbeanProxiesEuler$sets)
south_upset_df <- sets_to_binary_df(AfroSouthAmericanProxies$sets)
 

aa_upset_df_pretty <- aa_upset_df %>%
  dplyr::rename(!!!label_map_afr_am_upset)

aa_sets_pretty <- colnames(aa_upset_df_pretty)
 
# choose the set columns (these must match your column names in aa_upset_df)
aa_sets <- colnames(aa_upset_df_pretty)

p_upset_aa <- upset(
  aa_upset_df_pretty,
  intersect   = aa_sets_pretty,
  min_size    = 1,
  width_ratio = 0.2,
  
  base_annotations = list(
    `Number of\nCommon Variants` = intersection_size(
      counts = FALSE    # <- don't print numbers on top of bars
    )
  ),
  
  set_sizes = FALSE
) +
  theme_classic(base_size = 11) +
  theme(
    axis.text  = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    plot.margin  = margin(5,5,5,5),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(
      size = 11,
      face = "bold",
      hjust = 0.5,
      color = "black"
    )
  ) +
  labs(title = "African-American")


p_upset_aa

#CARRIBEAN GGPLOT 

carib_upset_df <- sets_to_binary_df(AfroCaribbeanProxiesEuler$sets)


carib_upset_df_pretty <- carib_upset_df %>%
  dplyr::rename(!!!label_map_afr_car_upset)

carib_sets_pretty <- colnames(carib_upset_df_pretty)

# choose the set columns (these must match your column names in aa_upset_df)
carib_sets <- colnames(carib_upset_df_pretty)


p_upset_car <- upset(
  carib_upset_df_pretty,
  intersect   = carib_sets_pretty,
  min_size    = 1,
  width_ratio = 0.2,
  
  base_annotations = list(
    `Number of\nCommon Variants` = intersection_size(
      counts = FALSE    # <- don't print numbers on top of bars
    )
  ),
  
  set_sizes = FALSE
) +
  theme_classic(base_size = 11) +
  theme(
    axis.text  = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    plot.margin  = margin(5,5,5,5),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(
      size = 11,
      face = "bold",
      hjust = 0.5,
      color = "black"
    )
  ) +
  labs(title = "Afro-Carribean")
p_upset_car 



#AFRO SOUTH AM UPSET PLOT
south_upset_df <- sets_to_binary_df(AfroSouthAmericanProxies$sets)

south_upset_df_pretty <- south_upset_df %>%
  dplyr::rename(!!!label_map_afro_south_upset)

south_sets_pretty <- colnames(south_upset_df_pretty)

# choose the set columns (these must match your column names in aa_upset_df)
south_sets <- colnames(south_upset_df_pretty)


p_upset_south <-upset(
  south_upset_df_pretty,
  intersect   = south_sets_pretty,
  min_size    = 1,
  width_ratio = 0.2,
  
  base_annotations = list(
    `Number of\nCommon Variants` = intersection_size(
      counts = FALSE    # <- don't print numbers on top of bars
    )
  ),
  
  set_sizes = FALSE
) +
  theme_classic(base_size = 11) +
  theme(
    axis.text  = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    plot.margin  = margin(5,5,5,5),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(
      size = 11,
      face = "bold",
      hjust = 0.5,
      color = "black"
    )
  ) +
  labs(title = "Afro-South American")
p_upset_south


#KS componest no mwest afr PLOT
ks2_upset_df <- sets_to_binary_df(KS.Components.Euler2$sets)

ks2_upset_df_pretty <- ks2_upset_df %>%
  dplyr::rename(!!!label_map_KS2)

ks2_sets_pretty <- colnames(ks2_upset_df_pretty)


# choose the set columnsks2 (these must match your column names in aa_upset_df)
ks2_sets <- colnames(ks2_upset_df_pretty)


p_upset_KS2 <-upset(
  ks2_upset_df_pretty,
  intersect   = ks2_sets_pretty,
  min_size    = 1,
  width_ratio = 0.2,
  
  base_annotations = list(
    `Number of\nCommon Variants` = intersection_size(
      counts = FALSE    # <- don't print numbers on top of bars
    )
  ),
  
  set_sizes = FALSE
) +
  theme_classic(base_size = 11) +
  theme(
    axis.text  = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    plot.margin  = margin(5,5,5,5),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(
      size = 11,
      face = "bold",
      hjust = 0.5,
      color = "black"
    )
  ) +
  labs(title = "Khoe San Ancestry Components")
p_upset_KS2


###

#KS componest no mwest afr PLOT
major_afr_df <- sets_to_binary_df(Major.Continental.Africa.Euler$sets)

major_afr_upset_df_pretty <-Major.Continental.Africa_df %>%
  dplyr::rename(!!!label_map_Major.Continental.Africa)

 
major_afr_sets_pretty <- colnames(major_afr_upset_df_pretty)


# choose the set columnsks2 (these must match your column names in aa_upset_df)
major_afr_sets <- colnames(major_afr_upset_df_pretty)


p_upset_major_afr <-upset(
  major_afr_upset_df_pretty,
  intersect   = major_afr_sets_pretty,
  min_size    = 1,
  width_ratio = 0.2,
  
  base_annotations = list(
    `Number of\nCommon Variants` = intersection_size(
      counts = FALSE    # <- don't print numbers on top of bars
    )
  ),
  
  set_sizes = FALSE
) +
  theme_classic(base_size = 11) +
  theme(
    axis.text  = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    plot.margin  = margin(5,5,5,5),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(
      size = 11,
      face = "bold",
      hjust = 0.5,
      color = "black"
    )
  ) +
  labs(title = "Major Continental Ancestry Components")
p_upset_major_afr

ggsave(
  "major_cont_africa2.png",
  p_upset_major_afr,
  width  = 12,
  height = 10,
  dpi    = 600,
  bg     = "white"
)



ggsave(
  "KS_no_west_africa2.png",
  p_upset_KS2,
  width  = 12,
  height = 10,
  dpi    = 600,
  bg     = "white"
)

###


library(patchwork)



p_upset_aa_wrapped    <- wrap_elements(full = p_upset_aa)
p_upset_car_wrapped   <- wrap_elements(full = p_upset_car)
p_upset_south_wrapped <- wrap_elements(full = p_upset_south)


combo_upset <- (
  p_upset_aa_wrapped |
    p_upset_car_wrapped |
    p_upset_south_wrapped
) +
  plot_layout(
    widths = c(1, 1, 1)  # you can tweak this if one needs more horizontal space
  )

combo_upset

ggsave(
  filename = "upset_all_groups.pdf",
  plot     = combo_upset,
  width    = 12,
  height   = 4,
  bg       = "white"
)

combo_upset_vertical <- (
  p_upset_aa_wrapped /
    p_upset_car_wrapped /
    p_upset_south_wrapped
) +
  plot_layout(
    heights = c(1, 1, 1)  # you can tweak this if one panel needs more vertical real estate
  ) &
  theme(
    plot.margin = margin(5,10,5,10)  # a little breathing room around every panel
  )

# 3. Preview
combo_upset_vertical


ggsave(
  filename = "upset_all_groups2.pdf",
  plot     = combo_upset_vertical,
  width    = 12,
  height   = 12,
  bg       = "white"
)

ggsave(
  "upset_all_groups_vertical.png",
  combo_upset_vertical,
  width  = 12,
  height = 12,
  dpi    = 600,
  bg     = "white"
)
