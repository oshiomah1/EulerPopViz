
library(readr)
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)

# read the table
trait_df_raw <- read_csv(
  "/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/dbNSFP_v4.9/data/14snps.csv"
)  %>%
  slice(c(1:4, 6:15))


trait_df <- trait_df_raw %>%
  mutate(
    # trim spaces in predictor calls
    BayesDel_noAF_pred  = str_trim(BayesDel_noAF_pred),
    BayesDel_addAF_pred = str_trim(BayesDel_addAF_pred),
    ESM1b_pred          = str_trim(ESM1b_pred),
    
    # turn "." or "" into NA for cleaner plotting
    BayesDel_noAF_pred  = na_if(BayesDel_noAF_pred,  "."),
    BayesDel_addAF_pred = na_if(BayesDel_addAF_pred, "."),
    ESM1b_pred          = na_if(ESM1b_pred,         "."),
    
    BayesDel_noAF_pred  = na_if(BayesDel_noAF_pred,  ""),
    BayesDel_addAF_pred = na_if(BayesDel_addAF_pred, ""),
    ESM1b_pred          = na_if(ESM1b_pred,          ""),
    
    # order GERP_35_scaled categories
    GERP_35_scaled = factor(
      GERP_35_scaled,
      levels = c("neutral", "large", "extreme")  # you can change this order if you want
    ),
    
    # make Trait a factor so order left→right is stable
    Trait = factor(Trait,
                   levels = c("Growth", "Hypertension", "Kidney disease",
                              "Malaria", "Pigmentation"))
  )

trait_long <- trait_df %>%
  pivot_longer(
    cols = c(BayesDel_noAF_pred, BayesDel_addAF_pred, ESM1b_pred),
    names_to = "Predictor",
    values_to = "Call"
  ) %>%
  mutate(
    Predictor = factor(
      Predictor,
      levels = c( "ESM1b_pred", "BayesDel_addAF_pred","BayesDel_noAF_pred" ),
      labels = c("ESM1b", "BayesDel (addAF)","BayesDel (noAF)")
    ),
    Call = factor(Call, levels = c("D", "T"))  # red = damaging (D), teal = tolerated (T)
  )

# p <- ggplot(
#   trait_long,
#   aes(
#     x = Trait,
#     y = GERP_35_scaled
#   )
# ) +
#   geom_point(
#     size = 3,
#     stroke = 0.4,
#     shape = 21,
#     aes(fill = Call),
#     color = "black"
#   ) +
#   facet_wrap(~ Predictor, nrow = 1) +
#   scale_fill_manual(
#     name = "Prediction",
#     values = c(
#       "D" = "#FFb683",  # damaging
#       "T" = "#77dd76"   # tolerated
#     ),
#     drop = FALSE
#   ) +
#   labs(
#     x = "Trait / phenotype association",
#     y = "GERP35 category",
#     title = "Trait-linked SNPs by predicted functional impact"
#   ) +
#   theme_classic(base_size = 12) +
#   theme(
#     axis.text.x = element_text(
#       angle = 45,
#       hjust = 1,
#       color = "black"
#     ),
#     axis.text.y = element_text(color = "black"),
#     axis.title  = element_text(color = "black"),
#     strip.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
#     strip.text = element_text(size = 11, face = "bold"),
#     legend.position = "right",
#     legend.title = element_text(size = 10, face = "bold"),
#     legend.text  = element_text(size = 10),
#     panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
#     plot.margin = margin(5,5,5,5)
#   )
# p


 

# 1. collapse to counts
trait_collapsed <- trait_long %>%
  count(Trait, GERP_35_scaled, Predictor, Call)

# 2. plot with dodge + size + tiny jitter
p <- ggplot(
  trait_collapsed,
  aes(
    x = Trait,
    y = GERP_35_scaled,
    fill = Call,
    size = n  # how many SNPs in that bin
  )
) +
  geom_point(
    stroke = 0.25,
    shape = 23,
    color = "black",
    position = position_jitterdodge(
      jitter.width = 0.95,   # tiny horizontal wiggle
      jitter.height = 0.0,  # tiny vertical wiggle so points don't sit perfectly on top
      dodge.width = 0.5      # separates D vs T
    )
  ) +
  facet_wrap(~ Predictor, nrow = 1) +
  scale_fill_manual(
    name = "Prediction",
    values = c(
      "D" = "#FFb683",
      "T" = "#77dd76"
    ),
    drop = FALSE
  ) +
  scale_size_continuous(
    name  = "SNP count",
    range = c(2, 6)  # min/max point size on plot
  ) +
  labs(
    x = "Trait / phenotype association",
    y = "GERP35 category",
    title = "Trait-linked SNPs by predicted functional impact"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(color = "black"),
    axis.title  = element_text(color = "black"),
    strip.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.4
    ),
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text  = element_text(size = 10),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.4
    ),
    plot.margin = margin(5,5,5,5)
  )

p





#PNG  
ggsave(
  filename = "trait_linked_SNPs_by_predictor.png",
  plot = p,
  width = 8,        # inches
  height = 3,       # inches
  dpi = 600,        # 600 is extra crisp
  bg = "white"      # make sure background isn't transparent
)
