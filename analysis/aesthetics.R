# Copyright (c) [2023] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' This script contains aesthetic utilities for general plots
#'

# Aesthetics
study_colors <- c("#68affc", "#afe642",
                  "#7446f3", "#e78361",
                  "#e78361", "#9c3b68",
                  "red")

names(study_colors) <- c("Chaffin2022_DCM", "Koenig2022_DCM",
                         "Kuppe2022_MI", "Reichart2022_DCM",
                         "LVReichart2022_DCM", "Simonson2023_ICM",
                         "Armute2023_LVAD")

col_pub <- readRDS("./data/color_list_figures.rds")

col_pub$etiologies_colors["DCM_rec"] <- "darkslategray3"

etiology_colors <- list(NF = "darkgrey",
                        ARVC = "maroon3",
                        NCC = "#ff69b4",
                        ICM = "gold",
                        HCM = "tan",
                        DCM = "darkgreen",
                        DCM_rec = "lightgreen",
                        MI_BZ = "#00ffff",
                        MI_RZ = "blue",
                        MI_IZ = "#6495ed",
                        MI_FZ = "#87CEFA",
                        FETAL = "#ff69b4") %>%
  unlist()

etiology_colors <- col_pub$etiologies_colors
