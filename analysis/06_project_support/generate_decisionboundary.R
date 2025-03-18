# Copyright (c) [2024] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Generate a linear discriminant from 2 factors

library(MOFAcellulaR)
library(tidyverse)
library(cowplot)
library(MASS)
source("./code/reheat2_pilot/aesthetics.R")

# Color definition of the projections based on Jan's colors
etiology_colors <- list(Sarcoidosis = "maroon3",
                        HF_HLHS = "orange") %>%
  unlist()

etiology_colors <- c(etiology_colors, col_pub$etiologies_colors)

# Main #########################################################################
metamofacell <- MOFA2::load_model("./results/meta_mofacell/metamodel_mofa.hdf5")
meta_data <- read_csv("./results/meta_mofacell/metamodel_meta.csv")[,c("study", "sample_id", "heart_failure")]

# Get all the tidy factors from the model
# fit a logistic regression model
# set a decision boundary

meta_factors <- MOFAcellulaR::get_tidy_factors(metamofacell,
                                               factor = "all",
                                               metadata = meta_data,
                                               sample_id_column = "sample_id",
                                               group = T)

meta_factors <- dplyr::filter(meta_factors,
              Factor %in% c("Factor1", "Factor2")) %>%
  dplyr::select(sample, heart_failure, Factor, value) %>%
  pivot_wider(names_from = Factor, values_from = value) %>%
  dplyr::mutate(heart_failure_bin = ifelse(heart_failure == "HF", 1, 0))

mdl <- glm(heart_failure_bin ~ .,
           data = meta_factors[,c("heart_failure_bin",
                                  "Factor1",
                                  "Factor2")],
           family=binomial)

slope <- coef(mdl)[2]/(-coef(mdl)[3])
intercept <- coef(mdl)[1]/(-coef(mdl)[3])

test <- stats::predict(mdl, meta_factors[,c("sample", "Factor1", "Factor2")],
                       type = "response")

original_decision_plt_grey <- ggplot(meta_factors, aes(x = Factor1, y = Factor2)) +
  geom_point(alpha = 0.3, color = "grey", size = 3) +
  geom_abline(slope = slope, intercept = intercept,
              color = "red", linetype = 2) +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank()) # Remove minor grid lines)

################################################################################

# Now get all projected data in the new space
# MI
MI_coords <- read_csv("./results/predict_project/MI/MI_factor_coords.csv")  %>%
  dplyr::select(sample_id, Factor1, Factor2, disease_code)

# LVAD
LVAD_coords <- read_csv("./results/predict_project/LVAD/LVAD_factor_coords.csv")  %>%
  dplyr::select(sample_id, Factor1, Factor2, disease_code) %>%
  dplyr::filter(disease_code != "NF")

# Pediatric - Fetal
Ped_coords <- read_csv("./results/predict_project/PED/ped_factor_coords.csv") %>%
  dplyr::rename("disease_code" = "group")  %>%
  dplyr::select(sample_id, Factor1, Factor2, disease_code)

# Pediatric
PedHill_coords <- read_csv("./results/predict_project/Hill/Hill_factor_coords.csv") %>%
  dplyr::rename("disease_code" = "Diagnosis") %>%
  dplyr::mutate(disease_code = ifelse(disease_code == "Donor",
                                      "NF",
                                      disease_code)) %>%
  dplyr::select(sample_id, Factor1, Factor2, disease_code)

# Sarcoidosis
Sarcoidosis_coords <- read_csv("./results/predict_project/HCA_reference/Sarcoidosis_factor_coords.csv") %>%
  dplyr::rename("disease_code" = "disease_state") %>%
  dplyr::select(sample_id, Factor1, Factor2, disease_code)

# HCM
HCM_coords <- read_csv("./results/predict_project/HCA_reference/HCM_Nicin_factor_coords.csv") %>%
  dplyr::mutate(disease_code = "HCM",
                sample_id = as.character(sample_id)) %>%
  dplyr::select(sample_id, Factor1, Factor2, disease_code)

# Building external dataset

external_wref <- bind_rows(MI_coords, LVAD_coords,
                           Ped_coords, PedHill_coords) %>%
  dplyr::mutate(type = "w/ ref")

external_woref <- bind_rows(Sarcoidosis_coords, HCM_coords) %>%
  dplyr::mutate(type = "wo/ ref")

external <- bind_rows(external_wref, external_woref)

external %>%
  write_csv("./Revision/figures/Figure6/Figure6B.csv")

# Plotting

lim_val <- max(abs(external$Factor2))

projection_plts <- original_decision_plt_grey +
  geom_point(data = external, aes(x = Factor1,
                                  y = Factor2,
                                  color = disease_code,
                                  shape = type),
             size = 3) +
  scale_color_manual(values = etiology_colors[unique(external$disease_code)]) +
  theme(legend.position = "right",
        axis.text = element_text(size =14),
        axis.title = element_text(size =14),
        legend.title = element_blank(),
        legend.text = element_text(size = 11),        # Adjust legend text size
        legend.key.size = unit(1.5, "lines")          # Adjust size of the legend keys
  ) +
  guides(
    color = guide_legend(nrow = 5) # Set the number of rows in the legend
  ) +
  xlim(lim_val * -1, lim_val) +
  ylim(lim_val * -1, lim_val) +
  coord_equal() +
  xlab("MCP1") +
  ylab("MCP2")

pdf("./results/predict_project/patient_map_proj.pdf", height = 6, width = 6)

plot(projection_plts)

dev.off()


ggplot_build(projection_plts)
# Misclassification test

external_noref <- bind_rows(Sarcoidosis_coords, HCM_coords)

test <- stats::predict(mdl, external_noref[,c("sample_id", "Factor1", "Factor2")],
                       type = "response")

external <- external_noref %>%
  dplyr::mutate(HF_true = test > 0.5,
                disease_code_resp = ifelse(disease_code == "NF", 0, 1)) %>%
  dplyr::mutate(HF_true = ifelse(HF_true == TRUE, 1, 0))

CER <- 1 - sum(external$HF_true == external$disease_code_resp)/nrow(external)


meta_factors


