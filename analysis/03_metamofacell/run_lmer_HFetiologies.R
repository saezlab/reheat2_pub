# Copyright (c) [2025] [Ricardo O. Ramirez Flores]
# flores@ebi.ac.uk

#' Runs specific linear mixed model with disease code
library(tidyverse)
library(lmerTest)

source("./code/reheat2_pilot/aesthetics.R")

# For association with HF we don't need these disease codes from Kuppe et. al.
input_df <- tibble(file = list.files("./data/pbulk/")) %>%
  dplyr::mutate(study = gsub("_pbulk[.]csv","",file),
                coldata_file =  paste0("./data/coldata/", gsub("_pbulk[.]csv","_coldata.csv",file)),
                metadata_file =  paste0("./data/metadata_ext/", gsub("_pbulk[.]csv","_metadata.csv",file))) %>%
  dplyr::filter(study != "Reichart2022_DCM") %>%
  dplyr::mutate(file = paste0("./data/pbulk/",file)) %>%
  dplyr::select(file, study, coldata_file, metadata_file)

# Generating the meta-data for all
all_meta <-  input_df %>%
  dplyr::mutate(metadata = map(metadata_file, read_csv)) %>%
  dplyr::select(study, metadata) %>%
  unnest()

# Extracting the info from the model
metamofacell <- MOFA2::load_model("./results/meta_mofacell/metamodel_mofa.hdf5")

assocs_df <- map(set_names(c("Factor1", "Factor2")), function(Fact) {

  factor_info <- MOFAcellulaR::get_tidy_factors(metamofacell,
                                                all_meta,
                                                sample_id_column = "sample_id",
                                                factor = Fact,
                                                group = T)

  factor_info <- factor_info %>%
    dplyr::filter(heart_failure == "HF") %>%
    dplyr::filter(disease_code != "ARVC")

  factor_model <- lmerTest::lmer("value ~ 1 + disease_code +  (1 | study)",
                                 data = factor_info)

  return(broom::tidy(anova(factor_model)))

})


NF_assocs_df <- map(set_names(c("Factor1", "Factor2")), function(Fact) {

  factor_info <- MOFAcellulaR::get_tidy_factors(metamofacell,
                                                all_meta,
                                                sample_id_column = "sample_id",
                                                factor = Fact,
                                                group = T)

  factor_info <- factor_info %>%
    dplyr::filter(heart_failure == "NF")

  factor_model <- lm("value ~ 1 + study",
                                 data = factor_info)

  return(broom::tidy(anova(factor_model)))

})

factor_info <- MOFAcellulaR::get_tidy_factors(metamofacell,
                                              all_meta,
                                              sample_id_column = "sample_id",
                                              factor = "all",
                                              group = T)

HF_dist_disease <- factor_info %>%
  dplyr::filter(Factor %in% c("Factor1", "Factor2"),
                heart_failure == "HF") %>%
  ggplot(aes(x =  study,fill = disease_code, y = value)) +
  geom_boxplot() +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust =0.5)) +
  facet_wrap(.~Factor, scales = "free_y") +
  scale_fill_manual(values = etiology_colors[factor_info$disease_code %>% unique()])

NF_dist_disease <- factor_info %>%
  dplyr::filter(Factor %in% c("Factor1", "Factor2"),
                heart_failure == "NF") %>%
  ggplot(aes(x =  study, y = value, fill = disease_code)) +
  geom_boxplot() +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust =0.5)) +
  facet_wrap(.~Factor, scales = "free_y") +
  scale_fill_manual(values = etiology_colors[factor_info$disease_code %>% unique()])


cowplot::plot_grid(HF_dist_disease, NF_dist_disease, nrow = 2)

