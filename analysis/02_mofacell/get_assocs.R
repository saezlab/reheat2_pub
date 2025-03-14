# Copyright (c) [2024] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we test extra associations of the models

library(MOFAcellulaR)
library(tidyverse)
library(cowplot)
# Aesthetics - colors for plotting
source("./code/reheat2_pilot/aesthetics.R")

# Main ---------------------------------------------------------------
setwd("/Users/ricardoramirez/Dropbox/PostDoc/Research/ReHeaT2/")
model_stats_plt_file <- "./results/models/model_stats_covars.pdf"
model_stats_file <- "./results/models/model_stats_covars.csv"
hf_model_stats_plt_file <- "./results/models/model_stats_covars_hf.pdf"
hf_model_stats_file <- "./results/models/model_stats_covars_hf.csv"
# This is a function to estimate R2

estimate_covariate_imp <- function(model,
                                   metadata,
                                   sample_id_column,
                                   test_variable,
                                   test_type,
                                   categorical_type,
                                   group,
                                   pthresh = 0.05) {

  var_explained <- model@cache$variance_explained$r2_total[[1]] %>%
    enframe("view","R2")

  # First check that the covariate is present

  if(test_variable %in% colnames(metadata)) {

    ezy_meta <- metadata[,c(sample_id_column, test_variable)] %>%
      na.omit()

    # Test that the categorical value has at least 2 groups with 3
    # samples ------------------------------------------------

    if(test_type == "categorical") {

      useful_classes <- ezy_meta %>%
        group_by_at(test_variable) %>%
        summarise(n_samples = n()) %>%
        dplyr::filter(n_samples > 2) %>%
        pull(test_variable)

      # If there is more than 1 class then keep those classes
      if(length(useful_classes )> 1) {

        ezy_meta <- ezy_meta[ezy_meta[, test_variable][[1]] %in% useful_classes, ]

      } else {

        covar_R2 <- tibble(view = var_explained$view)
        covar_R2[[test_variable]] = NA

        return(covar_R2)

      }

    }

    if(nrow(ezy_meta) >= 3) {

      # Test association - this is one output
      expl_var_testvar <- MOFAcellulaR::get_associations(model = model,
                                                         metadata = ezy_meta,
                                                         sample_id_column = sample_id_column,
                                                         test_variable = test_variable,
                                                         test_type = test_type,
                                                         categorical_type = categorical_type,
                                                         group = group)

      # Get useful factors
      useful_factors <- expl_var_testvar %>%
        dplyr::filter(adj_pvalue <= pthresh) %>%
        pull(Factor)

      # Make the data frames based on association


      if(rlang::is_empty(useful_factors)) {
        covar_R2 <- tibble(view = var_explained$view)
        covar_R2[[test_variable]] = 0

        return(covar_R2)

      } else {
        covar_R2 <- model@cache$variance_explained$r2_per_factor[[1]][useful_factors,, drop = F] %>%
          colSums() %>%
          enframe("view", test_variable)

        return(covar_R2)

      }


    } else {

      covar_R2 <- tibble(view = var_explained$view)
      covar_R2[[test_variable]] = NA

      return(covar_R2)

    }


  } else {

    covar_R2 <- tibble(view = var_explained$view)
    covar_R2[[test_variable]] = NA

    return(covar_R2)

  }

}

# Prepare input data frames
input_df <- tibble(file = list.files("./data/pbulk/")) %>%
  dplyr::mutate(study = gsub("_pbulk[.]csv","",file),
                metadata_file =  paste0("./data/metadata_ext/", gsub("_pbulk[.]csv","_metadata.csv",file)),
                mofa_out = paste0("./data/models/", gsub("_pbulk[.]csv","_mofa.hdf5",file))) %>%
  dplyr::select(-file)


# Run associations with MOFAcell model loadings

model_outs <- pmap(input_df, function(study, metadata_file, mofa_out) {

  print(study)

  # Importing model
  model <- MOFA2::load_model(mofa_out)

  # Importing meta_data
  meta_data <- read_csv(metadata_file,
                        show_col_types = FALSE)

  # R2 - of model


  # R2 of covariates of interest ---------------------------------------

  cat_vars <- c("sex")

  cat_assocs <- map(cat_vars, function(x) {

    estimate_covariate_imp(model = model,
                           metadata = meta_data,
                           sample_id_column = "sample_id",
                           test_variable = x,
                           test_type = "categorical",
                           categorical_type = "parametric",
                           group = F,
                           pthresh = 0.05) %>%
      pivot_longer(-view, names_to = "covar", values_to = "R2")



  }) %>%
    enframe() %>%
    unnest() %>%
    dplyr::select(-name)

  cont_vars <- c("age", "LVEF", "BMI")

  cont_assocs <- map(cont_vars, function(x) {

    estimate_covariate_imp(model = model,
                           metadata = meta_data,
                           sample_id_column = "sample_id",
                           test_variable = x,
                           test_type = "continous",
                           group = F,
                           pthresh = 0.05) %>%
      pivot_longer(-view, names_to = "covar", values_to = "R2")



  }) %>%
    enframe() %>%
    unnest() %>%
    dplyr::select(-name)

  model_summ <- bind_rows(cat_assocs, cont_assocs) %>%
    dplyr::mutate(study = study)

  # R2 of covariates of interest ---------------------------------------

  meta_data_hf <- meta_data %>%
    dplyr::filter(heart_failure == "HF")

  cat_assocs_hf <- map(cat_vars, function(x) {

    estimate_covariate_imp(model = model,
                           metadata = meta_data_hf,
                           sample_id_column = "sample_id",
                           test_variable = x,
                           test_type = "categorical",
                           categorical_type = "parametric",
                           group = F,
                           pthresh = 0.05) %>%
      pivot_longer(-view, names_to = "covar", values_to = "R2")



  }) %>%
    enframe() %>%
    unnest() %>%
    dplyr::select(-name)

  cont_assocs_hf <- map(cont_vars, function(x) {

    estimate_covariate_imp(model = model,
                           metadata = meta_data_hf,
                           sample_id_column = "sample_id",
                           test_variable = x,
                           test_type = "continous",
                           group = F,
                           pthresh = 0.05) %>%
      pivot_longer(-view, names_to = "covar", values_to = "R2")



  }) %>%
    enframe() %>%
    unnest() %>%
    dplyr::select(-name)

  model_summ_hf <- bind_rows(cat_assocs_hf, cont_assocs_hf) %>%
    dplyr::mutate(study = study)

  return(list("all" = model_summ,
              "hf" = model_summ_hf))

})

# Generate final map
model_stats <- map(model_outs, ~ .x[["all"]] )

model_stats <- enframe(model_stats) %>%
  unnest(c(value)) %>%
  dplyr::select(-name)

write_csv(model_stats, model_stats_file)

model_stats_plt <- model_stats %>%
  ggplot(aes(x = study, y = R2)) +
  geom_boxplot() +
  facet_wrap(covar ~ ., nrow = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

pdf(model_stats_plt_file, height = 3.3, width = 2.7)

plot(model_stats_plt)

dev.off()


#
model_stats_hf <- map(model_outs, ~ .x[["hf"]] )

model_stats_hf <- enframe(model_stats_hf) %>%
  unnest(c(value)) %>%
  dplyr::select(-name)

write_csv(model_stats_hf, hf_model_stats_file)

model_stats_hf_plt <- model_stats_hf %>%
  ggplot(aes(x = study, y = R2)) +
  geom_boxplot() +
  facet_wrap(covar ~ ., nrow = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("R2 in HF samples")


pdf(hf_model_stats_plt_file, height = 3.3, width = 2.7)

plot(model_stats_hf_plt)

dev.off()

# Paper stats

model_stats <- read_csv(model_stats_file)


model_stats %>%
  dplyr::filter(study != "Reichart2022_DCM") %>%
  group_by(study) %>%
  summarize(mean_r2 = mean(R2, na.rm = T))


model_stats %>%
  dplyr::filter(study != "Reichart2022_DCM") %>%
  group_by(covar) %>%
  summarize(mean_r2 = mean(R2, na.rm = T))




