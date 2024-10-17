# Copyright (c) [2023] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we train NF sample classifiers
#' using random forests, weighted by umbalanced data
library(MASS)
library(MOFA2)
library(tidyverse)

# Main ---------------------------------------------------------------
setwd("/Users/ricardoramirez/Dropbox/PostDoc/Research/ReHeaT2/")
model_out <- "./results/LDA/model_stats.rds"
set.seed(1908)

input_df <- tibble(file = list.files("./data/pbulk/")) %>%
  dplyr::mutate(study = gsub("_pbulk[.]csv","",file),
                metadata_file =  paste0("./data/metadata_ext/", gsub("_pbulk[.]csv","_metadata.csv",file)),
                mofa_out = paste0("./data/models/", gsub("_pbulk[.]csv","_mofa.hdf5",file))) %>%
  dplyr::mutate(file = paste0("./data/pbulk/",file))

model_outs <- pmap(input_df, function(study, file,metadata_file, mofa_out) {

  print(study)

  model <- MOFA2::load_model(mofa_out)

  meta_data <- read_csv(metadata_file,
                        show_col_types = FALSE)

  LDA_dat <- MOFA2::get_factors(model)[[1]] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample_id") %>%
    left_join(meta_data[, c("sample_id", "heart_failure")], by = "sample_id") %>%
    dplyr::mutate(HF = ifelse(heart_failure == "NF", "no", "yes")) %>%
    dplyr::mutate(HF = factor(HF)) %>%
    dplyr::select(-sample_id) %>%
    dplyr::select(-heart_failure)

  LDA_model <- lda(HF ~ ., data = LDA_dat)

  return(LDA_model)

})

# Organize results

model_outs <- model_outs %>%
  enframe() %>%
  dplyr::mutate(name = input_df$study)

saveRDS(model_outs, model_out)




