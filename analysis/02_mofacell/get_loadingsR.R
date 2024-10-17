# Copyright (c) [2024] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Generates loadings of first factor for all models

library(MOFAcellulaR)
library(tidyverse)

# These are the training objects
mofa_models_obj <-  tibble(file = list.files("./data/pbulk/")) %>%
  dplyr::mutate(study = gsub("_pbulk[.]csv","",file),
                mofa_out = paste0("./data/models/", gsub("_pbulk[.]csv","_mofa.hdf5",file))) %>%
  dplyr::filter(study != "Armute2023_LVAD") %>%
  dplyr::select(-file) %>%
  dplyr::mutate(mofa_model = map(mofa_out, MOFA2::load_model)) %>%
  dplyr::select(study, mofa_model)

mofa_loadings <- mofa_models_obj %>%
  dplyr::mutate(F1_loadings = map(mofa_model,  MOFAcellulaR::get_geneweights, factor = "Factor1")) %>%
  dplyr::select(study, F1_loadings) %>%
  unnest(cols = c(F1_loadings))


write_csv(mofa_loadings, "./results/models/all_F1_loadings.csv")






