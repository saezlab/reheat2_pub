# Copyright (c) [2024] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Generate reference of HCA

library(MOFAcellulaR)
library(tidyverse)

HCA <- readRDS("./data/val_studies/misc/HCA_pb.rds")

pb_data <- HCA$pb

coldat <- HCA$target  %>%
  column_to_rownames("pb_id") %>%
  dplyr::rename(ncells = "cell_count",
                cell_type = "cell_type_translated") %>%
  dplyr::filter(!is.na(cell_type))

pb_data <- pb_data[,rownames(coldat)]

# Processing the test data
# Defining cts
cts <- coldat$cell_type %>%
  unique() %>%
  set_names()

# Creating MOFAcell object to make projection
MOFAcell_obj <- MOFAcellulaR::create_init_exp(pb_data, coldat) %>%
  MOFAcellulaR::filt_profiles(pb_dat = .,
                              cts = cts,
                              ncells = 20,
                              counts_col = "ncells",
                              ct_col = "cell_type") %>%
  MOFAcellulaR::filt_gex_byexpr(pb_dat_list = .,
                                min.count = 20,
                                min.prop = 0.4) %>%
  MOFAcellulaR::tmm_trns(pb_dat_list = .,
                         scale_factor = 1000000) %>%
  MOFAcellulaR::pb_dat2MOFA(pb_dat_list = .,
                            sample_column = "donor") %>%
  group_by(feature) %>%
  dplyr::summarise(mu_HCA = mean(value),
                   sd_HCA = sd(value))

write_csv(MOFAcell_obj, "./data/val_studies/misc/HCA_ref_genes.csv")







