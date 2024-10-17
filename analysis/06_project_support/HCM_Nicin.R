# Copyright (c) [2024] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Project data: test on previous data

library(MOFAcellulaR)
library(tidyverse)
library(cowplot)
source("./code/reheat2_pilot/aesthetics.R")

projection_table <- "./results/predict_project/HCA_reference/HCM_Nicin.csv"
projection_table_2join <- "./results/predict_project/HCA_reference/HCM_Nicin_factor_coords.csv"

# Multi-cell metamodel
mcell_reheat <- MOFA2::load_model("./results/meta_mofacell/metamodel_mofa.hdf5")

# HCA reference
HCA_reference <- read_csv("./data/val_studies/misc/HCA_ref_genes.csv")

# New data
all_data <- readRDS("./data/val_studies/misc/HCM_dimmeler_pb.rds")

# Counts
pb_data <- all_data$pb

# Col-data
coldat <- all_data$target  %>%
  column_to_rownames("pb_id") %>%
  dplyr::rename(ncells = "cell_count",
                cell_type = "cell_type_translated",
                sample = "orig.ident") %>%
  dplyr::filter(!is.na(cell_type))

pb_data <- pb_data[, rownames(coldat)]

# Meta data
meta_data <- coldat[, c("sample", "condition")] %>%
  unique() %>%
  as_tibble()

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
  MOFAcellulaR::tmm_trns(pb_dat_list = .,
                         scale_factor = 1000000) %>%
  MOFAcellulaR::pb_dat2MOFA(pb_dat_list = .,
                            sample_column = "sample") %>%
  left_join(HCA_reference, by = "feature") %>%
  na.omit() %>%
  dplyr::mutate(value = (value - mu_HCA)/sd_HCA) %>%
  dplyr::select(-mu_HCA, -sd_HCA)


test_projection <- project_data(model = mcell_reheat,
                                test_data = MOFAcell_obj)

test_projection[,c("Factor1", "Factor2")]  %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  left_join(meta_data, by = c("sample_id" = "sample")) %>%
  write_csv(projection_table_2join)


test_projection <- test_projection %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample,
               names_to =  "Factor",
               values_to = "Factor_score") %>%
  dplyr::filter(Factor %in% c("Factor1", "Factor2")) %>%
  dplyr::mutate(study_label = "HCM_Nicin")

write_csv(test_projection, projection_table)
