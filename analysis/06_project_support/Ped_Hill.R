# Copyright (c) [2024] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Project data: test on previous data

library(MOFAcellulaR)
library(tidyverse)
library(cowplot)
source("./code/reheat2_pilot/aesthetics.R")

# Multi-cell model
mcell_reheat <- MOFA2::load_model("./results/meta_mofacell/metamodel_mofa.hdf5")

# outs
factor_stats <- "./results/predict_project/Hill/Hill_assocs.csv"
factor_bplots <- "./results/predict_project/Hill/Hill_factors.pdf"
factor_coords  <- "./results/predict_project/Hill/Hill_factor_coords.csv"

# New data
all_data <- readRDS("./data/val_studies/misc/CHD_HILL_pb.rds")

# Counts
pb_data <- all_data$pb

# Col-data
coldat <- all_data$target  %>%
  column_to_rownames("pb_id") %>%
  dplyr::rename(ncells = "cell_count",
                cell_type = "cell_type_translated") %>%
  dplyr::filter(!is.na(cell_type)) %>%
  dplyr::rename(sample = "orig_ident") %>%
  dplyr::filter(region == "LV")

pb_data <- pb_data[, rownames(coldat)]

# Meta data
meta_data <- coldat[, c("sample", "age", "gender",
                        "Diagnosis", "procedure", "DEid", "region")] %>%
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
  MOFAcellulaR::center_views(pb_dat_list = .) %>%
  MOFAcellulaR::pb_dat2MOFA(pb_dat_list = .,
                            sample_column = "sample")

# Projection in meta model:
test_projection <- project_data(model = mcell_reheat,
                                test_data = MOFAcell_obj)

out_test <- test_projection[,c("Factor1", "Factor2")] %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  left_join(meta_data, by = c("sample_id" = "sample"))

write_csv(out_test, factor_coords)

# Test association with condition

test_assocs_t <- test_projection %>%
  as.data.frame(.) %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "Factor", values_to = "Score") %>%
  left_join(meta_data) %>%
  dplyr::filter(Factor %in% c("Factor1", "Factor2")) %>%
  group_by(Factor) %>%
  nest() %>%
  mutate(pair_ttests = map(data, function(dat) {
    t.test(Score ~ DEid, dat) %>%
      broom::tidy()
  })) %>%
  dplyr::select(Factor, pair_ttests) %>%
  unnest() %>%
  ungroup()

write_csv(test_assocs_t, factor_stats)

test_assocs_aov_diag <- test_projection %>%
  as.data.frame(.) %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "Factor", values_to = "Score") %>%
  left_join(meta_data) %>%
  dplyr::filter(Factor %in% c("Factor1", "Factor2")) %>%
  group_by(Factor) %>%
  nest() %>%
  mutate(pair_ttests = map(data, function(dat) {
    aov(Score ~ Diagnosis, dat) %>%
      broom::tidy()
  })) %>%
  dplyr::select(Factor, pair_ttests) %>%
  unnest() %>%
  ungroup()

f_interest <- c("Factor1", "Factor2")

plots_diagnosis <- map(set_names(f_interest), function(fact) {

  useful_factors <- test_projection %>%
    as.data.frame(.) %>%
    rownames_to_column("sample") %>%
    pivot_longer(-sample, names_to = "Factor", values_to = "Score") %>%
    left_join(meta_data) %>%
    dplyr::filter(Factor %in% c(fact)) %>%
    ggplot(aes(x = Diagnosis, y = Score)) +
    geom_boxplot() +
    geom_point(size = 3, alpha = 0.8) +
    theme_bw() +
    theme(axis.text = element_text(size =12),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ylab(fact)

})

plots_diagnosis_simpf <- map(set_names(f_interest), function(fact) {

  useful_factors <- test_projection %>%
    as.data.frame(.) %>%
    rownames_to_column("sample") %>%
    pivot_longer(-sample, names_to = "Factor", values_to = "Score") %>%
    left_join(meta_data) %>%
    dplyr::filter(Factor %in% c(fact)) %>%
    ggplot(aes(x = DEid, y = Score)) +
    geom_boxplot() +
    geom_point(size = 3, alpha = 0.8) +
    theme_bw() +
    theme(axis.text = element_text(size =12),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ylab(fact)

})

pdf(factor_bplots, width = 2, height = 2.5)
plot(plots_diagnosis_simpf[[1]])
plot(plots_diagnosis_simpf[[2]]) + ylab("MCP2") + ylim(-.4,.6) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot(plots_diagnosis[[1]])
plot(plots_diagnosis[[2]])
dev.off()









