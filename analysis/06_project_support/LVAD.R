# Copyright (c) [2023] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Check if the multicellular signature relates to
#' cardiac reverse remodelling
#'

# Main ##############################################

library(MOFAcellulaR)
library(tidyverse)
library(cowplot)
source("./code/reheat2_pilot/aesthetics.R")

# Multi-cell model
mcell_reheat <- MOFA2::load_model("./results/meta_mofacell/metamodel_mofa.hdf5")
mds_out <- "./results/predict_project/LVAD/LVAD_fplot.pdf"
assoc_rec_out <- "./results/predict_project/LVAD/LVAD_ttest_rec_res.csv"
assoc_prog_out <- "./results/predict_project/LVAD/LVAD_ttest_prog_res.csv"
factor_bplots <- "./results/predict_project/LVAD/LVAD_bplots.pdf"
proj_w_file  <- "./results/predict_project/LVAD/LVAD_proj_w.csv"
factor_coords  <- "./results/predict_project/LVAD/LVAD_factor_coords.csv"

# Data to be processed
meta_data <- read_csv("./data/val_studies/metadata/Armute2023_LVAD_metadata.csv",
                      show_col_types = FALSE)[,-1] %>%
  dplyr::mutate(patient_id = strsplit(sample_id,"-[A-Z]") %>%
                  map_chr(., ~.x[[1]]))

# Project data to multicellular space

pb_data <- read_csv("./data/val_studies/pbulk/Armute2023_LVAD_pbulk.csv",
                    show_col_types = FALSE)

colnames(pb_data)[1] <- "sample_id"

pb_data <- pb_data %>%
  column_to_rownames("sample_id") %>%
  as.matrix() %>%
  t()

# Importing coldata of the matrices - ignoring not annotated cell types
coldat <- read_csv("./data/val_studies/coldata/Armute2023_LVAD_coldata.csv",
                   show_col_types = FALSE)[,-1]  %>%
  column_to_rownames("colname") %>%
  dplyr::rename(ncells = "counts") %>%
  dplyr::filter(cell_type != "none")

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
  MOFAcellulaR::tmm_trns(pb_dat_list = .,
                         scale_factor = 1000000) %>%
  MOFAcellulaR::center_views(pb_dat_list = .) %>%
  MOFAcellulaR::pb_dat2MOFA(pb_dat_list = .,
                            sample_column = "sample_id")

test_projection <- project_data(model = mcell_reheat,
                                test_data = MOFAcell_obj)

out_test <- test_projection[,c("Factor1", "Factor2")] %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  left_join(meta_data)

write_csv(out_test, factor_coords)

# Test for association with recovery status
# Do recovered patients see a difference in factors?

assoc_recovery <- test_projection %>%
  as.data.frame(.) %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(-sample_id, names_to = "Factor", values_to = "Score") %>%
  left_join(meta_data) %>%
  dplyr::filter(Factor %in% c("Factor1", "Factor2")) %>%
  dplyr::filter(heart_failure == "HF") %>%
  group_by(Factor, response) %>%
  nest() %>%
  dplyr::mutate(test_res = map(data, function(dat) {

    r_data <- dat %>%
      dplyr::select(patient_id, biopsy, Score) %>%
      pivot_wider(values_from = Score, names_from = biopsy)

    t.test(x = r_data$post, y = r_data$pre, paired = T,
           alternative = "two.sided") %>%
      broom::tidy()

  })) %>%
  unnest(test_res) %>%
  group_by(response) %>%
  dplyr::mutate(adj_p = p.adjust(p.value)) %>%
  dplyr::select(Factor, response, adj_p)

write_csv(assoc_recovery, assoc_rec_out)

# Test for association in pre-samples
# Do we see a difference between recovered and not recovered?
assoc_prognosis <- test_projection %>%
  as.data.frame(.) %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(-sample_id, names_to = "Factor", values_to = "Score") %>%
  left_join(meta_data) %>%
  dplyr::filter(Factor %in% c("Factor1", "Factor2")) %>%
  dplyr::filter(heart_failure == "HF") %>%
  group_by(Factor, biopsy) %>%
  nest() %>%
  dplyr::mutate(test_res = map(data, function(dat) {

    t.test(Score ~ response, data = dat,
           alternative = "two.sided") %>% broom::tidy()

  })) %>%
  unnest(test_res) %>%
  group_by(biopsy) %>%
  dplyr::mutate(adj_p = p.adjust(p.value)) %>%
  dplyr::select(Factor, biopsy, adj_p)

write_csv(assoc_prognosis, assoc_prog_out)

# Plot of axis of recovery

# First comparison: F1 separates pre-post

bplot <- test_projection[, "Factor1", drop = F] %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  left_join(meta_data, by = "sample_id") %>%
  dplyr::filter(heart_failure == "HF") %>%
  dplyr::mutate(biopsy = factor(biopsy, levels = c("pre", "post"))) %>%
  ggplot(aes(x = biopsy, y = Factor1, fill = response)) +
  geom_boxplot() +
  #geom_point(size = 3, alpha = 0.8) +
  theme_bw() +
  theme(axis.text = element_text(size =12),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c("#387780","darkslategray3")) +
  ylab("MCP1") +
  xlab("Biopsy type")


pdf(factor_bplots, width = 3.5, height = 2.5)
plot(bplot)
dev.off()

# Distribution across factors

fplot <- test_projection %>%
  as.data.frame(.) %>%
  rownames_to_column("sample_id") %>%
  left_join(meta_data) %>%
  ggplot2::ggplot(aes(x = Factor1, y = Factor2, color = disease_code)) +
  ggplot2::geom_point() +
  ggplot2::theme_classic() +
  ggplot2::theme(axis.text = ggplot2::element_text(size = 10)) +
  scale_color_manual(values = etiology_colors[meta_data$disease_code %>% unique])

pdf(mds_out, height = 2, width = 3.5)

plot(fplot)

dev.off()

# Generate projections

# projection_weights <- map(set_names(c("Factor1", "Factor2")), function(Fnum) {
#
#   F_weights <- get_projection_weights(model = mcell_reheat,
#                                       test_data = MOFAcell_obj,
#                                       sample_id_column = "sample_id",
#                                       factor = Fnum)
#
# }) %>%
#   bind_rows() %>%
#   left_join(meta_data)
#
# projection_weights %>%
#   group_by(ctype, gene, factor_name, condition) %>%
#   dplyr::summarise(med_w = median(value)) %>%
#   arrange(factor_name, desc(abs(med_w))) %>%
#   write_csv(proj_w_file)




