# Copyright (c) [2024] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Check if the multicellular signature relates to
#' myocardial infarction
#'

library(MOFAcellulaR)
library(tidyverse)
source("./code/reheat2_pilot/aesthetics.R")

# MAIN ########################################################

# Multi-cell model
factor_bplots <- "./results/predict_project/MI/MI_bplots.pdf"
mds_out <- "./results/predict_project/MI/MI_fplot.pdf"
proj_w_file <- "./results/predict_project/MI/MI_proj_w.csv"
ttest_file <- "./results/predict_project/MI/MI_ttest_res.csv"
aov_file <- "./results/predict_project/MI/MI_aov_res.csv"
factor_coords  <- "./results/predict_project/MI/MI_factor_coords.csv"

# Part 1: Create data to be projected via MOFAcell

# Data to be processed
meta_data <- read_csv("./data/val_studies/metadata/Kuppe2022_MI_metadata.csv",
                      show_col_types = FALSE)[,-1] %>%
  dplyr::mutate(patient_id = strsplit(sample_id,"-[A-Z]") %>%
                  map_chr(., ~.x[[1]])) %>%
  dplyr::filter(disease_code %in% c("NF", "MI_FZ", "MI_IZ"))

# Project
pb_data <- read_csv("./data/val_studies/pbulk/Kuppe2022_MI_pbulk.csv",
                    show_col_types = FALSE)

colnames(pb_data)[1] <- "sample_id"

pb_data <- pb_data %>%
  column_to_rownames("sample_id") %>%
  as.matrix() %>%
  t()

# Importing coldata of the matrices - ignoring not annotated cell types
coldat <- read_csv("./data/val_studies/coldata/Kuppe2022_MI_coldata.csv",
                   show_col_types = FALSE)[,-1]  %>%
  dplyr::rename(ncells = "counts") %>%
  dplyr::filter(cell_type != "none") %>%
  dplyr::filter(disease_code %in% c("NF", "MI_FZ", "MI_IZ")) %>%
  group_by(sample_id) %>%
  dplyr::mutate(sample_cells = sum(ncells)) %>%
  dplyr::filter(sample_cells >= 2000) %>%
  column_to_rownames("colname")


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

# Part 2: Project new data

model <- MOFA2::load_model("./results/meta_mofacell/metamodel_mofa.hdf5")

test_projection <- project_data(model = model,
                                test_data = MOFAcell_obj)

out_test <- test_projection[,c("Factor1", "Factor2")] %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  left_join(meta_data)

write_csv(out_test, factor_coords)

# Associations to disease
assocs_aov <- test_projection %>%
  as.data.frame(.) %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(-sample_id, names_to = "Factor", values_to = "Score") %>%
  left_join(meta_data) %>%
  dplyr::filter(Factor %in% c("Factor1", "Factor2")) %>%
  group_by(Factor) %>%
  nest() %>%
  dplyr::mutate(test_res = map(data, function(dat) {

    r_data <- dat %>%
      dplyr::select(patient_id, patient_group, Score)

    aov(Score ~ patient_group, r_data) %>% broom::tidy()

  })) %>%
  unnest(test_res) %>%
  dplyr::filter(term == "patient_group") %>%
  ungroup() %>%
  dplyr::mutate(adj_p = p.adjust(p.value)) %>%
  dplyr::select(Factor, term, p.value, adj_p)

assocs_ttest <- test_projection %>%
  as.data.frame(.) %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(-sample_id, names_to = "Factor", values_to = "Score") %>%
  left_join(meta_data) %>%
  dplyr::filter(Factor %in% c("Factor1", "Factor2")) %>%
  group_by(Factor) %>%
  nest() %>%
  dplyr::mutate(test_res = map(data, function(dat) {

    r_data <- dat %>%
      dplyr::select(patient_id, disease_code, Score)

    pairwise.t.test(r_data$Score, r_data$disease_code,
                    p.adjust.method = "BH") %>%
      broom::tidy()


  })) %>%
  unnest(test_res) %>%
  dplyr::select(-data)

write_csv(assocs_ttest, ttest_file)
write_csv(assocs_aov, aov_file)

# Plots

# Boxplots #############################################

F1_plot <- test_projection %>%
  as.data.frame(.) %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(-sample_id, names_to = "Factor", values_to = "Score") %>%
  left_join(meta_data) %>%
  dplyr::filter(Factor %in% c("Factor1")) %>%
  ggplot(aes(x = disease_code, y = Score, fill = disease_code)) +
  geom_boxplot() +
  geom_point(size = 3, alpha = 0.8) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust =0.5),
        axis.text = element_text(size =12),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = etiology_colors[meta_data$disease_code %>% unique]) +
  ylab("MCP1")

F2_plot <- test_projection %>%
  as.data.frame(.) %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(-sample_id, names_to = "Factor", values_to = "Score") %>%
  left_join(meta_data) %>%
  dplyr::filter(Factor %in% c("Factor2")) %>%
  ggplot(aes(x = disease_code, y = Score, fill = disease_code)) +
  geom_boxplot() +
  geom_point(size = 3, alpha = 0.8) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust =0.5),
        axis.text = element_text(size =12),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = etiology_colors[meta_data$disease_code %>% unique]) +
  ylim(-1.2, 1.8) +
  ylab("MCP2")

pdf(factor_bplots, width = 3.3, height = 2.5)
plot(F1_plot)
plot(F2_plot)
dev.off()

# Factor
fplot <- test_projection %>%
  as.data.frame(.) %>%
  rownames_to_column("sample_id") %>%
  left_join(meta_data) %>%
  ggplot2::ggplot(aes(x = Factor1, y = Factor2, color = disease_code)) +
  ggplot2::geom_point() +
  ggplot2::theme_classic() +
  ggplot2::theme(axis.text = ggplot2::element_text(size = 10)) +
  scale_color_manual(values = etiology_colors[meta_data$disease_code %>% unique])

pdf(mds_out, width = 3, height = 2.5)
plot(fplot)
dev.off()
