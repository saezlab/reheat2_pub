# Copyright (c) [2023] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we run MOFAcell models for each of the available datasets:
#' -Perform QC statistics

library(MOFAcellulaR)
library(tidyverse)
library(cowplot)
#reticulate::use_condaenv(condaenv = "/Users/ricardoramirez/mambaforge/envs/sc_base")
# Aesthetics - colors for plotting
source("./code/reheat2_pilot/aesthetics.R")

# Main ---------------------------------------------------------------
setwd("/Users/ricardoramirez/Dropbox/PostDoc/Research/ReHeaT2/")
model_stats_plt_file <- "./results/meta_mofacell/model_stats.pdf"
model_qc_plt_file <- "./results/meta_mofacell/model_input.pdf"
model_stats_file <- "./results/meta_mofacell/model_stats.csv"
mofa_out <- "./results/meta_mofacell/metamodel_mofa.hdf5"
meta_out <- "./results/meta_mofacell/metamodel_meta.csv"
mofa_pdf <- "./results/meta_mofacell/metamodel_mofa.pdf"
mds_pdf <- "./results/meta_mofacell/metamodel_mds.pdf"
mofa_loadings_file <- "./results/meta_mofacell/gene_loadings.csv"
violin_pdf <- "./results/meta_mofacell/metamodel_F1_violin.pdf"
model_stats_plt_R2file <-"./results/meta_mofacell/model_stats_R2HF.pdf"

# For association with HF we don't need these disease codes from Kuppe et. al.
input_df <- tibble(file = list.files("./data/pbulk/")) %>%
  dplyr::mutate(study = gsub("_pbulk[.]csv","",file),
                coldata_file =  paste0("./data/coldata/", gsub("_pbulk[.]csv","_coldata.csv",file)),
                metadata_file =  paste0("./data/metadata_ext/", gsub("_pbulk[.]csv","_metadata.csv",file)),
                marker_csv = paste0("./data/mrkrs/",
                                    gsub("_pbulk[.]csv","_mrkrs.csv",file))) %>%
  dplyr::mutate(file = paste0("./data/pbulk/",file)) %>%
  dplyr::filter(study %in% c("Chaffin2022_DCM", "Koenig2022_DCM",
                             "LVReichart2022_DCM", "Simonson2023_ICM"))

# We get all meta_data
all_meta <- tibble(file = list.files("./data/metadata_ext//")) %>%
  dplyr::mutate(study = gsub("_metadata[.]csv","",file)) %>%
  dplyr::filter(study %in% c("Chaffin2022_DCM", "Koenig2022_DCM",
                             "LVReichart2022_DCM", "Simonson2023_ICM")) %>%
  dplyr::mutate(file = paste0("./data/metadata_ext/",file)) %>%
  dplyr::mutate(metadata = map(file, ~ read_csv(.x, show_col_types = FALSE))) %>%
  unnest(metadata) %>%
  dplyr::select(-file)

write_csv(all_meta, meta_out)

all_samples <- all_meta$sample_id %>%
  unique() %>%
  length()

study_samples <- all_meta %>%
  dplyr::select(study, sample_id) %>%
  unique() %>%
  group_by(study) %>%
  summarize(n_samples = n())

# First we need to get all markers together [Consensus]
all_mrks <- read_csv("./results/consensus_markers.csv", show_col_types = F) %>%
  dplyr::filter(adj_Fisher_p <= 0.0001, mean_LFC > 2) %>%
  dplyr::select(cell_type, gene) %>%
  unique() %>%
  group_by(cell_type) %>%
  nest() %>%
  dplyr::mutate(data = map(data, ~.x[[1]])) %>%
  deframe()

# We need all pbulks together
# Processing will be handled for each study individually

input_df_pb <- input_df %>%
  dplyr::select(-marker_csv)

all_pbs <- pmap(input_df_pb, function(file, study, coldata_file, metadata_file) {

  print(study)

  # Importing pb data
  pb_data <- read_csv(file,
                      show_col_types = FALSE)

  colnames(pb_data)[1] <- "sample_id"

  pb_data <- pb_data %>%
    column_to_rownames("sample_id") %>%
    as.matrix() %>%
    t()

  # Importing meta_data
  meta_data <- read_csv(metadata_file,
                        show_col_types = FALSE)

  # Importing coldata of the matrices - ignoring not annotated cell types
  coldat <- read_csv(coldata_file,
                     show_col_types = FALSE)[,-1]  %>%
    column_to_rownames("colname") %>%
    dplyr::rename(ncells = "counts") %>%
    dplyr::filter(cell_type != "none") %>%
    dplyr::filter(sample_id %in% meta_data$sample_id) # This ensures that filtering happens

  pb_data <- pb_data[,rownames(coldat)]

  # Defining cts
  cts <- coldat$cell_type %>%
    unique() %>%
    set_names()

  # Defining parameters
  n_samples <- nrow(meta_data)
  min_samples <- (n_samples * 0.4) %>% floor()

  # Creating MOFAcell object to make manipulations
  MOFAcell_obj <- MOFAcellulaR::create_init_exp(pb_data, coldat) %>%
    MOFAcellulaR::filt_profiles(pb_dat = .,
                                cts = cts,
                                ncells = 20,
                                counts_col = "ncells",
                                ct_col = "cell_type") %>%
    MOFAcellulaR::filt_views_bysamples(pb_dat_list = .,
                                       nsamples = min_samples) %>%
    MOFAcellulaR::filt_gex_byexpr(pb_dat_list = .,
                                  min.count = 20,
                                  min.prop = 0.4) %>%
    MOFAcellulaR::filt_views_bygenes(pb_dat_list = .,
                                     ngenes = 50) %>%
    MOFAcellulaR::filt_samples_bycov(pb_dat_list = ., # Filtering of low quality samples
                                     prop_coverage = 0.97) %>%
    MOFAcellulaR::tmm_trns(pb_dat_list = .,
                           scale_factor = 1000000) %>%
    MOFAcellulaR::filt_gex_bybckgrnd(pb_dat_list = .,
                                     prior_mrks = all_mrks) %>%
    MOFAcellulaR::filt_views_bygenes(pb_dat_list = .,
                                     ngenes = 50) %>%
    MOFAcellulaR::pb_dat2MOFA(pb_dat_list = .,
                              sample_column = "sample_id") %>%
    dplyr::mutate(group = study)

})

all_pbs <- bind_rows(all_pbs)

# Model completion
# Generating QC of feature and sample coverage
view_completion <- all_pbs %>%
  group_by(view) %>%
  summarise(n_genes = length(unique(feature)),
            n_ind = length(unique(sample))) %>%
  dplyr::mutate(perc_samples = (n_ind/all_samples) * 100)

view_completion_study <- all_pbs %>%
  group_by(group, view) %>%
  summarise(n_genes = length(unique(feature)),
            n_ind = length(unique(sample)))  %>%
  left_join(study_samples, by = c("group" = "study")) %>%
  dplyr::mutate(perc_samples = (n_ind/n_samples) * 100) %>%
  dplyr::rename("study" = "group")

sample_qc <- ggplot(view_completion_study, aes(x = view, y = perc_samples)) +
  geom_bar(stat = "identity") +
  facet_wrap(.~study, ncol = 4) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("") +
  ylab("Percentage of samples")

genes_qc <- ggplot(view_completion_study, aes(x = view, y = n_genes)) +
  geom_bar(stat = "identity") +
  facet_wrap(.~study, ncol = 4) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("") +
  ylab("Number of genes")

# Panel with both type of stats
# -----------------------------------

qc_model_plt <- cowplot::plot_grid(sample_qc, genes_qc, ncol = 1, align = "hv")

pdf(model_qc_plt_file, height = 4.5, width = 6)

plot(qc_model_plt)

dev.off()

# Fitting the MOFA model

# Building MOFA model

MOFAobject <- create_mofa(all_pbs)

data_opts <- get_default_data_options(MOFAobject)

data_opts$center_groups <- TRUE

model_opts <- get_default_model_options(MOFAobject)

model_opts$num_factors <- 10

model_opts$spikeslab_weights <- FALSE

train_opts <- get_default_training_options(MOFAobject)

# Prepare MOFA model:
MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

# Train model:
model <- run_mofa(MOFAobject, mofa_out, use_basilisk = T)
model <- MOFA2::load_model(mofa_out)

# Variance Explained
var_explained <- model@cache$variance_explained$r2_total %>%
  map(., ~ enframe(.x, "view","R2_consensus")) %>%
  enframe("study") %>%
  unnest(c(value)) %>%
  dplyr::mutate(R2_consensus = as.numeric(R2_consensus))


# Association with HF (all studies)
# Association with HF
expl_var_HF <- MOFAcellulaR::get_associations(model = model,
                                              metadata = all_meta,
                                              sample_id_column = "sample_id",
                                              test_variable = "heart_failure",
                                              test_type = "categorical",
                                              categorical_type = "parametric",
                                              group = TRUE)

# Paper stats
expl_var_HF %>%
  dplyr::filter(adj_pvalue <= .05)

# Association per study

studies <- set_names(all_meta$study %>% unique)

assoc_list <- map(studies, function(s) {

  MOFAcellulaR::get_associations(model = model,
                                 metadata = all_meta %>%
                                   dplyr::filter(study == s),
                                 sample_id_column = "sample_id",
                                 test_variable = "heart_failure",
                                 test_type = "categorical",
                                 categorical_type = "parametric",
                                 group = TRUE)


})

enframe(assoc_list) %>%
  unnest() %>%
  dplyr::filter(adj_pvalue <= .05)

model@cache$variance_explained$r2_per_factor %>%
  map(., ~ .x[c("Factor1", "Factor2", "Factor5"),] %>%
        colSums() %>%
        mean()) %>%
  unlist() %>%
  mean()

# Get R2 associated to HF

R2_HF <- map(studies, function(s) {

  study_stat <- assoc_list[[s]]

  study_var_explained <- model@cache$variance_explained$r2_per_factor[[s]]

  useful_factors <- study_stat %>%
    dplyr::filter(adj_pvalue <= 0.05) %>%
    pull(Factor)


  covar_R2 <- study_var_explained[useful_factors,, drop = F] %>%
    colSums() %>%
    enframe("view", "R2_HF_consensus")

  return(covar_R2)

}) %>%
  enframe("study") %>%
  unnest(c(value))

# Summarize all the model stats

model_stats <- left_join(var_explained, R2_HF, by = c("study", "view")) %>%
  left_join(view_completion_study, by = c("study", "view"))

write_csv(model_stats, model_stats_file)

# Visualize?

assoc_list[["All"]] <- expl_var_HF

col_list <- list(heart_failure = c("NF" = "darkgrey",
                                   "HF" = "black"),
                 disease_code = etiology_colors[all_meta$disease_code %>% unique],
                 study = col_pub$study_colors_sc[all_meta$study %>% unique])

scores_hmap <- plot_MOFA_hmap(model = model,
                                            group = TRUE,
                                            metadata = all_meta,
                                            sample_id_column = "sample_id",
                                            sample_anns = c("heart_failure", "disease_code",
                                                            "study"),
                                            assoc_list = assoc_list,
                                            col_rows = col_list)


pdf(mofa_pdf, height = 12, width = 7)

draw(scores_hmap)

dev.off()


# MDS
# Generate MDS
MDS_dat <- MOFAcellulaR::plot_sample_2D(model,
                                        metadata = all_meta,
                                        color_by = "heart_failure",
                                        sample_id_column = "sample_id",
                                        group = T,method = "MDS")

MDS_plt <- last_plot()

pdf(mds_pdf, height = 3, width = 4.2)

plot(MDS_plt +
       scale_color_manual(values = c("black","darkgrey")) +
       coord_equal())

dev.off()

# R2 associated with disease
ct_order <- model_stats %>%
  group_by(view) %>%
  summarize(med_r2 = median(R2_HF_consensus)) %>%
  arrange(desc(med_r2)) %>%
  pull(view)

study_order <- model_stats %>%
  group_by(study) %>%
  summarize(med_r2 = mean(R2_HF_consensus)) %>%
  arrange(desc(med_r2)) %>%
  pull(study)

ct_r2_hf_plt <- model_stats %>%
  dplyr::mutate(view = factor(view, levels = ct_order)) %>%
  ggplot(aes(x = view, y = R2_HF_consensus)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("") +
  ylab("R2 associated \n with HF")

study_r2_hf_plt <- model_stats %>%
  dplyr::mutate(study = factor(study, levels = study_order)) %>%
  ggplot(aes(x = study, y = R2_HF_consensus)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("") +
  ylab("R2 associated \n with HF")


R2_disease <- cowplot::plot_grid(plotlist = list(ct_r2_hf_plt, study_r2_hf_plt),
                                 nrow = 1, align = "hv", rel_widths = c(1,0.65))

pdf(model_stats_plt_R2file, height = 3, width = 4)

plot(R2_disease)

dev.off()

# WRITING LOADINGS

F_names <- set_names(factors_names(model))

all_loadings <- map(F_names, function(factor_name){
  MOFAcellulaR::get_geneweights(model, factor_name)
}) %>%
  enframe("Factor") %>%
  unnest()

write_csv(all_loadings, mofa_loadings_file)
