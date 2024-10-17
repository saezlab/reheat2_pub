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
model_stats_plt_file <- "./results/models/model_stats.pdf"
model_stats_file <- "./results/models/model_stats.csv"
model_stats_plt_R2file <- "./results/models/model_stats_R2HF.pdf"


# For association with HF we don't need these disease codes from Kuppe et. al.
input_df <- tibble(file = list.files("./data/pbulk/")) %>%
  dplyr::mutate(study = gsub("_pbulk[.]csv","",file),
                coldata_file =  paste0("./data/coldata/", gsub("_pbulk[.]csv","_coldata.csv",file)),
                metadata_file =  paste0("./data/metadata_ext/", gsub("_pbulk[.]csv","_metadata.csv",file)),
                marker_csv = paste0("./data/mrkrs/",
                                    gsub("_pbulk[.]csv","_mrkrs.csv",file)),
                mofa_out = paste0("./data/models/", gsub("_pbulk[.]csv","_mofa.hdf5",file)),
                mofa_pdf = paste0("./results/models/", gsub("_pbulk[.]csv","_mofa.pdf",file)),
                factor_pdf = paste0("./results/models/", gsub("_pbulk[.]csv","_mofa_factor.pdf",file)),
                umap_pdf = paste0("./results/models/", gsub("_pbulk[.]csv","_mofa_umap.pdf",file))) %>%
  dplyr::mutate(file = paste0("./data/pbulk/",file))

# Run MOFAcell models

model_outs <- pmap(input_df, function(study, file,
                                      coldata_file, metadata_file,
                                      marker_csv, mofa_out,
                                      mofa_pdf, factor_pdf, umap_pdf) {

  print(study)

  # Importing pb data
  pb_data <- read_csv(file,
                      show_col_types = FALSE)

  colnames(pb_data)[1] <- "sample_id"

  pb_data <- pb_data %>%
    column_to_rownames("sample_id") %>%
    as.matrix() %>%
    t()

  # Importing coldata of the matrices - ignoring not annotated cell types
  coldat <- read_csv(coldata_file,
                     show_col_types = FALSE)[,-1]  %>%
    column_to_rownames("colname") %>%
    dplyr::rename(ncells = "counts") %>%
    dplyr::filter(cell_type != "none")

  pb_data <- pb_data[,rownames(coldat)]

  # Importing meta_data

  meta_data <- read_csv(metadata_file,
                        show_col_types = FALSE)

  # Importing markers
  mrkrs <- read_csv(marker_csv, show_col_types = FALSE) %>%
    dplyr::filter(FDR < 0.01, logFC >= 2) %>%
    dplyr::rename("cell_type" = name) %>%
    group_by(cell_type) %>%
    nest() %>%
    dplyr::mutate(data = map(data, ~.x[[1]])) %>%
    deframe()

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
                                     prior_mrks = mrkrs) %>%
    MOFAcellulaR::filt_views_bygenes(pb_dat_list = .,
                                     ngenes = 50) %>%
    MOFAcellulaR::pb_dat2MOFA(pb_dat_list = .,
                              sample_column = "sample_id")

  # Generating QC of feature and sample coverage
  view_completion <- MOFAcell_obj %>%
    group_by(view) %>%
    summarise(n_genes = length(unique(feature)),
              n_ind = length(unique(sample))) %>%
    dplyr::mutate(perc_samples = (n_ind/n_samples) * 100)


  # Building MOFA model

  MOFAobject <- create_mofa(MOFAcell_obj)

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
  model <- run_mofa(MOFAobject, mofa_out,use_basilisk = T)

  # Summarize model

  # R2 - of model
  var_explained <- model@cache$variance_explained$r2_total[[1]] %>%
    enframe("view","R2")

  # Association with HF
  expl_var_HF <- MOFAcellulaR::get_associations(model = model,
                                                metadata = meta_data,
                                                sample_id_column = "sample_id",
                                                test_variable = "disease_code",
                                                test_type = "categorical",
                                                categorical_type = "parametric",
                                                group = FALSE)

  # Complete summary table
  useful_factors <- expl_var_HF %>%
    dplyr::filter(adj_pvalue < 0.05) %>%
    pull(Factor)

  if(rlang::is_empty(useful_factors)) {
    var_explained <- var_explained %>%
      dplyr::mutate(R2_HF = 0)

  } else {
    R2_HF <- model@cache$variance_explained$r2_per_factor[[1]][useful_factors,, drop = F] %>%
      colSums() %>%
      enframe("view","R2_HF")

    var_explained <- var_explained %>%
      left_join(R2_HF, by = "view")

  }

  model_summ <- left_join(view_completion,
                          var_explained,
                          by = "view") %>%
    dplyr::mutate(study = study)

  # Visualize?

  assoc_list <- list(disease = expl_var_HF)

  col_list <- list(heart_failure = c("NF" = "darkgrey",
                                "HF" = "black"),
                   disease_code = etiology_colors[meta_data$disease_code %>% unique])

  scores_hmap <- MOFAcellulaR::plot_MOFA_hmap(model = model,
                                              group = FALSE,
                                              metadata = meta_data,
                                              sample_id_column = "sample_id",
                                              sample_anns = c("heart_failure", "disease_code"),
                                              assoc_list = assoc_list,
                                              col_rows = col_list)


  pdf(mofa_pdf, height = 8, width = 4.5)

  draw(scores_hmap)

  dev.off()


  # Generate UMAP
  UMAP_dat <- MOFAcellulaR::plot_sample_2D(model,
                                           method = "UMAP",
                                           metadata = meta_data,
                                           color_by = "disease_code",
                                           sample_id_column = "sample_id")

  UMAP_plt <- last_plot()


  pdf(umap_pdf, height = 3, width = 3.9)

  plot(UMAP_plt +
         scale_color_manual(values = etiology_colors[meta_data$disease_code %>% unique])  +
         ggtitle(study) )

  dev.off()

  # Etiology
  main_facts_plt_etiol <- MOFAcellulaR::get_tidy_factors(model = model,
                                                         factor = "all",
                                                         metadata = meta_data,
                                                         group = F,
                                                         sample_id_column = "sample_id") %>%
    dplyr::filter(Factor %in% c("Factor1", "Factor2")) %>%
    pivot_wider(values_from = value, names_from = Factor) %>%
    ggplot(aes(x = Factor1, y = Factor2, color = disease_code)) +
    geom_point(size = 1.5) +
    theme_bw() +
    theme(axis.text = element_text(size =10)) +
    scale_color_manual(values = etiology_colors[meta_data$disease_code %>% unique])


  pdf(factor_pdf, height = 3, width = 3.9)

  plot(main_facts_plt_etiol)

  dev.off()



  return(model_summ)

})

# Generate final map

model_stats <- enframe(model_outs) %>%
  unnest(c(value)) %>%
  dplyr::select(-name) %>%
  dplyr::mutate(R2 = as.numeric(R2))

write_csv(model_stats, model_stats_file)

# Cell type

R2_ct_plt <- ggplot(model_stats %>%
                      filter(study !="Reichart2022_DCM"), aes(x = view, y = R2_HF)) +
  geom_boxplot() +
  geom_point() +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("") +
  ylab("R2 associated \n with HF")

R2_study_plt <- ggplot(model_stats %>%
                         filter(study !="Reichart2022_DCM"), aes(x = study, y = R2_HF)) +
  geom_boxplot() +
  geom_point() +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("") +
  ylab("R2 associated \n with HF")

R2_disease <- cowplot::plot_grid(plotlist = list(R2_ct_plt, R2_study_plt), nrow = 1, align = "hv", rel_widths = c(1,0.65))

pdf(model_stats_plt_R2file, height = 3, width = 4)

plot(R2_disease)

dev.off()

R2_qc <- model_stats %>%
  ggplot(aes(x = study, y = R2)) +
  geom_bar(stat = "identity") +
  facet_wrap(. ~ view, nrow = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("") +
  ylab("Model R2")

genes_qc <- model_stats %>%
  ggplot(aes(x = study, y = n_genes)) +
  geom_bar(stat = "identity") +
  facet_wrap(. ~ view, nrow = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("") +
  ylab("Number of genes")

sample_qc <- model_stats %>%
  ggplot(aes(x = study, y = perc_samples)) +
  geom_bar(stat = "identity") +
  facet_wrap(. ~ view, nrow = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("") +
  ylab("% profiled samples")

R2_hf_qc <- model_stats %>%
  ggplot(aes(x = study, y = R2_HF)) +
  geom_bar(stat = "identity") +
  facet_wrap(. ~ view, nrow = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("") +
  ylab("Model R2 disease")

model_stats_plt <- cowplot::plot_grid(genes_qc, sample_qc, R2_qc, R2_hf_qc, ncol = 2, nrow = 2, align = "hv")

pdf(model_stats_plt_file, height = 8, width = 10)

plot(model_stats_plt)

dev.off()

#

model_stats <- read_csv(model_stats_file)

model_stats %>%
  filter(study !="Reichart2022_DCM") %>%
  pull(R2_HF) %>%
  mean()

model_stats %>%
  filter(study !="Reichart2022_DCM") %>%
  group_by(view) %>%
  summarize(mean(R2_HF))


model_stats %>%
  filter(study !="Reichart2022_DCM") %>%
  dplyr::mutate(reference = ifelse(view == 'Fib', "ref", "rest")) %>%
  t.test(R2_HF ~ reference, alternative = "greater", data = .) %>%
  broom::tidy()

model_stats %>%
  filter(study !="Reichart2022_DCM") %>%
  dplyr::mutate(reference = ifelse(view == 'Lymphoid', "ref", "rest")) %>%
  t.test(R2_HF ~ reference, alternative = "less", data = .) %>%
  broom::tidy()





