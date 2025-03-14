# Copyright (c) [2023] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we calculate the general QC of distinct studies
#' This includes number of cells, genes, and contamination score

library(tidyverse)
library(MOFAcellulaR)
library(cowplot)

# This calculates the number of non-zero genes
nview_genes <- function(pb_dat_list) {

  pb_dat_red <- purrr::map(pb_dat_list, function(x) {

    ct_mat <- SummarizedExperiment::assay(x, "counts")

    n_genes_df <- colSums(ct_mat != 0) %>%
      enframe(name = "colname", value = "n_genes") %>%
      column_to_rownames("colname")

    colData(x)[,"psbulk_n_genes"] <-  n_genes_df[rownames(colData(x)), "n_genes"]

    return(x)
  })

  return(pb_dat_red)
}

# This calculates the contamination score
get_contamination <- function(pb_dat_list, mrkrs_df, cell_type_id = "cell_type") {

  # get ct names from the actual data
  analyzed_cts <- names(pb_dat_list)
  # keep markers only of analyzed data
  mrkrs_filt <- mrkrs_df[mrkrs_df[,cell_type_id][[1]] %in% analyzed_cts,]

  pb_dat_red <- purrr::map(pb_dat_list, function(x) {

    # Get ct from the object
    obj_ct <- colData(x)[,"cell_type"] %>% unique()
    cont_cts <- analyzed_cts[!analyzed_cts %in% obj_ct]

    ct_mrkrs <- mrkrs_filt[mrkrs_filt[cell_type_id][[1]] %in% obj_ct, ] %>% pull(gene)

    cont_mrks <- mrkrs_filt[mrkrs_filt[cell_type_id][[1]] %in% cont_cts, ] %>% pull(gene)

    cont_mrks <- cont_mrks[!(cont_mrks %in% ct_mrkrs)]

    ct_mat <- SummarizedExperiment::assay(x, "counts")

    ct_reads <- colSums(ct_mat[ct_mrkrs, ])
    cont_reads <- colSums(ct_mat[cont_mrks,] )

    cont_score <- cont_reads/ct_reads

    colData(x)[,"contamination"] <-  cont_score[rownames(colData(x))]

    cont_score_norm <- (cont_reads/length(cont_mrks))/(ct_reads/length(ct_mrkrs))

    colData(x)[,"norm_contamination"] <-  cont_score_norm[rownames(colData(x))]

    # Now for each ct we need the contribution in the contamination score

    cont_summary <- map(set_names(cont_cts), function(cont_ct) {

      cont_ct_mrkrs <- mrkrs_filt[mrkrs_filt[cell_type_id][[1]] %in% cont_ct, ] %>% pull(gene)
      cont_ct_reads <- colSums(ct_mat[cont_ct_mrkrs, ])
      ct_contribution <- cont_ct_reads/cont_reads

    }) %>%
      bind_rows(.,.id = "cont_ct") %>%
      dplyr::mutate(cont_ct = paste0("cont_", cont_ct)) %>%
      as.data.frame() %>%
      column_to_rownames("cont_ct") %>%
      as.matrix() %>%
      t() %>%
      as.data.frame()

    colData(x) <- cbind(colData(x), DataFrame(cont_summary)[rownames(colData(x)),])


    return(x)
  })

  return(pb_dat_red)

}


# Main ---------------------------------------------------------------
setwd("/Users/ricardoramirez/Dropbox/PostDoc/Research/ReHeaT2/")
outplot_file <- "./results/qc/qc_stats.pdf"
outqc_file <- "./results/qc/qc_stats.csv"

input_df <- tibble(file = list.files("./data/pbulk/")) %>%
  dplyr::mutate(study = gsub("_pbulk[.]csv","",file),
                coldata_file =  paste0("./data/coldata/", gsub("_pbulk[.]csv","_coldata.csv",file)),
                marker_csv = paste0("./data/mrkrs/",
                                    gsub("_pbulk[.]csv","_mrkrs.csv",file))) %>%
  dplyr::mutate(file = paste0("./data/pbulk/",file)) %>%
  dplyr::select(study, file, coldata_file, marker_csv) %>%
  dplyr::filter(study != "Reichart2022_DCM")
#study = "Chaffin2022_DCM"
#file = "./data/pbulk/Chaffin2022_DCM_pbulk.csv"
#coldata_file = "./data/coldata/Chaffin2022_DCM_coldata.csv"
#marker_csv = "./data/mrkrs/Chaffin2022_DCM_mrkrs.csv"

qc_list <- pmap(input_df, function(study, file, coldata_file, marker_csv) {

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

  # Importing markers
  mrkrs <- read_csv(marker_csv, show_col_types = FALSE) %>%
    dplyr::filter(FDR < 0.01, logFC > 2) %>%
    dplyr::rename("cell_type" = name)

  # Defining cts
  cts <- coldat$cell_type %>%
    unique() %>%
    set_names()

  # Creating MOFAcell object to make manipulations
  MOFAcell_obj <- MOFAcellulaR::create_init_exp(pb_data, coldat) %>%
    MOFAcellulaR::filt_profiles(pb_dat = .,
                                cts = cts,
                                ncells = 0,
                                counts_col = "ncells",
                                ct_col = "cell_type") %>%
    nview_genes(pb_dat_list = .) %>%
    get_contamination(pb_dat_list = .,
                      mrkrs_df = mrkrs,
                      cell_type_id = "cell_type")


  # Getting the QC df
  qc_df <- map(MOFAcell_obj, function(x) {
    colData(x) %>% as.data.frame()
  }) %>%
    enframe() %>%
    unnest(c(value))

  cont_columns <- colnames(qc_df)[grepl("cont_",colnames(qc_df))]

  qc_columns <- c("cell_type", "disease_code",
                  "heart_failure", "psbulk_n_cells",
                  "psbulk_counts", "psbulk_n_genes",
                  "contamination", "norm_contamination")


  qc_df <- qc_df %>%
    dplyr::select_at(c(qc_columns,cont_columns)) %>%
    dplyr::mutate(study_id = study)

  return(qc_df)

})


all_qc <- enframe(qc_list) %>%
  unnest(cols = c(value)) %>%
  dplyr::select(-name)

write_csv(all_qc, outqc_file)

# Number of cells
ncells_plt <- all_qc %>%
  ggplot(., aes(y = study_id, x = log10(psbulk_n_cells))) +
  geom_boxplot() +
  facet_grid(heart_failure ~ cell_type) +
  theme_bw() +
  ylab("") +
  xlab("log10(number of cells)")

# Number of reads
depth_plt <- all_qc %>%
  ggplot(., aes(y = study_id, x = log10(psbulk_counts))) +
  geom_boxplot() +
  facet_grid(heart_failure ~ cell_type) +
  theme_bw() +
  ylab("") +
  xlab("log10(number of reads)")

# Number of genes
ngenes_plt <- all_qc %>%
  ggplot(., aes(y = study_id, x = round(log10(psbulk_n_genes),1))) +
  geom_boxplot() +
  facet_grid(heart_failure ~ cell_type) +
  theme_bw() +
  ylab("") +
  xlab("log10(number of genes)") +
  scale_x_continuous(breaks=c(3.5,4,4.5))

# Contamination
cont_plt <- all_qc %>%
  ggplot(., aes(y = study_id, x = log1p(contamination))) +
  geom_boxplot() +
  facet_grid(heart_failure ~ cell_type) +
  geom_vline(xintercept = log1p(1)) +
  theme_bw() +
  ylab("") +
  xlab("Contamination") +
  scale_x_continuous(breaks=c(0,1,2))

norm_cont_plt <- all_qc %>%
  ggplot(., aes(y = study_id, x = log1p(norm_contamination))) +
  geom_boxplot() +
  facet_grid(heart_failure ~ cell_type) +
  geom_vline(xintercept = log1p(1)) +
  theme_bw() +
  ylab("") +
  xlab("Norm. Contamination")

# Contamination fraction

cont_columns <- colnames(all_qc)[grepl("cont_",colnames(all_qc))]

cont_frac_plt <- all_qc %>%
  dplyr::select_at(c("cell_type","heart_failure", "study_id", cont_columns)) %>%
  pivot_longer(-c(cell_type, heart_failure, study_id), names_to = "cont_ct", values_to = "cont_fraction") %>%
  dplyr::mutate(cont_ct = gsub("cont_","", cont_ct)) %>%
  ggplot(aes(x = cont_ct, y = cont_fraction)) +
  geom_boxplot() +
  facet_grid(heart_failure ~ cell_type) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("contamination fraction") +
  xlab("contaminating cell")


qc_plts <- cowplot::plot_grid(ncells_plt, depth_plt,
                   ngenes_plt, cont_plt, cont_frac_plt,
                   nrow = 5, ncol = 1, align = "hv")

pdf(outplot_file, height = 13, width = 8.5)
plot(qc_plts)
dev.off()


# This tests background correction

all_qc %>%
  mutate(bckground = ifelse(study_id %in% c("Simonson2023_ICM", "Chaffin2022_DCM"),
                            "yes", "no")) %>%
  group_by(cell_type) %>%
  nest() %>%
  mutate(comp_res = map(data, function(dat) {

    t.test(log1p(contamination) ~ bckground, data= dat) %>%
      broom::tidy()

  })) %>%
  unnest(comp_res) %>%
  ungroup() %>%
  dplyr::mutate(adj_p = p.adjust(p.value)) %>%
  dplyr::select(cell_type, statistic, adj_p) %>%
  dplyr::filter(adj_p <= 0.05)

# Making contamination plots

cellbender_effect_plt <- all_qc %>%
  mutate(bckground = ifelse(study_id %in% c("Simonson2023_ICM", "Chaffin2022_DCM"),
                            "yes", "no")) %>%
  #dplyr::filter(cell_type == "Endo" ) %>%
  ggplot(aes(x = bckground, y = log1p(contamination))) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(~ cell_type,scales = "free_y", ncol = 2) +
  theme_bw() +
  xlab("background correction")

pdf("./results/qc/background_correction.pdf", height = 4.5, width = 2.7)

plot(cellbender_effect_plt)

dev.off()


cont_by_study_plt <- ggplot(all_qc, aes(y = study_id, x = log1p(contamination))) +
  geom_boxplot() +
  geom_vline(xintercept = log1p(1)) +
  ylab("") + theme_bw()

cont_by_ct_plt <- ggplot(all_qc, aes(y = cell_type, x = log1p(contamination))) +
  geom_boxplot() +
  geom_vline(xintercept = log1p(1)) +
  ylab("") + theme_bw()

cont_dist <- cowplot::plot_grid(cont_by_study_plt, cont_by_ct_plt,
                   ncol = 1, rel_heights = c(.57,1),
                   align = "hv")

pdf("./results/qc/contamination_dists.pdf", height = 3.5, width = 3.5)

plot(cont_dist)

dev.off()


##################






