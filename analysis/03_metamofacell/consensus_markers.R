# Copyright (c) [2023] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we generate a consensus marker analysis
library(survcomp)
library(tidyverse)
setwd("/Users/ricardoramirez/Dropbox/PostDoc/Research/ReHeaT2/")
consensus_file <- "./results/consensus_markers.csv"

# Main directory
markers <- tibble(marker_csv = list.files("./data/mrkrs/")) %>%
  dplyr::mutate(study = gsub("_mrkrs[.]csv","",marker_csv)) %>%
  dplyr::mutate(marker_csv = paste0("./data/mrkrs/",marker_csv)) %>%
  dplyr::filter(study != "Reichart2022_DCM") %>%
  dplyr::mutate(markers = map(marker_csv, read_csv,
                              show_col_types = FALSE)) %>%
  dplyr::select(-marker_csv) %>%
  unnest(c(markers))

# Function to make consensus with Fisher test

run_fisher_meta <- function(pval_mat, logFC_mat, n_missing = 1){

  # Use only genes that are present in all experiments (missing in n at most)
  filt_pval_mat <- pval_mat[rowSums(is.na(pval_mat))<=n_missing,]

  # Fisher combined test
  fisher_pvals <- apply(pval_mat, 1, function(x){
    survcomp::combine.test(x, "fisher", na.rm = T)
  })

  fisher_pvals_adj <- sort(p.adjust(fisher_pvals,"BH"))

  mean_lfc <- rowMeans(logFC_mat[names(fisher_pvals_adj),],
                       na.rm = TRUE)

  consensus <- enframe(mean_lfc,
                       name = "gene",
                       value = "mean_LFC") %>%
    left_join(enframe(fisher_pvals_adj,
                      name = "gene",
                      value = "adj_Fisher_p"),
              by = "gene")

  return(consensus)
}

# Generating matrices for consensus and consensus results
markers <- markers %>%
  group_by(name) %>%
  tidyr::nest() %>%
  dplyr::mutate(pval_mat = map(data, function(dat) {

    dat %>%
      dplyr::select(study, gene, FDR) %>%
      pivot_wider(names_from = study,
                  values_from = FDR,
                  values_fill = NA) %>%
      column_to_rownames("gene") %>%
      as.matrix()


  })) %>%
  dplyr::mutate(logFC_mat = map(data, function(dat) {

    dat %>%
      dplyr::select(study, gene, logFC) %>%
      pivot_wider(names_from = study,
                  values_from = logFC,
                  values_fill = NA) %>%
      column_to_rownames("gene") %>%
      as.matrix()

  })) %>%
  dplyr::mutate(consensus_markers = map2(pval_mat, logFC_mat,
                                         run_fisher_meta, n_missing = 1)) %>%
  dplyr::select(name, consensus_markers) %>%
  tidyr::unnest(consensus_markers) %>%
  dplyr::rename("cell_type" = "name")

write_csv(markers, consensus_file)

# Maybe plot?

markers %>%
  dplyr::filter(adj_Fisher_p <= 0.0001, mean_LFC > 2) %>%
  group_by(cell_type) %>%
  summarise(n())
