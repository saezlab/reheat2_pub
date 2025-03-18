library(MOFAcellulaR)
library(tidyverse)
library(cowplot)
library(ComplexHeatmap)


# Main ---------------------------------------------------------------
#setwd("/Users/ricardoramirez/Dropbox/PostDoc/Research/ReHeaT2/")
# Aesthetics - colors for plotting
source("./code/reheat2_pilot/aesthetics.R")

# For association with HF we don't need these disease codes from Kuppe et. al.
input_df <- tibble(file = list.files("./data/pbulk/")) %>%
  dplyr::mutate(study = gsub("_pbulk[.]csv","",file),
                coldata_file =  paste0("./data/coldata/",
                                       gsub("_pbulk[.]csv","_coldata.csv",file)),
                metadata_file =  paste0("./data/metadata_ext/",
                                        gsub("_pbulk[.]csv","_metadata.csv",file)),
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
    MOFAcellulaR::center_views(pb_dat_list = .) %>%
    MOFAcellulaR::pb_dat2MOFA(pb_dat_list = .,
                              sample_column = "sample_id") %>%
    dplyr::mutate(group = study)

})

# Continue with the gene description
all_pbs <- bind_rows(all_pbs)

# Check genes
loadings <- read_csv("./results/meta_mofacell/gene_loadings.csv") %>%
  dplyr::filter(Factor %in% c("Factor1", "Factor2"))

# Check coordination of top genes in Fib, in F1 and F2

for(Fact in c("Factor1","Factor2")) {

  test_feats <- loadings %>%
    dplyr::mutate(feature = paste0(ctype,"_",feature)) %>%
    dplyr::arrange(Factor, ctype, desc(abs(value))) %>%
    dplyr::group_by(Factor, ctype) %>%
    dplyr::slice(1:100) %>%
    dplyr::filter(Factor == Fact)

  complete_feats <- all_pbs %>%
    dplyr::filter(feature %in% test_feats$feature) %>%
    dplyr::select(feature, group) %>%
    unique() %>%
    group_by(feature) %>%
    summarize(n_studies = n()) %>%
    dplyr::filter(n_studies == 4,
                  !grepl("LINC0",feature))

  test_feats <- test_feats %>%
    dplyr::filter(feature %in% complete_feats$feature) %>%
    dplyr::slice(1:5)

  ex_genes <- all_pbs %>%
    dplyr::filter(feature %in% test_feats$feature) %>%
    dplyr::rename("sample_id" = "sample",
                  "study" = "group") %>%
    left_join(all_meta, by = c("sample_id", "study")) %>%
    group_by(view, feature, heart_failure, study) %>%
    summarise(mean_expr = mean(value)) %>%
    dplyr::mutate(feature = strsplit(feature,"_") %>%
                    map_chr(., ~ .x[[2]]),
                  colname = paste0(heart_failure, study))

  if(Fact == "Factor1") {

    write_csv(ex_genes, "./Revision/figures/Figure3/Figure3C.csv")

  }

  cts <- ex_genes$view %>% unique() %>% set_names()

  col_fun <- colorRamp2(c(min(ex_genes$mean_expr), 0, max(ex_genes$mean_expr)),
                        c("#1C1C5C", "white", "#d62728"))

  expr_mats <- map(cts, function(ct) {

    expr_mat <- ex_genes %>%
      dplyr::filter(view == ct) %>%
      ungroup() %>%
      dplyr::select(feature, colname, mean_expr) %>%
      pivot_wider(names_from = colname, values_from = mean_expr) %>%
      column_to_rownames("feature") %>%
      as.matrix()

    hmap <- Heatmap(t(expr_mat),
                    name = "norm_expr",
                    border = T,
                    show_column_dend = F,
                    show_row_names = F,
                    show_row_dend = F,
                    row_names_gp = gpar(fontsize = 7),
                    column_names_gp = gpar(fontsize = 7),
                    col = col_fun)

  })

  # Annotation heatmap

  col_list <- list(heart_failure = c("HF" = "black",
                                     "NF" = "darkgrey"),
                   study = col_pub$study_colors_sc)

  row_anns <- unique(ex_genes$colname) %>%
    tibble::enframe(value = "colname") %>%
    dplyr::select_at("colname") %>%
    dplyr::left_join(ex_genes %>%
                       ungroup %>%
                       dplyr::select(colname, heart_failure, study) %>%
                       unique(), by = "colname") %>%
    dplyr::select(colname, heart_failure, study) %>%
    as.data.frame() %>%
    column_to_rownames("colname")

  row_cols <- Heatmap(row_anns,
                      name = "anns",
                      border = T,
                      show_column_dend = F,
                      show_row_dend = F,
                      show_row_names = F,
                      row_names_gp = gpar(fontsize = 7),
                      column_names_gp = gpar(fontsize = 7),
                      col = c("HF" = "black",
                              "NF" = "darkgrey",
                              col_pub$study_colors_sc))

  pdf_out <- paste0("results/meta_mofacell/gene_expression_example_", Fact, ".pdf")

  pdf(pdf_out, height = 2.5, width = 4.6)

  draw(expr_mats[[1]] + expr_mats[[2]] + expr_mats[[3]] + expr_mats[[4]] +
         expr_mats[[5]] + expr_mats[[6]] + expr_mats[[7]] + row_cols,
       heatmap_legend_side = "bottom", annotation_legend_side = "bottom")

  dev.off()

}


