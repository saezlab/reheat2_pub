# Copyright (c) [2023] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Check if the multicellular signature relates to
#' development
#'

library(MOFAcellulaR)
library(tidyverse)
library(cowplot)
library(biomaRt)
source("./code/reheat2_pilot/aesthetics.R")

# outs
factor_stats_time <- "./results/predict_project/TAC/TAC_assocs_time.csv"
factor_stats_hf <- "./results/predict_project/TAC/TAC_assocs.csv"
factor_bplots <- "./results/predict_project/TAC/TAC_factors.pdf"

#' Basic function to convert human to mouse gene names
convertMouseGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

  genesV2 = getLDS(attributes = c("mgi_symbol"),
                   filters = "mgi_symbol",
                   values = x , mart = mouse,
                   attributesL = c("hgnc_symbol"),
                   martL = human, uniqueRows=T)

  humanx <- unique(genesV2[, 2])

  # Print the first 6 genes found to the screen
  return(genesV2)
}

# New data
dev_data <- readRDS("./data/val_studies/mice/GSE1200064_pb.rds")

# Counts
pb_data <- dev_data$pb
mice_genes <- convertMouseGeneList(x = rownames(pb_data))

# Duplicated mice and human genes, are deleted
repeated_m <- mice_genes$MGI.symbol[which(duplicated(mice_genes$MGI.symbol))] %>% unique()
repeated_h <- mice_genes$MGI.symbol[which(duplicated(mice_genes$HGNC.symbol))] %>% unique()
out_genes <- unique(repeated_m, repeated_h)

# useful genes
useful_genes <- rownames(pb_data)[!(rownames(pb_data) %in% out_genes)]

mice_genes <- mice_genes %>%
  dplyr::filter(MGI.symbol %in% useful_genes) %>%
  dplyr::filter(MGI.symbol %in% rownames(pb_data)) %>%
  column_to_rownames("MGI.symbol")

mice_genes <- na.omit(mice_genes[useful_genes, , drop = F])

mice_genes <- mice_genes[rownames(mice_genes) %in% rownames(pb_data), , drop =F]

# Finally filter and rename
pb_data <- pb_data[rownames(mice_genes),]
rownames(pb_data) <- mice_genes$HGNC.symbol

# Col-data
coldat <- dev_data$target  %>%
  column_to_rownames("pb_id") %>%
  dplyr::rename(ncells = "n") %>%
  dplyr::filter(celltype != "none")

pb_data <- pb_data[,rownames(coldat)]

# Processing the test data
# Defining cts
cts <- coldat$celltype %>%
  unique() %>%
  set_names()

# Creating MOFAcell object to make projection
MOFAcell_obj <- MOFAcellulaR::create_init_exp(pb_data, coldat) %>%
  MOFAcellulaR::filt_profiles(pb_dat = .,
                              cts = cts,
                              ncells = 20,
                              counts_col = "ncells",
                              ct_col = "celltype") %>%
  MOFAcellulaR::tmm_trns(pb_dat_list = .,
                         scale_factor = 1000000) %>%
  MOFAcellulaR::center_views(pb_dat_list = .) %>%
  MOFAcellulaR::pb_dat2MOFA(pb_dat_list = .,
                            sample_column = "sample")


# Meta data
meta_data <- coldat[, c("sample", "condition")] %>%
  unique() %>%
  as_tibble() %>%
  dplyr::mutate(condition = factor(condition,
                                   levels = c("0w", "2w", "5w", "8w", "11w"))) %>%
  dplyr::mutate(time = as.character(condition) %>%
                  gsub("w","",.) %>%
                  as.numeric(),
                simpf_condition = ifelse(condition == "0w", "NF", "HF"))

# Projection in all datasets
reheat_input = tibble(study = "reheat2",
                      mofa_out = "./results/meta_mofacell/metamodel_mofa.hdf5")


# Projection in meta model:

model <- MOFA2::load_model("./results/meta_mofacell/metamodel_mofa.hdf5")

test_projection <- project_data(model = model,
                                test_data = MOFAcell_obj)

time_assocs <- test_projection %>%
  as.data.frame(.) %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "Factor", values_to = "Score") %>%
  left_join(meta_data) %>%
  dplyr::filter(Factor %in% c("Factor1", "Factor2")) %>%
  group_by(Factor) %>%
  nest() %>%
  mutate(pair_ttests = map(data, function(dat) {
    lm(Score ~ 1 + time, dat) %>%
      broom::tidy()
  })) %>%
  dplyr::select(Factor, pair_ttests) %>%
  unnest() %>%
  dplyr::filter(term == "time")

hf_assocs <- test_projection %>%
  as.data.frame(.) %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "Factor", values_to = "Score") %>%
  left_join(meta_data) %>%
  dplyr::filter(Factor %in% c("Factor1", "Factor2")) %>%
  group_by(Factor) %>%
  nest() %>%
  mutate(pair_ttests = map(data, function(dat) {
    t.test(Score ~ simpf_condition, dat) %>%
      broom::tidy()
  })) %>%
  dplyr::select(Factor, pair_ttests) %>%
  unnest() %>%
  ungroup() %>%
  dplyr::mutate(p_adj = p.adjust(p.value, "BH"))

write_csv(time_assocs, factor_stats_time)
write_csv(hf_assocs, factor_stats_hf)

# Plot bplots:

f_interest <- c("Factor1", "Factor2")

plots_hf <- map(set_names(f_interest), function(fact) {

  useful_factors <- test_projection %>%
    as.data.frame(.) %>%
    rownames_to_column("sample") %>%
    pivot_longer(-sample, names_to = "Factor", values_to = "Score") %>%
    left_join(meta_data) %>%
    dplyr::filter(Factor %in% c(fact)) %>%
    ggplot(aes(x = simpf_condition, y = Score, color = simpf_condition)) +
    geom_boxplot() +
    geom_point(size = 3, alpha = 0.8) +
    theme_bw() +
    theme(axis.text = element_text(size =12),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_color_manual(values = c("black", "darkgrey")) +
    ylab(fact)

})

plots_time <- map(set_names(f_interest), function(fact) {

  useful_factors <- test_projection %>%
    as.data.frame(.) %>%
    rownames_to_column("sample") %>%
    pivot_longer(-sample, names_to = "Factor", values_to = "Score") %>%
    left_join(meta_data) %>%
    dplyr::filter(Factor %in% c(fact)) %>%
    ggplot(aes(x = condition, y = Score)) +
    geom_boxplot() +
    geom_point(size = 3, alpha = 0.8) +
    theme_bw() +
    theme(axis.text = element_text(size =12),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ylab(fact)

})


pdf(factor_bplots, width = 2.1, height = 2.5)
plot(plots_hf[[1]])
plot(plots_hf[[2]])
plot(plots_time[[1]]) + ylab("MCP1")
plot(plots_time[[2]]) + ylab("MCP2")
dev.off()

