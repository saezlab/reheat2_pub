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
factor_stats <- "./results/predict_project/angii/ANGII_assocs.csv"
factor_bplots <- "./results/predict_project/angii/ANGII_factors.pdf"

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
mice_data <- readRDS("./data/val_studies/mice/AngII_pb_sample_celltype.rds")

# Counts
pb_data <- mice_data$pb
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
coldat <- mice_data$target  %>%
  column_to_rownames("pb_id") %>%
  dplyr::rename(ncells = "cell_count",
                sample = "orig.ident") %>%
  dplyr::filter(!is.na(cell_type))

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
                            sample_column = "sample")


# Meta data
meta_data <- coldat[, c("sample", "group")] %>%
  unique() %>%
  as_tibble()

# Projection in meta model:

model <- MOFA2::load_model("./results/meta_mofacell/metamodel_mofa.hdf5")

test_projection <- project_data(model = model,
                                test_data = MOFAcell_obj)

test_assocs_t <- test_projection %>%
  as.data.frame(.) %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "Factor", values_to = "Score") %>%
  left_join(meta_data) %>%
  dplyr::filter(Factor %in% c("Factor1", "Factor2")) %>%
  group_by(Factor) %>%
  nest() %>%
  mutate(pair_ttests = map(data, function(dat) {
    t.test(Score ~ group, dat) %>%
      broom::tidy()
  })) %>%
  dplyr::select(Factor, pair_ttests) %>%
  unnest()

write_csv(test_assocs_t, factor_stats)

# Plot bplots:

f_interest <- c("Factor1", "Factor2")

plots <- map(set_names(f_interest), function(fact) {

  useful_factors <- test_projection %>%
    as.data.frame(.) %>%
    rownames_to_column("sample") %>%
    pivot_longer(-sample, names_to = "Factor", values_to = "Score") %>%
    left_join(meta_data) %>%
    dplyr::filter(Factor %in% c(fact)) %>%
    ggplot(aes(x = group, y = Score, color = group)) +
    geom_boxplot() +
    geom_point(size = 3, alpha = 0.8) +
    theme_bw() +
    theme(axis.text = element_text(size =12),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_color_manual(values = c("black", "darkgrey")) +
    ylab(fact)

})

pdf(factor_bplots, width = 3.0, height = 2.5)
plot(plots[[1]]) + ylab("MCP1")
plot(plots[[2]]) + ylab("MCP2")
dev.off()
