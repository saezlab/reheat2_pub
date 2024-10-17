# Copyright (c) [2023] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we evaluate how much a gene associated with HF
#' can be linked to changes in cell compositions
#'
#' The assumption is simpler, the more a gene is a marker
#' from a cell-type, if there's a consistency between the
#' cellular composition change and the log-fold change
#' then the gene has a compositional dependency

library(tidyverse)
library(MOFAcellulaR)
library(cowplot)

out_pdf <- "./results/bulk_integration/mol_consistency.pdf"
out_csv <- "./results/bulk_integration/mol_consistency.csv"

# Read data
reheat_summary <- read_delim("./data/reheat/meta_analysis_summary.txt",
                             delim = "\t",show_col_types = F) %>%
  unique()



model <- MOFA2::load_model("./results/meta_mofacell/metamodel_mofa.hdf5")

meta_data <- read_csv("./results/meta_mofacell/metamodel_meta.csv",
                      show_col_types = F)

# Get associations with HF
expl_var_HF <- MOFAcellulaR::get_associations(model = model,
                                              metadata = meta_data,
                                              sample_id_column = "sample_id",
                                              test_variable = "heart_failure",
                                              test_type = "categorical",
                                              categorical_type = "parametric",
                                              group = TRUE)


# Let's start with the most associated Factor

# First get direction HF - NF (as compositions, and original ReHeaT)
hf_weight <- MOFAcellulaR::get_tidy_factors(model, meta_data, "Factor1", TRUE, sample_id_column = "sample_id") %>%
  group_by(heart_failure) %>%
  dplyr::summarise(mean_score = mean(value)) %>%
  pivot_wider(names_from = heart_failure, values_from = mean_score) %>%
  dplyr::mutate(hf_sign = sign(HF-NF)) %>%
  pull(hf_sign)

# Get a general molecular score (mean weights)

factor_loadings <- MOFAcellulaR::get_geneweights(model, "Factor1")

molecular_scores <- factor_loadings  %>%
  dplyr::mutate(corr_value = value * hf_weight) %>% # to the same direction of ReHeaT
  dplyr::select(feature, corr_value, ctype) %>%
  dplyr::rename("gene" = "feature") %>%
  group_by(gene) %>%
  summarise(mol_score = mean(corr_value))

reheat_summary <- reheat_summary %>%
  left_join(molecular_scores) %>%
  dplyr::mutate(mol_consistent = ifelse(sign(mean_t) == sign(mol_score),
                                         "yes", "no")) %>%
  dplyr::mutate(mol_consistent = ifelse(is.na(mol_consistent),
                                         "unknown", mol_consistent))
# Make plots and stats

pdf(out_pdf, height = 2.5, width = 7.3, onefile = TRUE)

for(n_genes in c(100,250,500, 1000, nrow(reheat_summary))){

  sliced_data <- reheat_summary %>%
    dplyr::slice(1:n_genes)

  prop_plt <- sliced_data %>%
    group_by(mol_consistent) %>%
    summarise(n_genes = n()) %>%
    mutate(prop_genes = n_genes/sum(n_genes)) %>%
    ggplot(aes(y = mol_consistent, x = prop_genes)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    theme(axis.text = element_text(size =12)) +
    xlab("Proportion of genes") +
    ylab("") +
    ggtitle(paste0("top ", n_genes))

  point_plt <- sliced_data %>%
    ggplot(aes(x = mol_score, y = mean_t, size = -log10(fisher_pvalue))) +
    theme_bw() +
    theme(axis.text = element_text(size =12)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0, color = "red") +
    geom_vline(xintercept = 0, color = "red") +
    ylab("Bulk t-value (HF - NF)") +
    xlab("Molecular score (HF-NF)")

  joint_plt <- cowplot::plot_grid(prop_plt, point_plt, nrow = 1, rel_widths = c(0.7,1))

  plot(joint_plt)

}

dev.off()

# Make final table

factor_loadings <- factor_loadings %>%
  pivot_wider(names_from = ctype, values_from = value) %>%
  dplyr::rename("gene" = feature) %>%
  column_to_rownames("gene") %>%
  as.matrix()

factor_loadings <- factor_loadings * hf_weight # to the same direction of ReHeaT

colnames(factor_loadings) <- paste0("mcell_HFscore_", colnames(factor_loadings))

factor_loadings <- factor_loadings %>%
  as.data.frame() %>%
  rownames_to_column("gene")

reheat_summary <- reheat_summary %>%
  left_join(factor_loadings, by = "gene") %>%
  unique()

write_csv(reheat_summary, out_csv)

top500 <- reheat_summary %>%
  dplyr::slice(1:500)

cor.test(top500$mean_t, top500$mol_score, method = "spearman")











