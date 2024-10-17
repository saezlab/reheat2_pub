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
library(cowplot)

out_pdf <- "./results/bulk_integration/comp_consistency.pdf"
out_csv <- "./results/bulk_integration/comp_consistency.csv"

# Read data
reheat_summary <- read_delim("./data/reheat/meta_analysis_summary.txt",
                             delim = "\t",show_col_types = F)

# Following the same filtering as in MOFAcell model
markers <- read_csv("./results/consensus_markers.csv", show_col_types = F) %>%
  dplyr::filter(adj_Fisher_p <= 0.0001, mean_LFC > 2) %>%
  dplyr::ungroup() %>%
  dplyr::select(-adj_Fisher_p) %>%
  pivot_wider(names_from =  cell_type, values_from = mean_LFC) %>%
  column_to_rownames("gene") %>%
  as.matrix()

cell_type_comps <- read_csv("./results/compositions/t_res.csv",show_col_types = F) %>%
  column_to_rownames("cell_type") %>%
  rowMeans()

cell_type_comps <- cell_type_comps[colnames(markers)]

# Here's where we do the multiplication
# meanLFC of marker expression * t-value of compositional change (size-effect)

comp_scores <- t(t(markers) * cell_type_comps)

comp_scores <- rowSums(comp_scores, na.rm = T) %>%
  enframe(name = "gene", value = "comp_score")

# Put them together

reheat_summary <- reheat_summary %>%
  left_join(comp_scores) %>%
  dplyr::mutate(comp_consistent = ifelse(sign(mean_t) == sign(comp_score),
                       "yes", "no")) %>%
  dplyr::mutate(comp_consistent = ifelse(is.na(comp_consistent),
                                         "unknown", comp_consistent))

# Make plots and stats

pdf(out_pdf, height = 2.5, width = 7.3, onefile = TRUE)

for(n_genes in c(100,250,500, 1000, nrow(reheat_summary))){

  sliced_data <- reheat_summary %>%
    dplyr::slice(1:n_genes)

  prop_plt <- sliced_data %>%
    group_by(comp_consistent) %>%
    summarise(n_genes = n()) %>%
    mutate(prop_genes = n_genes/sum(n_genes)) %>%
    ggplot(aes(y = comp_consistent, x = prop_genes)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    theme(axis.text = element_text(size =12)) +
    xlab("Proportion of genes") +
    ylab("") +
    ggtitle(paste0("top ", n_genes))

  point_plt <- sliced_data %>%
    dplyr::select(comp_score, mean_t, fisher_pvalue) %>%
    na.omit() %>%
    ggplot(aes(x = comp_score, y = mean_t, size = -log10(fisher_pvalue))) +
    theme_bw() +
    theme(axis.text = element_text(size =12)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0, color = "red") +
    geom_vline(xintercept = 0, color = "red") +
    ylab("Bulk t-value (HF - NF)") +
    xlab("Comp. score (HF-NF)")

  joint_plt <- cowplot::plot_grid(prop_plt, point_plt, nrow = 1, rel_widths = c(0.7,1))

  plot(joint_plt)

}

dev.off()

# Write final composition matrix

colnames(markers) <- paste0("mrkrLFC_", colnames(markers))
markers <- markers %>%
  as.data.frame() %>%
  rownames_to_column("gene")


reheat_summary <- reheat_summary %>%
  left_join(markers, by = "gene") %>%
  unique()

write_csv(reheat_summary, out_csv)


top500 <- reheat_summary %>%
  dplyr::slice(1:500)

cor.test(top500$mean_t, top500$comp_score, method = "spearman")

