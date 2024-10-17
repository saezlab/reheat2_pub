# Copyright (c) [2023] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we calculate an overall annotation of ReHeaT
#'

library(tidyverse)

mol_csv <- "./results/bulk_integration/mol_consistency.csv"
comp_csv <- "./results/bulk_integration/comp_consistency.csv"
out_pdf <- "./results/bulk_integration/scell_consistency.pdf"
outdot_pdf <- "./results/bulk_integration/scell_consistency_dots.pdf"
outbar_pdf <- "./results/bulk_integration/scell_consistency_percs.pdf"


# Full reheat + mofacell loadings

mol_table <- read_csv(mol_csv)

# Compositional info

comp_table <- read_csv(comp_csv)

comp_cols <- colnames(comp_table)[which(grepl("mrkr",colnames(comp_table)))]

comp_table <- comp_table %>%
  select_at(c("gene", "comp_score", "comp_consistent",comp_cols)) %>%
  unique()

# Full info
reheat <- mol_table %>%
  left_join(comp_table, by = "gene") %>%
  group_by(gene) %>%
  dplyr::slice(1)

# Make manual annotations
reheat <- reheat %>%
  dplyr::mutate(scell_consistency = ifelse(mol_consistent == "yes" & comp_consistent == "yes", "comp and mol", "unknown")) %>%
  dplyr::mutate(scell_consistency = ifelse(mol_consistent == "yes" &
                                 (comp_consistent == "no" | comp_consistent == "unknown"),
                               "molecular",
                               scell_consistency)) %>%
  dplyr::mutate(scell_consistency = ifelse(comp_consistent == "yes" &
                                             (mol_consistent == "no" | mol_consistent == "unknown"),
                                           "compositional",
                                           scell_consistency)) %>%
  dplyr::arrange(fisher_pvalue) %>%
  ungroup()

# Make figure panel

point_plt_comp <- reheat %>%
  #dplyr::slice(1:500) %>%
  dplyr::filter(fisher_pvalue < 0.05) %>%
  dplyr::select(comp_score, mol_score, mean_t, fisher_pvalue) %>%
  pivot_longer(-c(mean_t, fisher_pvalue),
               names_to = "score_type",
               values_to = "score") %>%
  dplyr::mutate(score_type = ifelse(score_type == "mol_score",
                                    "Molecular",
                                    "Compositional")) %>%
  ggplot(aes(y = score, x = mean_t)) +
  theme_minimal(base_size = 15) +
  # Add a border around each facet
  theme(panel.border=element_rect(fill=NA, colour="grey40"),
        axis.text = element_text(size =12),
        axis.title = element_text(size =12)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, color = "red") +
  geom_vline(xintercept = 0, color = "red") +
  facet_wrap(. ~ score_type,ncol = 1, scales = "free_y") +
  xlab("Bulk t-value \n (HF - NF)") +
  ylab("Scell score (HF-NF)")

pdf(outdot_pdf, height = 4, width = 2.2)

plot(point_plt_comp)

dev.off()


# Male summary plots

pdf(out_pdf, height = 3, width = 5, onefile = TRUE)

map(set_names(c(100,250,500, 1000, 2500, 5000,nrow(reheat))), function(n_gene){

  summary <- reheat %>%
    dplyr::slice(1:n_gene) %>%
    group_by(scell_consistency) %>%
    summarize(n_genes = n()) %>%
    dplyr::mutate(perc = n_genes/sum(n_genes)) %>%
    dplyr::mutate(profiled_genes = n_gene)

}) %>%
  enframe() %>%
  unnest(value) %>%
  ggplot(aes(y = scell_consistency, x = perc)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  ggplot2::facet_wrap(. ~ profiled_genes) +
  ylab("") +
  xlab("Proportion of genes")

dev.off()

write_csv(reheat, "./results/resource/reheat2_stats_bulk.csv")

reheat %>%
  dplyr::filter(scell_consistency %in% c("compositional", "comp and mol")) %>%
  select_at(c("gene", comp_cols)) %>%
  pivot_longer(-gene) %>%
  group_by(name) %>%
  na.omit() %>%
  summarise(n())

# Paper stats

# Correlations

top500 <- reheat %>%
  dplyr::filter(fisher_pvalue < 0.05)

cor.test(top500$mean_t, top500$mol_score, method = "spearman")
cor.test(top500$mean_t, top500$comp_score, method = "spearman")

#
summary <- reheat %>%
  dplyr::filter(fisher_pvalue < 0.05) %>%
  group_by(scell_consistency) %>%
  summarize(n_genes = n()) %>%
  dplyr::mutate(perc = n_genes/sum(n_genes))


styling <- theme_void() +
  theme(
    axis.ticks.y.right = element_line(),  # Add ticks on the right
    axis.ticks.length.y.right = unit(1, 'mm'),  # Control the length of the right ticks
    axis.text.y.right = element_text(size = 13.5)  # Style the right y-axis text
  )

p1 <- ggplot(summary,
             aes(
               x = 0,
               y = perc,
               fill = scell_consistency
             )) +
  geom_col(position = 'stack', color = 'black') +
  scale_fill_manual(values = c("#FF8383", "#A899FF",
                               "#0073AE", "darkgrey")) +
  scale_y_continuous(labels = scales::label_percent(1),
                     sec.axis = dup_axis()) +  # Duplicate the y-axis on the right side
  styling


pdf(outbar_pdf, height = 4, width = 2.3)

plot(p1)

dev.off()
