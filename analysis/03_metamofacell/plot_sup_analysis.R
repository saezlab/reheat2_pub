# Copyright (c) [2025] [Ricardo O. Ramirez Flores]
# flores@ebi.ac.uk

#' Visualization of results
#'
#'
library(tidyverse)
source("./code/reheat2_pilot/aesthetics.R")

# Sup. Fig 1
# Number of interacting genes with age and sex and comorbidities

agesex_res <- read_csv("./results/meta_mofacell/interacting_sexage.csv")

n_agesex <- agesex_res %>%
  group_by(name, test) %>%
  summarise(n_genes = sum(padj <= 0.05))

age_plt <- ggplot(n_agesex, aes(x = name, fill = test, y = n_genes)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(size =10),
        legend.title = element_blank()) +
  ylab("Diff. expressed \n genes") +
  xlab("") +
  scale_fill_manual(values = c("darkgreen", "skyblue"))

comorb_res <- read_csv("./results/meta_mofacell/interacting_etiology.csv")

n_comorb <- comorb_res %>%
  dplyr::filter(term %in% c("ICM", "HCM")) %>%
  group_by(name, term) %>%
  summarize(n_degs = sum(padj <= 0.05))

com_plt <- ggplot(n_comorb, aes(x = name, fill = term, y = n_degs)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(size =10),
        legend.title = element_blank()) +
  ylab("Diff. expressed \n genes") +
  xlab("") +
  scale_fill_manual(values = etiology_colors[unique(n_comorb$term)])

cowplot::plot_grid(age_plt, com_plt, nrow = 1, align = "hv")

# Sup. Fig 2

lmer_res <- read_csv("./results/meta_mofacell/HF_lmer_covars.csv")

lmer_res_plt <- lmer_res %>%
  dplyr::mutate(Factor = ifelse(Factor == "Factor1", "MCP1", "MCP2")) %>%
  ggplot(aes(y = Factor, x = term, fill = adj_p)) +
  geom_tile(color = "black") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        axis.text = element_text(size =10)) +
  xlab("") +
  ylab("") + coord_equal()


pdf("./results/meta_mofacell/covars_etio_plots.pdf", height = 3, width = 9)
cowplot::plot_grid(age_plt, lmer_res_plt, com_plt,
                   nrow = 1 ,align = "hv", axis = "tblr")
dev.off()


n_comorb %>%
  dplyr::group_by(term) %>%
  summarize(mean(n_degs))

n_agesex %>%
  dplyr::group_by(test) %>%
  summarize(mean(n_genes))

