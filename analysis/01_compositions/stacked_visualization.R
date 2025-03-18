# Copyright (c) [2024] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we associate each covariate with
#' MCPs in the consensus analysis

library(tidyverse)
library(cowplot)

source("./code/reheat2_pilot/aesthetics.R")

# First use composition matrix to define order of patients

sample_compositions <- read_csv("./results/compositions/sample_comps.csv") %>%
  column_to_rownames("sample_id") %>%
  as.matrix()

# CT order
ct_order <- hclust(dist(t(sample_compositions)))
ct_order <- colnames(sample_compositions)[ct_order$order]
# Patient order
pat_order <- hclust(dist((sample_compositions)))
pat_order <- rownames(sample_compositions)[pat_order$order]

# Make heatmap of compositions
sample_compositions <- read_csv("./results/compositions/sample_comps.csv") %>%
  pivot_longer(-sample_id, names_to = "cell_type", values_to = "props")  %>%
  dplyr::mutate(sample_id = factor(sample_id, levels = pat_order),
                cell_type = factor(cell_type, levels = ct_order))

comps_plot <- ggplot(sample_compositions, aes(fill=cell_type, y=props, x=sample_id)) +
  geom_bar(position="stack", stat="identity", width= 1.1,  colour = NA) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(0, 2, 0, 0), 'cm'),
        legend.position = "top") +
  xlab("") +
  ylab("compositions") +
  scale_fill_manual(values = col_pub$ctype_colors)

# Make meta-data of patients

all_meta <- read_csv("./results/meta_mofacell/metamodel_meta.csv") %>%
  mutate(sample_id = factor(sample_id, levels = pat_order)) %>%
  dplyr::select(sample_id, heart_failure, disease_code, study)

study_plt <-ggplot(all_meta, aes(fill=study, y=1, x=sample_id)) +
  geom_tile() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(0, 2, 0, 0), 'cm'),
        legend.position = "top") +
  xlab("") +
  ylab("") +
  scale_fill_manual(values = col_pub$study_colors_sc)

hf_plt <- ggplot(all_meta, aes(fill=heart_failure, y=1, x=sample_id)) +
  geom_tile() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(0, 2, 0, 0), 'cm'),
        legend.position = "top") +
  xlab("") +
  ylab("") +
  scale_fill_manual(values = c("black","darkgrey"))

stacked_plt <- cowplot::plot_grid(comps_plot, study_plt, hf_plt,
                                  rel_heights = c(1, 0.4,0.4),
                                  ncol = 1, align = "hv")

pdf("results/compositions/stacked_comps.pdf", height = 6, width = 6)
plot(stacked_plt)
dev.off()

# Source

read_csv("./results/compositions/sample_comps.csv") %>%
  left_join(all_meta, by = "sample_id") %>%
  write_csv("./Revision/figures/Figure2/Figure2D.csv")


