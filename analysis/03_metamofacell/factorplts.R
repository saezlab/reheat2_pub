# Copyright (c) [2024] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Summarize the results of the meta model

library(tidyverse)
library(MOFAcellulaR)

mofa_out <- "./results/meta_mofacell/metamodel_mofa.hdf5"
source("./code/reheat2_pilot/aesthetics.R")

###
patient_map <- "./results/meta_mofacell/pat_map.pdf"
factor_dists <- "./results/meta_mofacell/dist_factors.pdf"
r2_ctypes <- "./results/meta_mofacell/celltype_R2.pdf"

# Meta model
model <- MOFA2::load_model(mofa_out)

# We get all meta_data
all_meta <- tibble(file = list.files("./data/metadata_ext//")) %>%
  dplyr::mutate(study = gsub("_metadata[.]csv","",file)) %>%
  dplyr::filter(study %in% c("Chaffin2022_DCM", "Koenig2022_DCM",
                             "LVReichart2022_DCM", "Simonson2023_ICM")) %>%
  dplyr::mutate(file = paste0("./data/metadata_ext/",file)) %>%
  dplyr::mutate(metadata = map(file, ~ read_csv(.x, show_col_types = FALSE))) %>%
  unnest(metadata) %>%
  dplyr::select(-file)

# Plot multicellular space

# Etiology
main_facts_plt_etiol <- MOFAcellulaR::get_tidy_factors(model = model,
                                                 factor = "all",
                                                 metadata = all_meta,
                                                 group = T,
                                                 sample_id_column = "sample_id") %>%
  dplyr::filter(Factor %in% c("Factor1", "Factor2")) %>%
  pivot_wider(values_from = value, names_from = Factor) %>%
  ggplot(aes(x = Factor1, y = Factor2, color = disease_code)) +
  geom_point() +
  theme_bw() +
  theme(axis.text = element_text(size =10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_manual(values = etiology_colors[all_meta$disease_code %>% unique]) +
  ylab("MCP2") + xlab("MCP1")

# Heart Failure
main_facts_plt_hf <- MOFAcellulaR::get_tidy_factors(model = model,
                                                       factor = "all",
                                                       metadata = all_meta,
                                                       group = T,
                                                       sample_id_column = "sample_id") %>%
  dplyr::filter(Factor %in% c("Factor1", "Factor2")) %>%
  pivot_wider(values_from = value, names_from = Factor) %>%
  ggplot(aes(x = Factor1, y = Factor2, color = heart_failure)) +
  geom_point() +
  theme_bw() +
  theme(axis.text = element_text(size =10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_manual(values = c("black","darkgrey")) +
  ylab("MCP2") + xlab("MCP1")

# Study
main_facts_plt_study <- MOFAcellulaR::get_tidy_factors(model = model,
                                                    factor = "all",
                                                    metadata = all_meta,
                                                    group = T,
                                                    sample_id_column = "sample_id") %>%
  dplyr::filter(Factor %in% c("Factor1", "Factor2")) %>%
  pivot_wider(values_from = value, names_from = Factor) %>%
  ggplot(aes(x = Factor1, y = Factor2, color = study)) +
  geom_point() +
  theme_bw() +
  theme(axis.text = element_text(size =10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_manual(values = study_colors[all_meta$study %>% unique]) +
  ylab("MCP2") + xlab("MCP1")


pmap_plot <- cowplot::plot_grid(main_facts_plt_etiol, main_facts_plt_hf, main_facts_plt_study, nrow = 3, ncol = 1, align = "hv")

pdf(patient_map, height = 5.5, width = 4.)

plot(pmap_plot)

dev.off()

####################################
# Distribution of scores

F1_plt <- MOFAcellulaR::get_tidy_factors(model = model,
                                         factor = "all",
                                         metadata = all_meta,
                                         group = T,
                                         sample_id_column = "sample_id") %>%
  dplyr::filter(Factor %in% c("Factor1")) %>%
  pivot_wider(values_from = value, names_from = Factor) %>%
  ggplot(aes(x = Factor1, y = disease_code, color = disease_code)) +
  geom_violin(width=1, alpha = 0.1, aes(fill = disease_code)) +
  geom_boxplot(width=0.15) +
  geom_point(size =1) +
  theme_bw() +
  scale_color_manual(values = etiology_colors[all_meta$disease_code %>% unique]) +
  scale_fill_manual(values = etiology_colors[all_meta$disease_code %>% unique]) +
  guides(fill = guide_legend(nrow = 2), color = guide_legend(nrow = 2)) +
  ylab("")


F2_plt <- MOFAcellulaR::get_tidy_factors(model = model,
                                         factor = "all",
                                         metadata = all_meta,
                                         group = T,
                                         sample_id_column = "sample_id") %>%
  dplyr::filter(Factor %in% c("Factor2")) %>%
  pivot_wider(values_from = value, names_from = Factor) %>%
  ggplot(aes(x = Factor2, y = disease_code, color = disease_code)) +
  geom_violin(width=1, alpha = 0.1, aes(fill = disease_code)) +
  geom_boxplot(width=0.15) +
  geom_point(size =1) +
  theme_bw() +
  scale_color_manual(values = etiology_colors[all_meta$disease_code %>% unique]) +
  scale_fill_manual(values = etiology_colors[all_meta$disease_code %>% unique]) +
  guides(fill = guide_legend(nrow = 2), color = guide_legend(nrow = 2)) +
  ylab("")


F5_plt <- MOFAcellulaR::get_tidy_factors(model = model,
                                         factor = "all",
                                         metadata = all_meta,
                                         group = T,
                                         sample_id_column = "sample_id") %>%
  dplyr::filter(Factor %in% c("Factor5")) %>%
  pivot_wider(values_from = value, names_from = Factor) %>%
  ggplot(aes(x = Factor5, y = disease_code, color = disease_code)) +
  geom_violin(width=1, alpha = 0.1, aes(fill = disease_code)) +
  geom_boxplot(width=0.15) +
  geom_point(size =1) +
  theme_bw() +
  scale_color_manual(values = etiology_colors[all_meta$disease_code %>% unique]) +
  scale_fill_manual(values = etiology_colors[all_meta$disease_code %>% unique]) +
  guides(fill = guide_legend(nrow = 2), color = guide_legend(nrow = 2)) +
  ylab("")

dist_plot <- cowplot::plot_grid(F1_plt, F2_plt, F5_plt, nrow = 3, ncol = 1, align = "hv")

pdf(factor_dists, height = 5.5, width = 4.4)

plot(dist_plot)

dev.off()

# Plot explained variance per cell-type

r2_ctype_plt <- model@cache$variance_explained$r2_per_factor %>%
  map(., ~ .x[c("Factor1","Factor2"), ] %>%
        as.data.frame() %>%
        rownames_to_column("Factor") %>%
        pivot_longer(-Factor,names_to = "cell_type")) %>%
  enframe() %>%
  unnest() %>%
  ggplot(aes(x = cell_type, y = value)) +
  geom_boxplot() +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5)) +
  facet_wrap(.~Factor)


pdf(r2_ctypes, height = 2.5, width = 4)

plot(r2_ctype_plt)

dev.off()

## Generating associations

# Association with HF (all studies)
# Association with HF
expl_var_HF <- MOFAcellulaR::get_associations(model = model,
                                              metadata = all_meta,
                                              sample_id_column = "sample_id",
                                              test_variable = "heart_failure",
                                              test_type = "categorical",
                                              categorical_type = "parametric",
                                              group = TRUE)


studies <- set_names(all_meta$study %>% unique)

assoc_list <- map(studies, function(s) {

  MOFAcellulaR::get_associations(model = model,
                                 metadata = all_meta %>%
                                   dplyr::filter(study == s),
                                 sample_id_column = "sample_id",
                                 test_variable = "heart_failure",
                                 test_type = "categorical",
                                 categorical_type = "parametric",
                                 group = TRUE)


})

assoc_list[["All"]] <- expl_var_HF

assoc_list <- enframe(assoc_list) %>%
  unnest()

pval_hmap <- assoc_list %>%
  dplyr::mutate(name = strsplit(name,"_") %>%
                  map_chr(., ~.x[[1]])) %>%
  dplyr::mutate(Factor = factor(Factor,
                                   levels = paste0("Factor", seq(1:10)))) %>%
  ggplot(aes(x = Factor, y = name, fill = -log10(adj_pvalue))) +
  geom_tile(color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5),
        legend.position = "bottom") +
  scale_fill_gradient2(low = "white", mid = "white", high = "#d62728",
                       midpoint = -log10(0.15)) +
  coord_equal() +
  xlab("") + ylab("")


pdf("./results/meta_mofacell/assocs_hmap.pdf", height = 2.5, width = 3.5)
plot(pval_hmap)
dev.off()

# Make heatmaps of R2

R2_list <- model@cache$variance_explained$r2_per_factor
col_fun <- colorRamp2(c(0, 30),
                      c("white", "#d62728"))


expr_mats <- map(R2_list, function(r2_mat) {

  hmap <- Heatmap(r2_mat,
                  name = "R2",
                  border = T,
                  show_column_dend = F,
                  show_row_names = T,
                  show_row_dend = F,
                  cluster_rows = F,
                  cluster_columns = F,
                  row_names_gp = gpar(fontsize = 7),
                  column_names_gp = gpar(fontsize = 7),
                  col = col_fun)
})

pdf("./results/meta_mofacell/celltype_R2_hmap.pdf", height = 1.5, width = 3)
draw(expr_mats[[1]] + expr_mats[[2]])
draw(expr_mats[[3]] + expr_mats[[4]])
dev.off()


# Source data

MOFAcellulaR::get_tidy_factors(model = model,
                               factor = "all",
                               metadata = all_meta,
                               group = T,
                               sample_id_column = "sample_id") %>%
  dplyr::filter(Factor %in% c("Factor1", "Factor2")) %>%
  write_csv("./Revision/figures/Figure3/Figure3B.csv")
