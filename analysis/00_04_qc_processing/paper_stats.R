# Copyright (c) [2024] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we calculate the general QC of distinct studies
#' for paper text


outqc_file <- "./results/qc/qc_stats.csv"
qc <- read_csv(outqc_file)


qc %>%
  dplyr::select(study_id, cell_type, psbulk_n_genes, psbulk_n_cells) %>%
  unique() %>%
  pull(psbulk_n_genes) %>%
  median()

qc %>%
  dplyr::select(study_id, cell_type, psbulk_n_genes, psbulk_n_cells) %>%
  unique() %>%
  pull(psbulk_n_cells) %>%
  median()

# Testing differences in number of genes

qc %>%
  dplyr::select(study_id, cell_type, psbulk_n_genes, psbulk_n_cells) %>%
  unique() %>%
  group_by(cell_type) %>%
  nest() %>%
  mutate(n_genes_comp = map(data, function(dat) {

    aov(psbulk_n_genes ~ study_id, dat) %>%
      broom::tidy() %>%
      dplyr::filter(term == "study_id")

  })) %>%
  dplyr::select(n_genes_comp) %>%
  unnest() %>%
  ungroup() %>%
  mutate(adj_p = p.adjust(p.value))

qc %>%
  dplyr::select(study_id, cell_type, psbulk_n_genes, psbulk_n_cells) %>%
  unique() %>%
  group_by(cell_type) %>%
  nest() %>%
  mutate(n_genes_comp = map(data, function(dat) {

    pairwise.t.test(dat$psbulk_n_genes,
                    dat$study_id,
                    p.adjust.method = "BH") %>%
      broom::tidy()
  })) %>%
  dplyr::select(n_genes_comp) %>%
  unnest() %>%
  ungroup()

# Testing in number of cells

qc %>%
  dplyr::select(study_id, cell_type, psbulk_n_genes, psbulk_n_cells) %>%
  unique() %>%
  group_by(cell_type) %>%
  nest() %>%
  mutate(n_cells_comp = map(data, function(dat) {

    aov(psbulk_n_cells ~ study_id, dat) %>%
      broom::tidy() %>%
      dplyr::filter(term == "study_id")

  })) %>%
  dplyr::select(n_cells_comp) %>%
  unnest() %>%
  ungroup() %>%
  mutate(adj_p = p.adjust(p.value))

qc %>%
  dplyr::select(study_id, cell_type, psbulk_n_genes, psbulk_n_cells) %>%
  unique() %>%
  group_by(cell_type) %>%
  nest() %>%
  mutate(n_cells_comp = map(data, function(dat) {

    pairwise.t.test(dat$psbulk_n_cells,
                    dat$study_id,
                    p.adjust.method = "BH") %>%
      broom::tidy()
  })) %>%
  dplyr::select(n_cells_comp) %>%
  unnest() %>%
  dplyr::filter(p.value < 0.05) %>%
  print(n = 50)


pairwise.t.test(qc$psbulk_n_cells, qc$cell_type, p.adjust.method = "BH") %>%
  broom::tidy()


# Contamination statistics:

qc %>%
  group_by(cell_type) %>%
  nest() %>%
  mutate(cont_sign = map(data, function(dat) {
    t.test(dat$contamination, mu = 0, alternative = "greater") %>%
      broom::tidy()
  })) %>%
  dplyr::select(cont_sign) %>%
  unnest() %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p.value))  %>%
  dplyr::filter(p_adj < 0.05) %>%
  dplyr::arrange(cell_type, study_id)

qc %>%
  group_by(cell_type, study_id) %>%
  nest() %>%
  mutate(cont_sign = map(data, function(dat) {
    t.test(dat$contamination, mu = 1, alternative = "greater") %>%
      broom::tidy()
  })) %>%
  dplyr::select(cont_sign) %>%
  unnest() %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p.value))  %>%
  dplyr::filter(p_adj < 0.05) %>%
  dplyr::arrange(cell_type)


pairwise.t.test(qc$contamination, qc$cell_type, p.adjust.method = "BH") %>%
  broom::tidy() %>%
  dplyr::filter(p.value < 0.05)


pairwise.t.test(qc$contamination, qc$cell_type, p.adjust.method = "BH") %>%
  broom::tidy() %>%
  dplyr::filter(p.value < 0.05)

qc %>%
  ggplot(aes(x = cell_type, y = log1p(contamination))) +
  geom_boxplot()


cont_cols <- colnames(qc)[grepl("cont_", colnames(qc))]

cont_fracs <- qc[,c("cell_type", cont_cols)] %>%
  pivot_longer(-cell_type) %>%
  na.omit()


pairwise.t.test(cont_fracs$value, cont_fracs$name, p.adjust.method = "BH") %>%
  broom::tidy() %>%
  dplyr::filter(p.value < 0.05)


qc[,cont_cols] %>% as.matrix() %>%
  colMeans(na.rm = T)


qc %>%
  mutate(bckground = ifelse(study_id %in% c("Simonson2023_ICM", "Chaffin2022_DCM"),
                            "yes", "no")) %>%
  group_by(cell_type) %>%
  nest() %>%
  mutate(comp_res = map(data, function(dat) {

    t.test(contamination ~ bckground, data= dat) %>%
      broom::tidy()

  })) %>%
  unnest(comp_res) %>%
  ungroup() %>%
  dplyr::mutate(adj_p = p.adjust(p.value)) %>%
  dplyr::select(cell_type, statistic, adj_p)

cellbender_effect_plt <- qc %>%
  mutate(bckground = ifelse(study_id %in% c("Simonson2023_ICM", "Chaffin2022_DCM"),
                            "yes", "no")) %>%
  #dplyr::filter(cell_type == "Endo" ) %>%
  ggplot(aes(x = bckground, y = log1p(contamination))) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(~ cell_type,scales = "free_y") +
  theme_bw() +
  xlab("background correction")


qc %>%
  dplyr::select(cell_type, contamination) %>%
  group_by(cell_type) %>%
  summarize(n = mean(contamination))

