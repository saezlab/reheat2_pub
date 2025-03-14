# Here we generate numbers of the publication
library(tidyverse)

meta_model_stats <- read_csv("./results/meta_mofacell/model_stats.csv")

meta_model_stats %>%
  dplyr::mutate(high = ifelse(view == "Fib", view, "rest")) %>%
  t.test(R2_HF_consensus ~ high, data = .)

meta_model_stats %>%
  dplyr::mutate(high = ifelse(view == "Fib", view, "rest")) %>%
  t.test(n_genes ~ high, data = ., alternative = "greater")

meta_model_stats %>%
  ggplot(aes(x = view, y = n_genes)) +
  geom_boxplot()


meta_model_stats %>%
  dplyr::group_by(view) %>%
  summarize(mean(R2_HF_consensus))


meta_model_stats %>%
  dplyr::filter(! view %in% c("Fib", "vSMCs")) %>%
  aov(R2_HF_consensus ~ view, data = .) %>%
  broom::tidy()
  summarize(mean(R2_HF_consensus))


####
mofamodel_file <- paste0("./results/meta_mofacell/metamodel_mofa.hdf5")
model <- MOFA2::load_model(mofamodel_file)

MOFAcellulaR::get_geneweights(model = model,factor = "Factor1")


## Mcell percentage
loadings <- read_csv("./results/meta_mofacell/gene_loadings.csv") %>%
  dplyr::filter(Factor %in% c("Factor1", "Factor2"))

loadings %>%
  dplyr::filter(abs(value) >= 0.1) %>%
  dplyr::select(Factor, feature) %>%
  group_by(Factor, feature) %>%
  summarize(n_cells = n()) %>%
  dplyr::mutate(total_genes = n()) %>%
  ungroup() %>%
  dplyr::mutate(n_cells = ifelse(n_cells > 1, "mcell", "ucell")) %>%
  dplyr::select(-feature) %>%
  dplyr::group_by(Factor, n_cells) %>%
  dplyr::mutate(prop_genes = n()) %>%
  unique() %>%
  dplyr::mutate(perc = prop_genes/total_genes)
