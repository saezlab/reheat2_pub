# Copyright (c) [2025] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Calculating differences on hypertrophic and fibrotic markers
library(tidyverse)

qPCR_dat <- read_csv("./data/ExpData/hypertroph_fibro_data.csv")

meta_data <- read.csv("./data/ExpData/hypertroph_metadata_table.csv")

# Get control data

ctrl_dat <- qPCR_dat %>%
  dplyr::filter(grepl("Ctrl", sample_id)) %>%
  dplyr::select(-sheet) %>%
  unique()

ctrl_meta_data <- meta_data %>%
  dplyr::filter(grepl("Ctrl", sample_id)) %>%
  dplyr::select(-sheet) %>%
  unique()

# Get PE data

PE_dat <- qPCR_dat %>%
  dplyr::filter(grepl("PE", sample_id)) %>%
  dplyr::select(-sheet) %>%
  unique()

PE_meta_data <- meta_data %>%
  dplyr::filter(grepl("PE", sample_id)) %>%
  dplyr::select(-sheet) %>%
  unique()

# Data without PE and ctrl

qPCR_dat <- qPCR_dat %>%
  dplyr::filter(!grepl("PE", sample_id)) %>%
  dplyr::filter(!grepl("Ctrl", sample_id)) %>%
  dplyr::select(-sheet) %>%
  unique()

qPCR_dat <- bind_rows(qPCR_dat, PE_dat)

meta_data <- meta_data %>%
  dplyr::filter(!grepl("PE", sample_id)) %>%
  dplyr::filter(!grepl("Ctrl", sample_id)) %>%
  dplyr::select(-sheet) %>%
  unique()

meta_data <- bind_rows(meta_data, PE_meta_data)

# We make long formats
qPCR_dat <- pivot_longer(qPCR_dat, -sample_id,
                         names_to = "gene",
                         values_to = "expression") %>%
  left_join(meta_data)

ctrl_dat <- pivot_longer(ctrl_dat, -sample_id,
                         names_to = "gene",
                         values_to = "expression") %>%
  left_join(ctrl_meta_data)


test_dat = qPCR_dat %>%
  group_by(gene, condition) %>%
  nest() %>%
  dplyr::mutate(test_vs_ctrl = map2(gene, data, function(g, dat){

    ctrl_d <- ctrl_dat %>%
      dplyr::filter(gene == g) %>%
      dplyr::select(-c(gene))

    all_dat <- bind_rows(dat %>%
                           dplyr::mutate(condition = "rest"),
                         ctrl_d) %>%
      na.omit()

    disease_model <- suppressMessages(lmerTest::lmer("expression ~ condition + (1 | experiment)",
                                                     data = all_dat,verbose = 0))

    sum <- summary(disease_model)

    sum <- sum$coefficients %>%
      as.data.frame() %>%
      rownames_to_column("term") %>%
      dplyr::rename("pval" = "Pr(>|t|)",
                    "tval" = "t value") %>%
      dplyr::filter(term == "conditionrest") %>%
      dplyr::select(term, Estimate, pval, tval)

    return(sum)

})) %>%
  dplyr::select(gene, condition, test_vs_ctrl) %>%
  unnest()

# Visualization

test_dat <- test_dat %>%
  ungroup() %>%
  dplyr::mutate(adj_p = p.adjust(pval)) %>%
  dplyr::arrange(adj_p) %>%
  dplyr::mutate(sign = ifelse(adj_p <= 0.05, "*",""))

stats_plt <- test_dat %>%
  dplyr::mutate(condition = factor(condition,
                                   levels = c("Ctrl", "PE",
                                              "BMP4", "MXRA5", "NRG1"))) %>%
  dplyr::mutate(gene = factor(gene,
                              levels = c("Nppa", "Nppb",
                                         "Col1a1", "Tgfb1"))) %>%
  ggplot(aes(x = condition, y = gene, fill = tval, label = sign)) +
  geom_tile() +
  geom_text(color = "#FFF1DB", size = 7) +
  theme_minimal() +
  coord_equal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,
                                   vjust = 0.5)) +
  scale_fill_gradient2(low = "#EF5A6F", high = "#536493")

expr_plt <- qPCR_dat %>%
  bind_rows(ctrl_dat) %>%
  dplyr::mutate(condition = factor(condition,
                                   levels = c("Ctrl", "PE",
                                              "BMP4", "MXRA5", "NRG1"))) %>%
  dplyr::mutate(gene = factor(gene,
                                   levels = c("Nppa", "Nppb",
                                              "Col1a1", "Tgfb1"))) %>%
  na.omit() %>%
  ggplot(aes(x = condition, y = expression, fill = condition)) +
  geom_boxplot() +
  facet_wrap(.~ gene, scales = "free") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,
                                   vjust = 0.5),
        legend.position = "none")

pdf("./results/qPCR/hyper_fib_expression.pdf", height = 3.5,width = 3.5)

expr_plt

dev.off()

pdf("./results/qPCR/hyper_fib_stats.pdf", height = 2.5,width = 2.5)

stats_plt

dev.off()

