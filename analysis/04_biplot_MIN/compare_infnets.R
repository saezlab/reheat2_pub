# Copyright (c) [2024] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Comparison of MCPs info nets

library(tidyverse)

MCP1_net <- read_csv("./results/MCP_characs/MIN_Factor1.csv") %>%
  dplyr::select(-c(R2,Estimate)) %>%
  pivot_wider(names_from = condition,
              values_from = cor_estimate) %>%
  dplyr::mutate(diff = abs(HF - NF))

MCP1_difftest <- MCP1_net %>%
  group_by(target) %>%
  nest() %>%
  dplyr::mutate(dif_test = map(data, function(dat) {

    t.test(dat$diff,mu = 0,alternative = "greater") %>%
      broom::tidy()

  })) %>%
  dplyr::select(target, dif_test) %>%
  unnest() %>%
  ungroup() %>%
  dplyr::mutate(padj = p.adjust(p.value, "BH")) %>%
  dplyr::arrange(padj, desc(estimate))

MCP1_net %>%
  arrange(desc(abs(diff))) %>%
  dplyr::mutate(label = ifelse(target == "CM", predictor, "")) %>%
ggplot(., aes(x = NF, y = HF, color = target, label = label)) +
  geom_point() +
  geom_text() +
  geom_abline(slope = 1, intercept = 0) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  xlim(-0.7, 0.7) +
  ylim(-0.7, 0.7)

cor.test(MCP1_net$NF, MCP1_net$HF, method = "spearman")

MCP1_cortest <- MCP1_net %>%
  group_by(target) %>%
  nest() %>%
  dplyr::mutate(dif_test = map(data, function(dat) {

    t.test(dat$diff,mu = 0,alternative = "greater") %>%
      broom::tidy()

  })) %>%
  dplyr::select(target, dif_test) %>%
  unnest() %>%
  ungroup() %>%
  dplyr::mutate(padj = p.adjust(p.value, "BH")) %>%
  dplyr::arrange(padj, desc(estimate))

t.test(MCP1_net$HF,MCP1_net$NF,paired = T, "greater")

MCP1_net %>%
  dplyr::filter(target == "CM") %>%
  arrange(desc(diff))



MCP1_overall <- MCP1_net %>%
  group_by(predictor, target) %>%
  dplyr::summarise(MCP1 = mean(cor_estimate)) %>%
  arrange(-abs(MCP1))

MCP2_net <- read_csv("./results/MCP_characs/MIN_Factor2.csv") %>%
  dplyr::select(-c(R2,Estimate)) %>%
  pivot_wider(names_from = condition,
              values_from = cor_estimate) %>%
  dplyr::mutate(diff = abs(HF - NF))

MCP2_difftest <- MCP2_net %>%
  group_by(target) %>%
  nest() %>%
  dplyr::mutate(dif_test = map(data, function(dat) {

    t.test(dat$diff,mu = 0,alternative = "greater") %>%
      broom::tidy()

  })) %>%
  dplyr::select(target, dif_test) %>%
  unnest() %>%
  ungroup() %>%
  dplyr::mutate(padj = p.adjust(p.value, "BH")) %>%
  dplyr::arrange(padj, desc(estimate))

cor.test(MCP2_net$NF, MCP2_net$HF)

MCP2_net %>%
  dplyr::filter(target == "Fib") %>%
  arrange(desc(diff))

