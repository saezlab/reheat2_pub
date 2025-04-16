# Copyright (c) [2025] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Calculating differences on response markers
library(tidyverse)

qPCR_dat <- read_csv("./data/ExpData/Mrkrs_Expression_Data.csv")

BMP4_t <- tibble(Gene = c("Gadd45b", "Sacs", "Ccnd1",
            "Ankrd1", "Ier3", "Prkd1",
            "Nav2", "Pros1", "Prom1", "Kcnmb4"),
            Condition = "BMP4",
            flag = 'keep')

MXRA5_t <- tibble(Gene = c("Gadd45b", "Ccnd1", "Ankrd1", "Ier3",
             "Mcl1", "Socs2", "Ahr", "Fgf9", "Ddit3"),
             Condition = "MXRA5",
             flag = 'keep')

NRG1_t <- tibble(Gene = c("Sacs", "Ccnd1", "Ankrd1", "Ier3",
            "Prickle2", "Lyar", "Ifnar2", "Bmal2", "Phlpp2"),
            Condition = "NRG1",
            flag = 'keep')

all_targets <- bind_rows(BMP4_t, MXRA5_t, NRG1_t)

ctrl_dat <- qPCR_dat %>%
  dplyr::filter(Condition == "CT")

qPCR_dat <- qPCR_dat %>%
  dplyr::filter(Condition != "CT")

test_dat = qPCR_dat %>%
  group_by(Gene, Condition) %>%
  nest() %>%
  dplyr::mutate(test_vs_ctrl = map2(Gene, data, function(g, dat){

    ctrl_d <- ctrl_dat %>%
      dplyr::filter(Gene == g) %>%
      dplyr::select(-c(Gene))

    all_dat <- bind_rows(dat %>%
                           dplyr::mutate(Condition = "rest"),
                         ctrl_d) %>%
      na.omit()

    disease_model <- suppressMessages(lmerTest::lmer("Expression ~ Condition + (1 | Experiment)",
                                                     data = all_dat,verbose = 0))

    sum <- summary(disease_model)

    sum <- sum$coefficients %>%
      as.data.frame() %>%
      rownames_to_column("term") %>%
      dplyr::rename("pval" = "Pr(>|t|)",
                    "tval" = "t value") %>%
      dplyr::filter(term == "Conditionrest") %>%
      dplyr::select(term, Estimate, pval, tval)

    return(sum)

  })) %>%
  dplyr::select(Gene, Condition, test_vs_ctrl) %>%
  unnest()

# Keep all the genes that are tested

PE_res <- test_dat %>%
  dplyr::filter(Condition == "PE",
                Gene %in% all_targets$Gene)

Lig_res <- test_dat %>%
  left_join(all_targets) %>%
  dplyr::filter(flag == "keep") %>%
  arrange(Gene) %>%
  dplyr::select(-flag)

all_res <- bind_rows(PE_res, Lig_res)  %>%
  ungroup() %>%
  dplyr::mutate(adj_p = p.adjust(pval)) %>%
  dplyr::arrange(adj_p) %>%
  dplyr::mutate(sign = ifelse(adj_p <= 0.05, "*","")) %>%
  dplyr::mutate(Condition = factor(Condition,
                                   levels = c("PE", "BMP4", "MXRA5", "NRG1"))) %>%
  complete(Gene, Condition)

stats_plt <- all_res %>%
ggplot(aes(x = Gene, y = Condition, fill = tval, label = sign)) +
  geom_tile(col = "black") +
  geom_text(color = "#FFF1DB", size = 7) +
  theme_bw() +
  coord_equal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,
                                   vjust = 0.5)) +
  scale_fill_gradient2(low = "#EF5A6F", high = "#536493", na.value = "grey") +
  ylab("condition") +
  xlab("gene")

pdf("./results/qPCR/mrkr_stats.pdf", height = 3.5,width = 5.0)

stats_plt

dev.off()








