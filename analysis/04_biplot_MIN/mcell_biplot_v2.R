# Overlap tables to see unique and overlapping features

library(tidyverse)
library(cowplot)
library(ComplexHeatmap)
library(biclust)
source("./code/reheat2_pilot/aesthetics.R")

# This generates the coordinates of all enriched terms
all_single_cell_F1 <- read_csv("./results/meta_mofacell/Factor1_single_cell_enrichment_v2.csv") %>%
  dplyr::mutate(inf_score = GenesInList/GenesInPathway,
                pval_score = 1 - corr_p_value) %>%
  dplyr::mutate(coord = pval_score * inf_score) %>%
  dplyr::select(direction, gset, ct, coord) %>%
  pivot_wider(names_from = direction, values_from = coord) %>%
  dplyr::mutate(F1_coord = pos - neg) %>%
  dplyr::select(gset, ct, F1_coord)

all_single_cell_F2 <- read_csv("./results/meta_mofacell/Factor2_single_cell_enrichment_v2.csv") %>%
  dplyr::mutate(inf_score = GenesInList/GenesInPathway,
                pval_score = 1 - corr_p_value) %>%
  dplyr::mutate(coord = pval_score * inf_score) %>%
  dplyr::select(direction, gset, ct, coord) %>%
  pivot_wider(names_from = direction, values_from = coord) %>%
  dplyr::mutate(F2_coord = pos - neg) %>%
  dplyr::select(gset, ct, F2_coord)

all_single_cell <- left_join(all_single_cell_F1, all_single_cell_F2,
                             by = c("gset", "ct"))

# This collects significant terms
sign_single_cell_F1 <- read_csv("./results/meta_mofacell/Factor1_single_cell_enrichment_v2.csv") %>%
  dplyr::filter(corr_p_value <= 0.1) %>%
  dplyr::select(gset, ct)

sign_single_cell_F2 <- read_csv("./results/meta_mofacell/Factor2_single_cell_enrichment_v2.csv") %>%
  dplyr::filter(corr_p_value <= 0.1) %>%
  dplyr::select(gset, ct)

sign_single_cell <- bind_rows(sign_single_cell_F1, sign_single_cell_F2) %>%
  unique() %>%
  dplyr::mutate(class = "significant")

all_single_cell <- left_join(all_single_cell, sign_single_cell,
                             by = c("gset", "ct")) %>%
  dplyr::filter(class == "significant")

all_single_cell %>%
  dplyr::mutate(top_val = ifelse(abs(F1_coord) > abs(F2_coord),
                                 abs(F1_coord), abs(F2_coord))) %>%
  dplyr::filter(top_val >= 0.1) %>%
  arrange(ct, desc(top_val)) %>%
  dplyr::select(-top_val) %>%
  write_csv("./results/MCP_characs/singlecell_sets_biplot.csv")


# Select manually sets of interest
gsets_scell <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                 "GOBP_SARCOMERE_ORGANIZATION",
                 "HALLMARK_MYOGENESIS",
                 "GOMF_FATZ_BINDING",
                 "HALLMARK_FATTY_ACID_METABOLISM",
                 "WP_ANGIOGENESIS",
                 "NABA_COLLAGENS",
                 "GOBP_RESPONSE_TO_NITROSATIVE_STRESS",
                 "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                 "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                 "BIOCARTA_VEGF_PATHWAY",
                 "GOBP_BLOOD_VESSEL_ENDOTHELIAL_CELL_DIFFERENTIATION")

# Jans curation
singlecell_pways= c(##CM
  #"HP_ABNORMAL_SARCOMERE_MORPHOLOGY",
  "GOBP_SARCOMERE_ORGANIZATION",
  "HP_ABNORMAL_Z_DISC_MORPHOLOGY",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "KEGG_FATTY_ACID_METABOLISM",
  "KEGG_CARDIAC_MUSCLE_CONTRACTION",
  "KEGG_CALCIUM_SIGNALING_PATHWAY",
  ##FIB
  "NABA_COLLAGENS",
  "BIOCARTA_RECK_PATHWAY",
  ## ENDO
  "GOBP_BLOOD_VESSEL_ENDOTHELIAL_CELL_DIFFERENTIATION",
  "GOBP_LYMPHANGIOGENESIS",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "BIOCARTA_NO1_PATHWAY",
  "WP_ANGIOGENESIS",
  ##PC
  "GOBP_RESPONSE_TO_NITROSATIVE_STRESS",
  ##Myeloid
  "HALLMARK_IL6_JAK_STAT3_SIGNALING"
)

gsets_scell <- union(gsets_scell, singlecell_pways)

gsets_scell <- singlecell_pways

temp <- all_single_cell %>%
  dplyr::mutate(top_val = ifelse(abs(F1_coord) > abs(F2_coord),
                                 abs(F1_coord), abs(F2_coord))) %>%
  dplyr::filter(gset %in% gsets_scell)


xmax <- abs(temp$F1_coord) %>% max()
ymax <- abs(temp$F2_coord) %>% max()

single_cell_plt <- all_single_cell %>%
  dplyr::filter(gset %in% gsets_scell) %>%
  dplyr::mutate(top_val = ifelse(abs(F1_coord) > abs(F2_coord),
                                 abs(F1_coord), abs(F2_coord))) %>%
  dplyr::filter(top_val >= 0.1) %>%
  dplyr::mutate(gset = sub("^[^_]*_", "", gset)) %>%
  dplyr::mutate(gset = gsub("_", " ", gset) %>%
                  tolower() %>%
                  stringr::str_to_title()) %>%
  dplyr::arrange(ct, gset) %>%
  dplyr::group_by(ct, gset) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(ct, top_val) %>%
  dplyr::group_by(ct) %>%
  #dplyr::slice(1:5) %>% # Max of 5 per cell-types
  ggplot(., aes(x = 0, y = 0)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_segment(aes(xend = F1_coord, yend = F2_coord, color = ct),
               arrow = arrow(length = unit(0.07, "inches"))) +
  # ggrepel::geom_label_repel(aes(x = F1_coord,
  #                               y = F2_coord,
  #                               label = gset),
  #                           size = 3.2, force = 2,
  #                           color = NA,
  #                           fill = "white",
  #                           alpha = 0.7) +
  ggrepel::geom_text_repel(aes(x = F1_coord,
                               y = F2_coord,
                               label = gset),
                           color = "black") +
  theme_bw() +
  #xlim(-1 * (xmax - 0.2), (xmax - 0.2)) +
  ylim(-1 * (ymax), (ymax)) +
  labs(title = "Cell-type specific", x = "MCP1", y = "MCP2") +
  theme(axis.text = element_text(size =13),
        axis.title = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_manual(values = col_pub$ctype_colors)

pdf("./results/MCP_characs/ct_biplot.pdf", height = 5.3, width = 6.5)
plot(single_cell_plt)
dev.off()

# Now repeat for multicell
# Multicell
all_multi_cell_F1 <- read_csv("./results/meta_mofacell/Factor1_multicell_enrichment_v2.csv") %>%
  dplyr::select(-collection) %>%
  unique() %>%
  dplyr::arrange(direction, corr_p_value) %>%
  dplyr::group_by(direction, gset) %>%
  dplyr::slice(1) %>%
  dplyr::mutate(inf_score = GenesInList/GenesInPathway,
                pval_score = 1 - corr_p_value) %>%
  dplyr::mutate(coord = pval_score * inf_score) %>%
  dplyr::select(direction, gset, coord) %>%
  #unique() %>%
  pivot_wider(names_from = direction, values_from = coord) %>%
  dplyr::mutate(F1_coord = pos - neg) %>%
  dplyr::select(gset, F1_coord)

all_multi_cell_F2 <- read_csv("./results/meta_mofacell/Factor2_multicell_enrichment_v2.csv") %>%
  dplyr::select(-collection) %>%
  unique() %>%
  dplyr::arrange(direction, corr_p_value) %>%
  dplyr::group_by(direction, gset) %>%
  dplyr::slice(1) %>%
  dplyr::mutate(inf_score = GenesInList/GenesInPathway,
                pval_score = 1 - corr_p_value) %>%
  dplyr::mutate(coord = pval_score * inf_score) %>%
  dplyr::select(direction, gset, coord) %>%
  #unique() %>%
  pivot_wider(names_from = direction, values_from = coord) %>%
  dplyr::mutate(F2_coord = pos - neg) %>%
  dplyr::select(gset, F2_coord)

all_multi_cell <- left_join(all_multi_cell_F1, all_multi_cell_F2,
                            by = c("gset"))

# This collects significant terms
sign_multi_cell_F1 <- read_csv("./results/meta_mofacell/Factor1_multicell_enrichment_v2.csv") %>%
  dplyr::select(-collection) %>%
  unique() %>%
  dplyr::arrange(direction, corr_p_value) %>%
  dplyr::group_by(direction, gset) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  dplyr::select(gset, corr_p_value) %>%
  dplyr::filter(corr_p_value <= 0.1)

sign_multi_cell_F2 <- read_csv("./results/meta_mofacell/Factor2_multicell_enrichment_v2.csv") %>%
  dplyr::select(-collection) %>%
  unique() %>%
  dplyr::arrange(direction, corr_p_value) %>%
  dplyr::group_by(direction, gset) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  dplyr::select(gset, corr_p_value) %>%
  dplyr::filter(corr_p_value <= 0.1)

sign_multi_cell <- bind_rows(sign_multi_cell_F1, sign_multi_cell_F2) %>%
  ungroup() %>%
  dplyr::select(gset) %>%
  unique() %>%
  dplyr::mutate(class = "significant")

all_multi_cell <- left_join(all_multi_cell, sign_multi_cell,
                            by = c("gset")) %>%
  dplyr::filter(class == "significant")

temp <- all_multi_cell %>%
  dplyr::mutate(top_val = ifelse(abs(F1_coord) > abs(F2_coord),
                                 abs(F1_coord), abs(F2_coord))) %>%
  arrange(desc(top_val))

temp %>% write_csv("./results/MCP_characs/multicell_sets_biplot.csv")


gsets_mcell <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                 "HALLMARK_MYC_TARGETS_V1",
                 "HALLMARK_FATTY_ACID_METABOLISM",
                 "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
                 "HALLMARK_TGF_BETA_SIGNALING",
                 "REACTOME_INTRINSIC_PATHWAY_FOR_APOPTOSIS",
                 "BIOCARTA_STRESS_PATHWAY",
                 "KEGG_CARDIAC_MUSCLE_CONTRACTION",
                 "GOBP_MUSCLE_HYPERTROPHY",
                 "HALLMARK_HYPOXIA",
                 "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                 "HALLMARK_INTERFERON_ALPHA_RESPONSE")

multicell_pways <- c(## fact1 & fact2
  "BIOCARTA_VIP_PATHWAY",
  "BIOCARTA_HDAC_PATHWAY",
  ## fact1
  "BIOCARTA_CALCINEURIN_PATHWAY",
  "GOBP_POSITIVE_REGULATION_OF_NITRIC_OXIDE_METABOLIC_PROCESS",
  #"HP_VENTRICULAR_HYPERTROPHY",
  "HALLMARK_FATTY_ACID_METABOLISM",
  ## fact2
  "HALLMARK_HYPOXIA",
  "BIOCARTA_IGF1MTOR_PATHWAY",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)

gsets_mcell <- union(gsets_mcell, multicell_pways)

temp <- all_multi_cell %>%
  dplyr::mutate(top_val = ifelse(abs(F1_coord) > abs(F2_coord),
                                 abs(F1_coord), abs(F2_coord))) %>%
  dplyr::filter(gset %in% gsets_mcell)

xmax <- abs(temp$F1_coord) %>% max()
ymax <- abs(temp$F2_coord) %>% max()

single_multicell_plt <- all_multi_cell %>%
  dplyr::filter(gset %in% gsets_mcell) %>%
  dplyr::mutate(gset = sub("^[^_]*_", "", gset)) %>%
  dplyr::mutate(gset = gsub("_", " ", gset) %>%
                  tolower() %>%
                  stringr::str_to_title()) %>%
  ggplot(., aes(x = 0, y = 0)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_segment(aes(xend = F1_coord, yend = F2_coord),
               arrow = arrow(length = unit(0.07, "inches"))) +
  # ggrepel::geom_label_repel(aes(x = F1_coord,
  #                               y = F2_coord,
  #                               label = gset),
  #                           size = 3.2, force = 2,
  #                           color = NA,
  #                           fill = "white",
  #                           alpha = 0.7) +
  ggrepel::geom_text_repel(aes(x = F1_coord,
                               y = F2_coord,
                               label = gset),
                           color = "black") +
  theme_bw() +
  xlim(-1 * (xmax + 0.05), (xmax + 0.05)) +
  ylim(-1 * (ymax - 0.1), (ymax + 0.05)) +
  labs(title = "Multicellular", x = "MCP1", y = "MCP2") +
  theme(axis.text = element_text(size =13),
        axis.title = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pdf("./results/MCP_characs/multicell_biplot.pdf", height = 5.3, width = 6.5)
plot(single_multicell_plt)
dev.off()

all_bis <- cowplot::plot_grid(single_cell_plt, single_multicell_plt,
                   align = "hv",
                   ncol = 2)

pdf("./results/MCP_characs/all_biplots.pdf", height = 5, width = 13)
plot(all_bis)
dev.off()

