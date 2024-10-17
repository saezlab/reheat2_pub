# Copyright (c) [2023] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we calculate the compositional similarities
#' across studies
library(tidyverse)
library(compositions)
library(cowplot)
library(ComplexHeatmap)
library(lmerTest)
library(cluster)
library(circlize)


setwd("/Users/ricardoramirez/Dropbox/PostDoc/Research/ReHeaT2/")
hmp_plt <- "./results/compositions/sample_comps.pdf"
si_plt_file <- "./results/compositions/si_scores.pdf"
meta_file <- "./results/compositions/comps_meta.csv"
comps_file <- "./results/compositions/sample_comps.csv"
tval_hmp_plt <- "./results/compositions/cntrst_comps.pdf"
dens_plt_file <- "./results/compositions/comps_dists.pdf"
lmer_file <- "./results/compositions/lmer_res.csv"
lmer_plt <- "./results/compositions/lmer_res.pdf"
t_file <- "./results/compositions/t_res.csv"
lmer_bars <- "./results/compositions/lmer_bars.pdf"

# Aesthetics - colors for plotting
source("./code/reheat2_pilot/aesthetics.R")

# This is combining all data to make

coldata <- tibble(file = list.files("./data/coldata/")) %>%
  dplyr::mutate(study = gsub("_coldata[.]csv","",file)) %>%
  dplyr::mutate(file = paste0("./data/coldata/",file)) %>%
  dplyr::mutate(coldata = map(file, ~ read_csv(.x, show_col_types = FALSE)[,-1])) %>%
  unnest(coldata) %>%
  dplyr::select(-file) %>%
  dplyr::group_by(study, sample_id) %>%
  dplyr::mutate(sample_cells = sum(psbulk_n_cells)) %>%
  ungroup() %>%
  dplyr::mutate(cell_type_prop = psbulk_n_cells/sample_cells) %>%
  dplyr::filter(study != "Reichart2022_DCM")

meta_data <- coldata %>%
  dplyr::select(study, sample_id, disease_code, heart_failure) %>%
  unique()

write_csv(meta_data, meta_file)

prop_data <- coldata %>%
  dplyr::select(sample_id, cell_type, cell_type_prop)

prop_data <- prop_data %>%
  pivot_wider(values_from = cell_type_prop, names_from = cell_type, values_fill = 0) %>%
  column_to_rownames("sample_id") %>%
  as.matrix()

write_csv(as.data.frame(prop_data) %>%
            rownames_to_column("sample_id"),
          comps_file)

# Here we can calculate the sillouhette score based on distances in the compositions

meta_info_ss <- meta_data %>%
  dplyr::mutate(hf_ix = .data[["heart_failure"]] %>%
                  as.factor() %>%
                  as.integer(),
                study_ix  = .data[["study"]] %>%
                  as.factor() %>%
                  as.integer())


patient_dists <- dist(prop_data[meta_info_ss$sample_id,],
                      method = "euclidean")

si_hf <- (silhouette(x = meta_info_ss$hf_ix, patient_dists)) %>%
  as.data.frame() %>%
  dplyr::select(cluster, sil_width) %>%
  left_join(unique(meta_info_ss[,c("heart_failure", "hf_ix")]),
            by = c("cluster" = "hf_ix"))  %>%
  dplyr::select(-cluster) %>%
  dplyr::mutate(covar = "Disease") %>%
  dplyr::rename("var_name" = heart_failure)

si_study <- (silhouette(x = meta_info_ss$study_ix, patient_dists)) %>%
  as.data.frame() %>%
  dplyr::select(cluster, sil_width) %>%
  left_join(unique(meta_info_ss[,c("study", "study_ix")]),
            by = c("cluster" = "study_ix")) %>%
  dplyr::select(-cluster) %>%
  dplyr::mutate(covar = "Study") %>%
  dplyr::rename("var_name" = study)

si <- bind_rows(si_hf, si_study)

# Plot

si_plt <- si %>%
  ggplot(aes(x = var_name, y = sil_width)) +
  geom_boxplot() +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10)) +
  ylab("silhouette score") +
  facet_grid(.~covar,scales = "free_x", space='free') +
  xlab("")

pdf(si_plt_file, height = 3.2, width = 2.6)

plot(si_plt)

dev.off()

# Finding which groups are influenced?

si %>%
  group_by(var_name) %>%
  nest() %>%
  mutate(zero_dif = map(data, function(dat){

    t.test(dat$sil_width, mu = 0, alternative = "greater") %>%
      broom::tidy()

  })) %>%
  dplyr::select(-data) %>%
  unnest() %>%
  ungroup() %>%
  dplyr::mutate(adj_pval = p.adjust(p.value)) %>%
  dplyr::filter(adj_pval <= 0.05)

# Numbers for the paper

si %>%
  group_by(var_name) %>%
  summarize(median(sil_width))


# Testing differences in HF
prop_data_clr <- prop_data %>%
  compositions::acomp() %>%
  compositions::clr()

# Compare the vector of comparisons
study_diff_stats <- prop_data_clr %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(-sample_id, names_to = "cell_type") %>%
  left_join(meta_data, by = "sample_id") %>%
  group_by(study, cell_type) %>%
  nest() %>%
  dplyr::mutate(HF_diff = map(data, function(dat) {

    t.test(value ~ heart_failure, data = dat) %>%
      broom::tidy()

  })) %>%
  dplyr::select(-data) %>%
  unnest(HF_diff) %>%
  dplyr::select(cell_type, study, statistic, p.value) %>%
  dplyr::mutate(study = strsplit(study,"_") %>%
                  map_chr(., ~.x[[1]]))

t_matrix <- study_diff_stats %>%
  dplyr::select(-p.value) %>%
  pivot_wider(names_from = study, values_from = statistic) %>%
  column_to_rownames("cell_type") %>%
  as.matrix()

write_csv(as.data.frame(t_matrix)  %>%
            rownames_to_column("cell_type"), t_file)

p_matrix <- study_diff_stats %>%
  dplyr::select(-statistic) %>%
  group_by(study) %>%
  dplyr::mutate(adj_p = p.adjust(p.value)) %>%
  dplyr::select(-p.value) %>%
  pivot_wider(names_from = study, values_from = adj_p) %>%
  column_to_rownames("cell_type") %>%
  as.matrix()

# Finally linear mixed model of clr
# Studies are random effects
study_diff_stats_lmer <- prop_data_clr %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(-sample_id, names_to = "cell_type") %>%
  left_join(meta_data, by = "sample_id") %>%
  group_by(cell_type) %>%
  nest() %>%
  dplyr::mutate(HF_diff = map(data, function(dat) {

    dat <- dat %>%
      mutate(heart_failure = factor(heart_failure,
                                    levels = c("NF","HF")))

    # Satterthwaite's degrees of freedom method
    # for significance of fixed effect
    model_res <- lmerTest::lmer(value ~ heart_failure + (1 | study),
                                data = dat) %>%
      summary()

    fixed_summ <- model_res$coefficients %>%
      as.data.frame() %>%
      rownames_to_column("term") %>%
      dplyr::filter(grepl(pattern = "heart_failure", term)) %>%
      dplyr::rename("p_val" = "Pr(>|t|)")

    random_summ <- model_res$varcor %>%
      as.data.frame() %>%
      dplyr::select(grp, vcov) %>%
      pivot_wider(names_from = grp, values_from = vcov) %>%
      dplyr::mutate(perc_studyvar = study/(study + Residual)) %>%
      pull(perc_studyvar)

    all_res <- fixed_summ %>%
      dplyr::mutate(perc_studyvar = random_summ)

    return(all_res)

  })) %>%
  dplyr::select(-data) %>%
  unnest(c(HF_diff)) %>%
  ungroup() %>%
  dplyr::mutate(adj_pval = p.adjust(p_val))

write_csv(study_diff_stats_lmer, lmer_file)

# Finally we make the complex heatmap

# barplot color function
# barplot color function

t_matrix = t(t_matrix)
p_matrix = t(p_matrix)

color_fun <- function(estimates, p_vals) {
  ifelse(
    p_vals < 0.05,
    ifelse(estimates > 0,
           col_pub$general_up_down_colors["up"],
           col_pub$general_up_down_colors["down"]
    ),  # Positive estimates are blue, negative red
    "grey"  # Use grey if p-value is not significant
  )
}
# heatmap color function
heatmap_color_fun <- colorRamp2(
  c(min(t_matrix), 0, max(t_matrix)),  # Define the range of the data
  c(col_pub$general_up_down_colors["down"], "white", col_pub$general_up_down_colors["up"])
)
# Create the bar plot annotation with color coding
row_ha <- rowAnnotation(
  estimate = anno_barplot(
    study_diff_stats_lmer$Estimate,
    gp = gpar(fill = color_fun(study_diff_stats_lmer$Estimate, study_diff_stats_lmer$p_val)),
    ylim = c(min(study_diff_stats_lmer$Estimate, 0), max(study_diff_stats_lmer$Estimate, 0))
  )
)
star_matrix <- t(ifelse(p_matrix < 0.05, "*", ""))

# Create the heatmap with the left annotation and stars
comp_changes_hmap <- ComplexHeatmap::Heatmap(
  t(t_matrix),
  name = "t-value\nHF-NF",
  rect_gp = gpar(col = "black", lwd = 1),
  left_annotation = row_ha,
  show_column_dend = F,
  show_row_dend = F,
  col = heatmap_color_fun,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(star_matrix[i, j], x, y, gp = gpar(fontsize = 10, col = "black"))
  }
)

pdf(tval_hmp_plt, height = 3, width = 3)

draw(comp_changes_hmap)

dev.off()

# Plot distributions
significant_cells <- study_diff_stats_lmer %>%
  dplyr::filter(adj_pval <= 0.05) %>%
  pull(cell_type)

dens_plt <- prop_data %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(-sample_id, names_to = "cell_type") %>%
  left_join(meta_data, by = "sample_id") %>%
  dplyr::filter(cell_type %in% significant_cells) %>%
  ggplot(aes(x = value, y = heart_failure, color = heart_failure)) +
  geom_violin(width=.8, alpha = 0.1, aes(fill = heart_failure)) +
  geom_boxplot(width=0.15) +
  facet_wrap(. ~ cell_type,
             ncol = 2,
             nrow = 4,
             scales = "free") +
  scale_fill_manual(values = c("black", "darkgrey")) +
  scale_color_manual(values = c("black", "darkgrey")) +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("") +
  xlab("proportions")

pdf(dens_plt_file, height = 4, width = 3.7)

plot(dens_plt)

dev.off()


### PLOT of lmer
study_var <- read_csv(lmer_file) %>%
  ggplot(aes(x = -log10(adj_pval),
             y = perc_studyvar,
             label = cell_type)) +
  geom_vline(xintercept = -log10(0.05)) +
  geom_point() +
  geom_label(size =3) +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12)) +
  ylab("% of var by study")

###
pdf(lmer_plt, height = 3, width = 3.3)
plot(study_var)
dev.off()
