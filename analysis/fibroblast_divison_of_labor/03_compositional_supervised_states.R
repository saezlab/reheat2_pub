# HEADER --------------------------------------------
#
# Author:     Jan D. Lanzer
# Copyright   Copyright 2023 - Jan D. Lanzer
# Email:      lanzerjan@gmail.com
#
# Date:     2023-10-24
#
# Script Name:    ~/R-projects/reheat2_pilot/sub_cell_type/compositional_supervised_states.R
#
# Script Description:
#
#' In this script we calculate the compositional similarities
#' across studies
#' 
library(tidyverse)
library(compositions)
library(cowplot)
library(ComplexHeatmap)
library(lmerTest)

source("aesthetics.R")

hmp_plt <- "output/fib_sub_analysis/figures/sample_comps.pdf"
tval_hmp_plt <- "output/fib_sub_analysis/figures/sample_comps_t.pdf"
prop_data= read.csv("output/fib_sub_analysis/cluster_composition_sub.csv")
meta_data= read.csv("output/fib_sub_analysis/obs_df_sub.csv")

meta_data= meta_data%>% select(sample_id, disease_code, heart_failure, batch, disease_code)%>% distinct()
prop_data
prop_data= prop_data %>% pivot_longer(cols = c(paste0("X", 0:5)), 
                             names_to = "cell_type", 
                             values_to= "cell_type_prop")%>%
  mutate(cell_type = str_replace_all(cell_type, "X", "Fib_"))

p.sample.comp= ggplot(prop_data, aes(x= cell_type, y= cell_type_prop,color= heart_failure))+
  geom_jitter(size= 0.8)+
  geom_boxplot(alpha= 0.8, outlier.color = NA)+
  facet_wrap(~batch)+
  theme_cowplot()+
  scale_color_manual( values = c("black", "darkgrey"))+
  labs(x= "",y= "composition (%)")+
  theme(axis.text.x = element_text(angle=90, vjust= 0.5))
p.sample.comp

pdf(file = "output/fib_sub_analysis/figures/sample_comps_boxplot.pdf",
    width= 7,height= 6)
p.sample.comp
dev.off()

prop_data <- prop_data %>%
  dplyr::select(sample_id, cell_type, cell_type_prop)

prop_data <- prop_data %>%
  pivot_wider(values_from = cell_type_prop, names_from = cell_type, values_fill = 0) %>%
  column_to_rownames("sample_id") %>%
  as.matrix()

# Check names in aesthetics.R!!!!!
col_list <- list(disease_code = etiology_colors,
                 heart_failure = c("HF" = "black",
                                   "NF" = "darkgrey"),
                 study = study_colors)

sample_anns <- c("heart_failure", "disease_code", "batch")

max_fact <- max(prop_data) %>% max() %>% abs()

col_fun_fact <- circlize::colorRamp2(seq((max_fact + 0.5) * 
                                           -1, max_fact + 0.5, length = 50), 
                                     grDevices::hcl.colors(50, "Blue-Red", rev = F))

row_anns <- rownames(prop_data) %>% 
  tibble::enframe(value = "sample_id") %>% 
  dplyr::select_at("sample_id") %>% 
  dplyr::left_join(meta_data, by = "sample_id") %>% 
  dplyr::select_at(sample_anns) %>% 
  as.data.frame()

row_ha <- ComplexHeatmap::rowAnnotation(df = as.data.frame(row_anns), 
                                        gap = grid::unit(2.5, "mm"), 
                                        border = TRUE,
                                        col = col_list)

comp_clustering <- ComplexHeatmap::Heatmap(scale(prop_data), 
                                           show_column_dend = F, 
                                           show_row_names = F, 
                                           right_annotation = row_ha,
                                           #col = col_fun_fact,
                                           name = "scl. props")

pdf(hmp_plt, height = 5, width = 5.3)

plot(comp_clustering)

dev.off()




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
  #dplyr::filter(!disease_code %in% test_exclusion) %>%
  group_by(batch, cell_type) %>%
  nest() %>%
  dplyr::mutate(HF_diff = map(data, function(dat) {
    
    t.test(value ~ heart_failure, data = dat) %>%
      broom::tidy()
    
  })) %>%
  dplyr::select(-data) %>%
  unnest(HF_diff) 

t_mat= study_diff_stats%>%
  dplyr::select(cell_type, batch, statistic) %>%
  pivot_wider(names_from = batch, values_from = statistic) %>%
  column_to_rownames("cell_type") %>%
  as.matrix()

p_mat= study_diff_stats%>%
  dplyr::select(cell_type, batch, p.value) %>%
  pivot_wider(names_from = batch, values_from = p.value) %>%
  column_to_rownames("cell_type") %>%
  as.matrix()

# Aesthetics

col_fun_fact <- circlize::colorRamp2(seq((max(t_mat) + 0.5) * 
                                           -1, max(t_mat) + 0.5, length = 25), 
                                     grDevices::hcl.colors(25, "Blue-Yellow 3", rev = F))

contrast_diff_plt= Heatmap(t(t_mat), cell_fun = function(j, i, x, y, w, h, fill) {
  if(t(p_mat)[i, j] < 0.001) {
    grid.text("***", x, y)
  } else if(t(p_mat)[i, j] < 0.01) {
    grid.text("**", x, y)
  }
  
},name = "t-value\nHF - NF",
show_column_dend = T,
rect_gp = gpar(col = "white", lwd = 2))

#col = col_fun_fact)
pdf(tval_hmp_plt, height = 2.4, width = 4)
plot(contrast_diff_plt)
dev.off()

saveRDS(study_diff_stats, file= "output/fib_sub_analysis/fib_state_composition_tvals.rds")


# joined model via lmm ----------------------------------------------------

prop_data_clr_tidy <- prop_data_clr %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(-sample_id, names_to = "cell_type") %>%
  left_join(meta_data, by = "sample_id") %>%
  mutate(heart_failure = factor(heart_failure, levels = c("NF", "HF")))

wanted_cells<- unique(prop_data_clr_tidy$cell_type) 

#run lmm per cell type
study_diff_stats_lmer <- map(wanted_cells, function(cell){
  model_df <- prop_data_clr_tidy%>% 
    filter(cell_type ==cell)
  
  fit= lmerTest::lmer(formula = as.formula(paste0("value ~ heart_failure + ( 1 | batch)")),
                      data= model_df)
  model_res <- fit%>% summary()
  
  fixed_summ <- model_res$coefficients %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    dplyr::filter(grepl(pattern = "heart_failureHF", term)) %>%
    dplyr::rename("p_val" = "Pr(>|t|)")
  
  random_summ <- model_res$varcor %>%
    as.data.frame() %>%
    dplyr::select(grp, vcov) %>%
    pivot_wider(names_from = grp, values_from = vcov) %>%
    dplyr::mutate(perc_studyvar = batch/(batch + Residual)) %>%
    pull(perc_studyvar)
  
  all_res <- fixed_summ %>%
    dplyr::mutate(perc_studyvar = random_summ)%>%
    dplyr::mutate(cell_type = cell)
  
  return(all_res)
}) %>%
  do.call(rbind,.)%>%
  ungroup() %>%
  dplyr::mutate(adj_pval = p.adjust(p_val))

study_diff_stats_lmer


# plot results together ---------------------------------------------------

library(ComplexHeatmap)
library(circlize)
library(grid)

study_diff_stats
study_diff_stats_lmer

stats2<-study_diff_stats %>%
  select(batch, cell_type, statistic)%>%
  #filter(signature_matrix_category==sig)%>%
  mutate(statistic= ifelse(is.na(statistic), 0, statistic))%>%
  pivot_wider(names_from = cell_type, values_from= statistic)%>%
  as.data.frame() %>%
  column_to_rownames("batch")%>%
  as.matrix()

p.vals<-study_diff_stats %>%
  dplyr::mutate(adj_pval = p.adjust(p.value))%>%
  select(batch, cell_type, adj_pval)%>%
  #filter(signature_matrix_category==sig)%>%
  mutate(adj_pval= ifelse(is.na(adj_pval), 1, adj_pval))%>%
  pivot_wider(names_from = cell_type, values_from= adj_pval)%>%
  as.data.frame() %>%
  column_to_rownames("batch")%>%
  as.matrix()

study_diff_stats_lmer
color_fun <- function(estimates, p_vals) {
  ifelse(
    p_vals < 0.05,  # Check if p-value is less than 0.05
    ifelse(estimates > 0, 
           color_list$general_up_down_colors["up"], 
           color_list$general_up_down_colors["down"]
    ),  # Positive estimates are blue, negative are red
    "grey"  # Use grey if p-value is not significant
  )
}
color_list<-readRDS("color_list_figures.rds")
heatmap_color_fun <- colorRamp2(
  c(min(stats2), 0, max(stats2)),  # Define the range of the data
  c(color_list$general_up_down_colors["down"], "white", color_list$general_up_down_colors["up"])
)

# Create the bar plot annotation with color coding
row_ha <- rowAnnotation(
  estimate = anno_barplot(
    study_diff_stats_lmer$Estimate,
    gp = gpar(fill = color_fun(study_diff_stats_lmer$Estimate, study_diff_stats_lmer$p_val)),
    ylim = c(min(study_diff_stats_lmer$Estimate, 0), max(study_diff_stats_lmer$Estimate, 0)),
    width = unit(1.8,"cm")
  )
)

star_matrix <- t(ifelse(p.vals < 0.05, "*", ""))

# Create the heatmap with the left annotation and stars
fib_comp_changes<- ComplexHeatmap::Heatmap(
  t(stats2),
  name = "t-value\nHF-NF",
  rect_gp = gpar(col = "black", lwd = 1),
  left_annotation = row_ha,
  col = heatmap_color_fun,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(star_matrix[i, j], x, y, gp = gpar(fontsize = 12, col = "black"))
  }, 
  show_column_dend = F, 
  show_row_dend = F
)
fib_comp_changes  

pdf("output/fib_sub_analysis/figures/sample_comps_t_w_lmm.pdf", 
    width= 2.8, height= 2.5)
  print(fib_comp_changes)
dev.off()
