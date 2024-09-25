# HEADER --------------------------------------------
#
# Author:     Jan D. Lanzer
# Copyright   Copyright 2024 - Jan D. Lanzer
# Email:      lanzerjan@gmail.com
#
# Date:     2024-01-12
#
# Script Name:    
#
# Script Description:
# Fibroblast state characterization
# 1. Use cell state marker for interpretation
# 2. Use the pseudobulked profiles of states to identify characteristics

library(msigdbr)
library(decoupleR)
library(tidyverse)
library(ComplexHeatmap)
library(cowplot)
library(tune)
library(lme4)
library(lmerTest)
library("r2glmm")
library(ggrepel)
library(rstatix)
library(msigdbr)

# data -------------------------------------------------------------------
#read.csv("data/metadata/Reichart2022_DCM_RVsamples.csv")
obs <- read.csv("output/fib_sub_analysis/all_state_coldata_sub_reg.csv")

pb = readRDS( "output/fib_sub_analysis/all_state_pb_processed_sub_reg.rds")

state.marker = read_csv("output/fib_sub_analysis/cluster_markers_processed.csv")

counts_cell = read_csv("output/fib_sub_analysis/cluster_counts_percluster_reg.csv")

cell_state= "leiden_0.5"

## calculate how many patients contributed cells to each state

counts_cell<-counts_cell%>% dplyr::rename(cell_state = "leiden_0.5")
summed_per_patient= counts_cell %>% distinct(sample_id, cell_state, cell_count)%>%
  group_by(sample_id, cell_state)%>%
  summarise(cell_count= sum(cell_count))

summed_per_patient%>%
  group_by(sample_id, cell_state)%>%
  summarise(s= cell_count==0)%>%
  ungroup()%>% pull(s)%>% table()

summed_per_patient%>%ggplot(aes(x= factor(cell_state), fill= cell_count))+
  geom_bar()


# prepare data -------------------------------------------------------------------------
# extract sample level target matrix for all patients:
pb_target <- pb %>% lapply(., function(x){
  x %>% 
    distinct(sample_state,sample_id, cell_state, heart_failure, disease_code, batch)
  
})%>% do.call(rbind, .)%>%
  mutate(cell_state= factor(paste0("Fib_", cell_state)))

#1 get exp matrices per study (different feature space for each study)
list.mats <- pb %>% lapply(., function(x){
  x %>% 
    select(feature, sample_state, value)%>%
    pivot_wider(names_from = sample_state, values_from= value)%>%
    as.data.frame()%>%
    #drop_na() %>%
    column_to_rownames("feature")%>% 
    as.matrix()
})

lapply(list.mats,dim)

#2 get exp joint matrix for all studies, (same feature space)
pb_long <- pb %>% lapply(., function(x){
  x %>% 
    select(feature, sample_state, value)
})%>% do.call(rbind, .)

mat <- pb_long%>%
  pivot_wider(names_from = sample_state, values_from= value)%>%
  as.data.frame()%>%
  drop_na() %>%
  column_to_rownames("feature")%>% 
  as.matrix()

dim(mat)

write.csv(mat, "output/pb_state_patient_profiles.csv")
write.csv(pb_target, "output/pb_state_patient_profiles_target.csv")
# load PK -----------------------------------------------------------------
# pways
net_prog  <- get_progeny(organism = 'human', top = 500)
#msig
net_cur <- msigdbr::msigdbr(category = "C2")%>%
  select(gs_name,human_gene_symbol )%>%
  dplyr::rename(source= gs_name, target= human_gene_symbol)%>%
  distinct()

net_HM <- msigdbr::msigdbr(category = "H")%>%
  select(gs_name,human_gene_symbol )%>%
  dplyr::rename(source= gs_name, target= human_gene_symbol)%>%
  distinct()%>% 
  mutate(source= str_replace_all(source, "HALLMARK_", ""))

net_naba<- net_cur %>% filter(grepl("NABA",source))%>%
  mutate(source = str_replace_all(source, "NABA_", ""))

#cytokines
cytosig = read_csv("data/prior_knowledge/immune.csv")
ILs <- cytosig$source %>% unique()
ILs <- enframe(ILs, value = "target", name= "source")%>%
  mutate(source= "Cytokine_espression", weight=0.5 )%>%
  select(target, source, weight)
net_cytosig<- rbind(cytosig, 
                ILs)

# -------------------------------------------------------------------------
# function to
# 1. enrich genes via ulm
# 2. select gsets with high variance between states (anova p-value), 
# 3. do a one vs. all t-test on those gsets, correct p-values, 
# 4. plot as hmap
##
get_enrich<- function(mat, alpha, net, ...){
  # 1. enrich genes via ulm
  res= run_ulm(mat,   net, ...)
  
  res <- res%>%
    left_join(., pb_target,by= join_by(condition == sample_state))
  
  res <- res %>% 
    mutate(p.adj= p.adjust(p_value))%>% 
    group_by(source, batch)%>% 
    mutate(score_scaled= scale(score))
  
  # 2. select gsets with high variance between states (anova p-value), 
  pways<- res %>%ungroup() %>%
    group_by(source)%>%
    anova_test(score_scaled~cell_state)%>%
    as_tibble()%>%
    filter(p<alpha)%>% 
    arrange(p)%>% 
    pull(source)
  
  print(length(pways))
  # 3. do a one vs. all t-test on those gsets, correct p-values, 
  t_test_results <- res %>%
    ungroup()%>%
    filter(source %in% pways)%>%
    distinct(cell_state, source) %>%
    pmap_df(~ {
      current_group <- res %>%
        filter(cell_state == ..1, source == ..2) %>%
        pull(score_scaled)
      
      rest_groups <- res %>%
        filter(cell_state != ..1 , source == ..2) %>%
        pull(score_scaled)
      
      t_test_result <- t.test(x= current_group,y = rest_groups, alternative = "greater",paired = F )
      
      tibble(
        cell_state = ..1,
        source = ..2,
        mean_group = mean(current_group),
        mean_rest = mean(rest_groups),
        mean_diff = mean_group-mean_rest, 
        t_value = t_test_result$statistic,
        p_value = t_test_result$p.value
      )
    })
  
  
  
  #return(res%>% mutate(significant= source %in% pways))
  return(t_test_results)
}

df.cytos <- get_enrich(mat, 0.05, net_cytosig, .mor ="weight")
df.prog<- get_enrich(mat, 0.05, net_prog, .mor ="weight")
df.cur<- get_enrich(mat, 0.05, net_cur)
df.h<- get_enrich(mat, 0.05, net_HM)
df.n<- get_enrich(mat, 0.05, net_naba)

x<-c(#naba
  "BASEMENT_MEMBRANES", "COLLAGENS", "CORE_MATRISOME", "SECRETED_FACTORS",
  #progeny
  "TGFb", "TNFa",
  #cytokines
  "BMP6","Activin A", "BMP4", "IL6", "IL1A", "IL1B", "IL22", "IL2",
  #hmarks
  "ANGIOGENESIS","HEDGEHOG_SIGNALING", "IL6_JAK_STAT3_SIGNALING", "PROTEIN_SECRETION", "COAGULATION"
  )
df.h %>% filter(grepl("HED", source))
df<- rbind(df.h, df.n, df.prog, df.cytos) %>% 
  mutate(p_adj= p.adjust(p_value))%>%
  filter(source %in% x)

plot_hmap<-function(df){
  
  mat2<- df%>%
    #mutate(star= ifelse(p_value<0.05 & t_value>0 , "*", ""))%>%
    select(cell_state, source, p_adj)%>%
    #mutate(p_value= p.adjust(p_value))%>%
    pivot_wider(names_from = source, values_from = p_adj)%>%
    as.data.frame() %>% column_to_rownames("cell_state")
  
  X <- 
    df %>%  
    ungroup()%>% 
    select(mean_group, source, cell_state )%>%
    group_by(source, cell_state)%>%
    pivot_wider(names_from = source, values_from = mean_group)%>%
    as.data.frame() %>% 
    column_to_rownames("cell_state")%>% 
    scale()
  
  hmap.pways<-Heatmap(X, border_gp = gpar(col = "black", lty = 1),
                      rect_gp = gpar(col = "white", lwd = 1), 
                      name= "Enrichment\nscore", 
                      cluster_rows = F,
                      show_row_dend = F,
                      show_column_dend = F,
                      #row_dend_side = "right", 
                      row_names_side = "left",
                      column_names_side = "top",
                      row_title_side = "left",
                      cell_fun = function(j, i, x, y, w, h, fill) {
                        if(mat2[i, j] < 0.001) {
                          grid.text("**", x, y)
                        } else if(mat2[i, j] < 0.01) {
                          grid.text("*", x, y)
                        }
                        
                      }
  )
  
  
  print(hmap.pways)
}

p.hmap_selected<-plot_hmap(df)
p.hmap_selected
pdf("output/figures/fibstates_selected_pways.pdf", 
    width= 5, height= 3.5)
p.hmap_selected
dev.off()


map(list(df.h, df.n, df.prog, df.cytos), plot_hmap)


state.marker%>% filter(cluster=="Fib_5")
