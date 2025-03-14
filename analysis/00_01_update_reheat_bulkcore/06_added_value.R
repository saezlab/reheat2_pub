# HEADER --------------------------------------------
#
# Author:     Jan D. Lanzer
# Copyright   Copyright 2024 - Jan D. Lanzer
# Email:      lanzerjan@gmail.com
#
# Date:     2024-06-28
#
# Script Name:    ~/R-projects/reheat2_pilot/update_reheat/added_value.R
#
# Script Description:
#
#' In this script we test the differences in the predictions between the signature observed in our meta-analysis
#' and individual signatures from experiments


library(ggpubr)
library(cowplot)
library(broom)
library(reshape2)
library(ROCR)
library(fgsea)
library(tidyverse)

source("../HF_meta-analysis/src/data_utils.R") #general functions 
source("../HF_meta-analysis/src/misc_utils.R")

#METAheart = readRDS(file = "../HF_meta-analysis/data/METAheart.rds") #main object

# METAheart = readRDS(file = "output/reheat1/METAheart2023.rds") #main object
# METAheart$Magnet<-c()
# METAheart$Hua19<-c()
meta2021 = readRDS("~/R-projects/Collaborations/cheerio/app_data/study_ranks.rds")
meta2023= read.csv("output/reheat1/meta_analysis_summary.txt", sep = '\t')%>% as_tibble()%>%
  mutate(rank = rank(fisher_pvalue, ties.method = "first"))%>%
  dplyr::select(rank, everything())
meta2023 %>% pull(gene)%>% unique() %>% length()
#library(WriteXLS)


# function ----------------------------------------------------------------

#' @description, this function runs two metrics for pairwise study comparison. 
#' first it assesses an AUC to test how well samples can be separated and second it
#' calculates ES to on top up or downregulated genes. 
#' Both metrics are calculated in two scenarios, A) by taking top genes from one study 
#' to another, and B) taking the top consensus genes. 
#' The top genes are passed as a named list that contain the ranks one wants to assess.
#' This enables to check also how lower ranked genes perform without using top performing genes.
#' However, this feature is not used for the pariwise study comparison, there the highes value of the
#' passed ranks is used to select 1:n genes.
#' 

calc_added_value <- function(METAheart, top_genes, meta_rank, experiments=NULL){
  
  if(is.null(experiments)){
    experiments = names(METAheart)
    names(experiments) = experiments
  }else{
    METAheart<- METAheart[experiments]
  }
  
  # For labeling
  t_matrix = get_all_limma(meta_list = METAheart,
                           limma_column = "t")
  
  
  #### PART 1  AUC calculation on disease score
  # calculate the pariwise disease classifier AUCs on the meta heart studies
  pairwise_AUC <- lapply(names(top_genes), function(n){
    print(paste0("calculate AUC for pairwise studies on top ",max(top_genes[[n]]), " genes.."))
    # Here we get AUC for pairwise classifiers
    pairwise_500 = pairwise_ds(experiments = experiments,
                               meta_list = METAheart,
                               t_matrix = t_matrix,
                               ngenes = max(top_genes[[n]])) %>%
      mutate(top_n = n)#Second page excel
    
  })%>% do.call(rbind, .)
  
  # #extract the fisher rank 
  # fisher_rank = run_fisher_meta(meta_list = METAheart,
  #                               n_missing = length(METAheart) - 10)
  # genes = names(fisher_rank)
  
  
  # calculate the disease classifier AUCs using the top n genes from consensus
  meta_pred <- lapply(names(top_genes), function(n){
    
    print(paste0("calculate AUC for meta score on top ",n, " genes.."))
    genes = meta_rank %>% filter(rank %in% top_genes[[n]])%>% pull(gene)
    print(length(genes))
    ds_top = getRisk_Stats_v2(Experiment_List = METAheart,
                              limma_t_mat = t_matrix, 
                              genes = as.character(genes)
                              )
    
    ds_top_predictions = enframe(lapply(ds_top, function(x) {
      enframe(x[["SingleAUC"]])
      })) %>% unnest()%>%
      mutate(top_n = n)
    
    colnames(ds_top_predictions) = c("PredictedExperiment", 
                                     "PredictorExperiment", 
                                     "meta_auc", 
                                     "top_n")
    return(ds_top_predictions)
  })%>% do.call(rbind, .)
  
  # Merge both: 
  comp_df = left_join(pairwise_AUC,
                      meta_pred) %>%
            dplyr::filter(PredictedExperiment != PredictorExperiment) %>%
    dplyr::select(PredictorExperiment, 
                  PredictedExperiment,
                  single, meta_auc, top_n)
  
  #### PART 2  ES score calculation for consistency
  ES_scores <- lapply(names(top_genes), function(n){
      
    print(paste0("calculate ES on top ",max(top_genes[[n]]), " genes.."))
    #x= METAheart$Kim16
    study_deg_list_up = lapply(METAheart, function(x){
      deg = x$HF_limma %>%  
        filter(t >0 )%>% 
        slice_min(order_by = adj.P.Val, n = max(top_genes[[n]]))
      return(deg[[1]])
    })
    
    study_deg_list_down = lapply(METAheart, function(x){
      # deg = dplyr::filter(x$HF_limma,
      #                     ID %in% names(fisher_rank[1:max(top_genes[[n]])])) %>%
      #   filter(t < 0)
      deg = x$HF_limma %>%  
        filter(t <0 )%>% 
        slice_min(order_by = adj.P.Val, n = max(top_genes[[n]]))
      return(deg[[1]])
    })
    
    # we will add the top n genes from the meta analysis: 
    
    study_deg_list_up$meta =  meta_rank %>% 
      filter(mean_t>0, 
             rank %in% top_genes[[n]])%>%
      pull(gene)
    
    study_deg_list_down$meta = meta_rank %>%
      filter(mean_t<0, 
             rank %in% top_genes[[n]])%>%
      pull(gene)
    
    # upregulation ES  on pairwise study comparison
    
    up_ES = lapply(experiments, function(x){
      
      stat_rank = METAheart[[x]][["HF_limma"]][["t"]]
      names(stat_rank) = METAheart[[x]][["HF_limma"]][["ID"]]
      stat_rank = sort(stat_rank)
      set.seed(1234)
      
      up_row = as_tibble(fgsea(pathways = study_deg_list_up,
                               stats = stat_rank, nperm = 1)) %>%
        dplyr::select(pathway,ES)
    })
    
    up_ES = up_ES %>% 
      enframe("Reference") %>% unnest()
    
    colnames(up_ES) = c("PredictorExperiment","PredictedExperiment","up_ES")
    
    down_ES = lapply(experiments, function(x){
      
      stat_rank = METAheart[[x]][["HF_limma"]][["t"]]
      names(stat_rank) = METAheart[[x]][["HF_limma"]][["ID"]]
      stat_rank = sort(stat_rank)
      set.seed(1234)
      
      up_row = as_tibble(fgsea(pathways = study_deg_list_down,
                               stats = stat_rank,nperm = 1)) %>%
        dplyr::select(pathway,ES)
    })
    
    down_ES = down_ES %>% 
      enframe("Reference") %>% unnest()
    
    colnames(down_ES) = c("PredictorExperiment","PredictedExperiment","down_ES")
    
    return(left_join(down_ES, up_ES)%>% mutate(top_n= n))
    
  })%>%
    do.call(rbind, .)
  
  
  res_df<- comp_df %>% 
    full_join(ES_scores, by= c("PredictorExperiment", 
                               "PredictedExperiment",
                               "top_n"))
  
  return(res_df)
}


# run added value for our two rankings (2021 and 2023) --------------------

#top_genes <- list(100, 500, 1000, 2000, 5000,10000)

top_genes <- list("1_500"= 1:500,
                  "500_2000" =500:2000, 
                  "2000_5000" = 2000:5000, 
                  "5000_10000" = 5000:10000)

# top_genes <- list("1_500"= 1:500,
#                   "5000_7000" = 5000:7000)

set.seed(1)
#results based on meta_ranking 2021

experiments = names(METAheart)
names(experiments) = experiments

# to avoid overfitting (testing old signature on unseen data)
# deselect new studies
experiments21 <- experiments[!names(experiments) %in% c("Wang22", "Forte22", "Rao21", "Flam19")]

res2021 <- calc_added_value(METAheart = METAheart,
                            top_genes = top_genes,
                            meta_rank=meta2021,
                            experiments =experiments21 )%>%
  mutate(version= "2021")

#results based on updated meta_ranking 2023
res2023 <- calc_added_value(METAheart = METAheart, 
                            top_genes = top_genes, 
                            meta_rank=meta2023, 
                            experiments = experiments)%>%
  mutate(version= "2024")

res_comp = rbind(res2021, res2023)


saveRDS(res_comp, "output/reheat1/added_value_combined_gene_ranged.rds")
res_comp <- readRDS("output/reheat1/added_value_combined_gene_ranged.rds")

# stats test for differences (improvements over old consensus) ------------

# For each question we have two metrics to test for, the ES and the AUC. 
# We will test first whether the new consensus improved over the old consensus
# for both metrics (ignoring the pairwise comparisons)

print("Are AUCs better?")
#test
res_comp <- res_comp %>%# filter(top_n != 5000)%>%
  mutate(version= factor(version, levels=c("2024", "2021")))%>% 
  mutate(top_n= factor(top_n, levels= names(top_genes)))%>% 
  dplyr::rename(down_regulated_genes = "down_ES", 
         up_regulated_genes= up_ES)

res_comp %>% 
  distinct(PredictorExperiment, PredictedExperiment, meta_auc, top_n, version)%>%
  dplyr::select(-PredictorExperiment, -PredictedExperiment)%>%
  drop_na()%>%
  nest(data = c(-top_n)) %>%
  mutate(data = map(data, ~ t.test(meta_auc ~ version, data = .x, alternative = "two.sided")),
         data = map(data, tidy)) %>% 
  unnest(data)
# nope p-values are >0.8

#plot=
p.ranking_comp_auc<- res_comp %>% #filter(top_n != 5000)%>%
  filter(PredictorExperiment!= PredictedExperiment)%>%
  distinct(PredictorExperiment, PredictedExperiment, meta_auc, top_n, version)%>%
  dplyr::select(-PredictorExperiment, -PredictedExperiment)%>%
  mutate(top_n= str_replace_all(top_n, "_", ":"),
         top_n = factor(top_n, 
                        levels =c("1:500", "500:2000", "2000:5000", "5000:10000"))
         )%>%
  ggplot(aes(x= top_n, color = version, y= meta_auc))+
  #geom_jitter()+
  geom_boxplot(width= 0.7)+
  #geom_violin()+
  scale_fill_brewer(type = "qual",palette = 2)+
  theme_cowplot()+
  #ggpubr::stat_compare_means(method = "wilcox.test", alternative = "two.sided")+
  scale_color_brewer(type="qual", palette=2)+
  labs(x= "gene ranking in consensus signature", 
       y= "Disease score AUC")+
  theme(axis.text.x = element_text(angle=40, hjust= 1))
p.ranking_comp_auc
print("Are ESs better?")

#stats 
res_comp %>% 
  distinct(PredictorExperiment, PredictedExperiment,
           down_regulated_genes,  top_n, version)%>%
  filter(PredictedExperiment=="meta")%>%
  dplyr::select(-PredictorExperiment, -PredictedExperiment)%>%
  nest(data = c(-top_n)) %>%
  mutate(data = map(data, ~ t.test(down_regulated_genes ~ version, data = .x)),
         data = map(data, tidy)) %>% 
  unnest(data)

res_comp %>% 
  distinct(PredictorExperiment, PredictedExperiment, up_regulated_genes, top_n, version)%>%
  filter(PredictedExperiment=="meta")%>%
  dplyr::select(-PredictorExperiment, -PredictedExperiment)%>%
  nest(data = c(-top_n)) %>%
  mutate(data = map(data, ~ t.test(up_regulated_genes ~ version, data = .x)),
         data = map(data, tidy)) %>% 
  unnest(data)

res_comp %>% group_by(version, top_n)%>% summarise(m= mean(meta_auc, na.rm= T), 
                                                   dn= mean(down_regulated_genes, na.rm= T),
                                                   up= mean(up_regulated_genes, na.rm= T))
#plot: 
p.ranking_comp_ES<- res_comp %>% 
  distinct(PredictorExperiment, PredictedExperiment, down_regulated_genes,up_regulated_genes, top_n, version)%>%
  pivot_longer(cols = c(down_regulated_genes, up_regulated_genes) , names_to = "ES_sign", values_to= "ES")%>%
  filter(PredictedExperiment=="meta")%>%
  mutate(top_n= str_replace_all(top_n, "_", ":"),
         top_n = factor(top_n, 
                        levels =c("1:500", "500:2000", "2000:5000", "5000:10000"))
  )%>%
  dplyr::select(-PredictorExperiment, -PredictedExperiment)%>%
  
  ggplot(aes(x= ES_sign, y= ES, color = version))+
  geom_boxplot()+
  facet_grid(~top_n)+
  theme_cowplot()+
  theme(axis.text.x= element_blank(),
        axis.ticks.x = element_blank())+
  scale_color_brewer(type="qual", palette=2)+
  labs(x= "up or down regulated genes")
  

p <- plot_grid(p.ranking_comp_auc+labs(x= "", y= "Disease Score\nAUROC"), p.ranking_comp_ES, ncol = 1, align = "v")

pdf("output/reheat1/ranking_comparison_es_auc.pdf", 
    width = 6, height= 5)
p
dev.off()
# compare number of DEG at alpha
alpha = 10^-5 
meta2021%>% filter(fisher_pvalue< alpha)%>% pull(gene)%>% length()
meta2023%>% filter(fisher_pvalue< alpha)%>% pull(gene)%>% length()

#compare coverage 
dim(meta2021)
dim(meta2023)
