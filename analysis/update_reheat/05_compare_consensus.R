# HEADER --------------------------------------------
#
# Author:     Jan D. Lanzer
# Copyright   Copyright 2023 - Jan D. Lanzer
# Email:      lanzerjan@gmail.com
#
# Date:     2023-09-21
#
# Script Name:    ~/R-projects/reheat2_pilot/update_reheat/compare_consensus.R
#
# Script Description:
# after update, assess amount of change between rankings

library(tidyverse)
library(cowplot)

meta2021 = readRDS("~/R-projects/Collaborations/cheerio/app_data/study_ranks.rds")
meta2023= read.csv("output/reheat1/meta_analysis_summary.txt", sep = '\t')%>% as_tibble()%>%
  mutate(rank = rank(fisher_pvalue, ties.method = "first"))

directed_signature = readRDS("app_data/signature.rds")

label_genes= c("NPPA", "MYH6","FNDC1", "COL1A1", "NOS3")
## compare pval distribution
p.ranks= rbind(meta2021 %>% mutate(v = "v2021"), 
      meta2023 %>% mutate(v= "v2024"))%>%
  mutate(gene2= as.character(gene),
    label = ifelse(gene2 %in% label_genes & v =="v2024", gene2, "") )%>%
  ggplot(., aes(y= -log10(fisher_pvalue), x= rank, color= v, label =label))+
  geom_point()+
  ggrepel::geom_text_repel(max.overlaps = 5000,
                           color="black",
                           size= 3, 
                           force_pull = 1,nudge_x = 10,
                           min.segment.length = 0.1,
                           show.legend = F)+
  geom_vline(xintercept = 500,lty= 2)+
  geom_hline(yintercept = -log10(10^-5), lty= 2)+
  theme_cowplot()+
  labs(color= "version")+
  scale_color_manual(values=c("grey", "#CD5C5C" ))+
  labs(x= "Consenus gene ranking")
p.ranks


pdf("output/reheat1/fisher_rank.pdf",
    width= 6, 
    height= 3.)
p.ranks
dev.off()

rank1= meta2021
rank2= meta2023 

  rank1 = rank1%>% dplyr::select(gene, rank )%>% dplyr::rename(rank1= rank)
  rank2 = rank2%>% dplyr::select(gene, rank)%>% dplyr::rename(rank2= rank)
  
  df= full_join(rank1, rank2, by= "gene")
  
  
  #test jaccard for the overlap of a sequence of top genes
  
  top_genes= c(50,500, 1000, 2000, 5000, 14000, dim(df)[1])
  x= 200000
  jaccs <- sapply(top_genes, function(x){
    print(x)
    l1 <- df %>% filter(rank1<= x)%>% pull(gene)
    l2 <- df %>% filter(rank2<= x)%>% pull(gene)
    print(length(l2))
    jacc = length(intersect(l1, l2)) / length(union(l1 , l2))
  })
  max(jaccs)
  
  p.ranking_comp <- data.frame(top_genes, jaccs)%>%
    ggplot( aes(x= as.factor(top_genes), y= jaccs))+
    geom_point(size= 3)+
    geom_col(width= 0.04, fill ="black", color="black")+
    labs(x= "n top genes", 
         y= "Jaccard's Index")+
    theme_linedraw()+
    ylim(c(0 , 1))+
    coord_flip()
    
  pdf("output/reheat1/ranking_comparison_jacc.pdf", 
      height= 2, width= 2)
  p.ranking_comp
  dev.off()
  
  ggplot(df %>% filter(rank1<500 | rank2 <500) , aes(x= rank1, y= rank2))+
    geom_point(size= 0.1)+
    geom_abline(slope= 1)+
    geom_hline(yintercept = 500)+
    geom_vline(xintercept = 500)
  
  res<-cor.test(df$rank1,df$rank2, method="kendall")
  res
  
