# HEADER --------------------------------------------
#
# Author:     Jan D. Lanzer
# Copyright   Copyright 2024 - Jan D. Lanzer
# Email:      lanzerjan@gmail.com
#
# Date:     2024-06-13
#
# Script Name:    
#
# Script Description:
# compare the differntial expressed genes of val studies with the factors 

library(msigdbr)
library(decoupleR)
library(cowplot)
library(tidyverse)
library(ComplexHeatmap)

toptables= rbind("LVAD"= readRDS("data/additional_studies/DGE_celltypes_AMRUTE.rds")%>%
                  mutate(dataset= "LVAD"), 
     "FET"= readRDS("data/additional_studies/DGE_celltypes_FETAL.rds")%>%
       mutate(dataset ="FETAL")
     )

g.load = read_csv("output/mofa/gene_loadings.csv")

g.load

j.toptables<- 
toptables%>% dplyr::rename(ctype = clust,feature= gene)%>%
  left_join(g.load%>% filter(Factor=="Factor1"), by= c("feature", "ctype"))

# explore associations  ---------------------------------------------------
# number of significant genes to calculate a corr coeff
n<- 10

result <- j.toptables %>%
  #filter(coef %in% c("FET_v_CT", "pre_r_v_pre_nr")) %>%
  filter(adj.P.Val < 0.1) %>% #filter for dge here 
  group_by(coef, ctype) %>%
  filter(n() > n) %>%  # filter groups with more than n data points
  summarise(cor = cor.test(value, t)$estimate, 
            cor.p = cor.test(value, t)$p.value)

p.corr_overview<- result %>%
  mutate(label=ifelse(cor.p<0.01, "*", ""), 
         data=ifelse(grepl("p", coef), "Amrute19", "Mehdiabadi22"))%>%
  ggplot(., aes(y= ctype, fill= cor, x = coef, label = label))+
  geom_tile(color="black")+
  geom_text()+
  facet_grid(~data,scales = "free")+
  theme_cowplot()+
  scale_fill_gradient2(low= col_list$general_up_down_colors[1], mid="white", high = col_list$general_up_down_colors[2])+
  theme(axis.text.x = element_text(angle= 45,hjust= 1), 
        strip.text=element_text(size=9))+
  labs(y="", x= "contrast")

pdf("output/figures/project_correlation_ctypes.pdf", 
    width = 4 ,
    height= 3 )
p.corr_overview
dev.off()


p.LVAD<- j.top %>%
  group_by(ctype)%>%
  mutate(fact1 = ifelse(sign(pre_r_v_pre_nr)!= sign(value) & abs(value)>0.1, "fact1", "other"))%>%
  ggplot(aes(x=value , y= pre_r_v_pre_nr, col= fact1))+
  geom_point(size= 0.4, alpha = 0.8,show.legend = F)+
  facet_grid(~ctype)+
  theme_cowplot()+
  scale_color_manual(values= c( "darkgreen","grey"))+
  labs(x= "Factor1 loading", y= "Rec vs Nonrec\n(t-statistic)")
p.LVAD_red<- j.top %>%
  group_by(ctype)%>%filter(ctype %in% c("CM", "Fib"))%>%
  mutate(fact1 = ifelse(sign(pre_r_v_pre_nr)!= sign(value) & abs(value)>0.1, "fact1", "other"))%>%
  ggplot(aes(x=value , y= pre_r_v_pre_nr))+
  geom_point(size= 0.4, alpha = 0.8,show.legend = F)+
  facet_grid(~ctype)+
  theme_cowplot()+
  stat_cor()+
  scale_color_manual(values= c( "darkgreen","grey"))+
  labs(x= "Factor1 loading", y= "Rec vs Nonrec\n(t-statistic)")
p.LVAD_red
p.FET<-j.top %>%
  group_by(ctype)%>%
  mutate(fact1 = ifelse(sign(FET_v_CT)!= sign(value) & abs(value)>0.1  , "fact1", "other"))%>%
  ggplot(aes(x=value , y= FET_v_CT, col= fact1))+
  geom_point(size= 0.4, alpha = 0.8,show.legend = F)+
  facet_grid(~ctype)+
  theme_cowplot()+
  scale_color_manual(values= c( "darkgreen","grey"))+
  labs(x= "Factor1 loading", y= "Fetal vs CT \n(t-statistic)")

sig.genes = j.toptables %>% filter(ctype %in% c("CM", "Fib"), coef== "FET_v_CT", adj.P.Val <0.01)%>% pull(feature)
p.FET_red <- j.top %>%
  group_by(ctype)%>%
  filter(ctype %in% c("CM", "Fib"))%>%
  mutate(fact1 = ifelse(sign(FET_v_CT)!= sign(value) & abs(value)>0.1  , "fact1", "other"))%>%
  filter(feature %in% sig.genes)%>%
  ggplot(aes(x=value , y= FET_v_CT, col= fact1))+
  geom_point(size= 0.4, alpha = 0.8,show.legend = F)+
  facet_grid(~ctype)+
  geom_smooth(method = "lm", color="black")+
  theme_cowplot()+
  scale_color_manual(values= c( "darkgreen","grey"))+
  labs(x= "Factor1 loading", y= "Fetal vs CT \n(t-statistic)")
p.FET_red

p.full<- plot_grid(p.LVAD+labs(x= ""), p.FET, nrow= 2)
p.red<- plot_grid(p.LVAD_red+labs(x= ""), p.FET_red, nrow= 2)
#calculate correlation per cell type
  #filter(ctype=="Myeloid", coef=="FET_v_DCM")%>%


result

# characterize vectors with gsea ----------------------------------------------------
source("utils1.R")
# get
net_list <- get_pk_nets()
net_list$tfs <- NULL
toptables$coef%>% unique()

names(net_list)
# function to filter on agreeing regulatory events of factor 

get_enrich_on_fact1<- function(toptable, coefficient, diagonal=TRUE){
  filtered.table <- 
    toptable %>% 
    filter(coef == coefficient)%>%
    dplyr::rename(ctype = clust,feature= gene)%>%
    left_join(g.load%>% filter(Factor=="Factor1"), by= c("feature", "ctype"))
  
  if(diagonal){
    filtered.table <- filtered.table %>%
      group_by(ctype)%>%
      mutate(fact1 = ifelse(sign(t)!= sign(value) & abs(value)>0.05, "fact1", "other"))%>%
      filter(abs(value)>0.1 & fact1=="fact1")
  }else{
    filtered.table <- filtered.table %>%
      group_by(ctype)%>%
      mutate(fact1 = ifelse(sign(t)== sign(value) & abs(value)>0.05, "fact1", "other"))%>%
      filter(abs(value)>0.1 & fact1=="fact1")
  }
  mat <- filtered.table%>% 
    select(feature, t, ctype)%>%
    pivot_wider(names_from = ctype, values_from = t, values_fill = 0)%>%
    as.data.frame()%>%
    column_to_rownames("feature")%>%as.matrix()
  lapply(names(net_list), function(netz){
    res <- run_ulm(mat, net_list[[netz]]%>% 
                     {if (!"weight" %in% colnames(.)) mutate(., weight = 1) else .} ,
                   .mor = "weight")%>%
      mutate(net = netz, p.adj= p.adjust(p_value))
    
  })%>% do.call(rbind, .)%>% mutate(Factor_agreement = diagonal)
}
#x2<-get_enrich_on_fact1(toptable = toptables, "pre_r_v_pre_nr")

fet_ct<-get_enrich_on_fact1(toptable = toptables, "FET_v_CT")%>% 
  mutate(coef ="FET_v_CT")

fet_ct2<-get_enrich_on_fact1(toptable = toptables, "FET_v_CT", diagonal = FALSE)%>% 
  mutate(coef ="FET_v_CT")

lvad<-get_enrich_on_fact1(toptable = toptables, "pre_r_v_pre_nr", diagonal = TRUE)%>% 
  mutate(coef ="pre_r_v_pre_nr")

lvad2<-get_enrich_on_fact1(toptable = toptables, "pre_r_v_pre_nr", diagonal = FALSE)%>% 
  mutate(coef ="pre_r_v_pre_nr")

#collect all in one df
df<- rbind(fet_ct, fet_ct2, lvad, lvad2)

y= "Myeloid"
## FET
plot_enriched<- function(coef_id= "FET_v_CT", alpha= 0.1 ){
  
  map(unique(df$condition), function(y){
    pways<- df %>%
      filter(condition==y)%>%
      ungroup() %>%
      filter(coef== coef_id & Factor_agreement, p.adj<alpha) %>%
      group_by(net, coef)%>%
      slice_max(order_by = abs(score), n= 4)%>% 
      pull(source)
    pways
    df %>% 
      filter(condition==y, source %in% pways)%>%
      filter(coef== coef_id & Factor_agreement) %>%
      ggplot(aes(y= reorder(source,score), x= score, fill = score))+
      geom_col()+
      scale_fill_gradient2()+
      ggtitle(y)+
      theme_cowplot()+
      labs(y= "", x="ES score")
  })
}

plot_enriched()
plot_enriched(coef_id = "pre_r_v_pre_nr")


CM_fetal<- c("GOBP_FATTY_ACID_BETA_OXIDATION","HALLMARK_OXIDATIVE_PHOSPHORYLATION", "BMP6", "PPARA", "PPARD", "NR3C1")
Fib_fetal<- c("NABA_BASEMENT_MEMBRANES", 
              "WP_HEDGEHOG_SIGNALING_PATHWAY_WP4249",
              "GOBP_COLLAGEN_FIBRIL_ORGANIZATION",
              "IL4",
              "TGFB1","TNFA","GLI3")

CM_lvad<- c("NO","IL17A", "GOBP_CELLULAR_RESPONSE_TO_NITROGEN_COMPOUND")
Fib_lvad<- c("NABA_CORE_MATRISOME",
             #"BMP6",
             "HALLMARK_INFLAMMATORY_RESPONSE",
             "IL1B",
             "EGF",
             "IL17A", "IL2", "TNFA" )
Endo_lvad<- c("HALLMARK_IL6_JAK_STAT3_SIGNALING",
             #"BMP6",
             "HALLMARK_INTERFERON_ALPHA_RESPONSE",
             "IL1B",
             "EGF",
             "IL17A", "IL2", "TNFA" )
Myeloid_lvad<- c("IL12", "IL15", "REACTOME_NEUTROPHIL_DEGRANULATION", "HALLMARK_INFLAMMATORY_RESPONSE")


# attempt at hmap
p.lvad.sets<- df %>% 
  filter(condition %in% c("Myeloid", "Endo","Fib"), 
         source %in% c(Myeloid_lvad,Endo_lvad,Fib_lvad), 
         Factor_agreement)%>%
  mutate(label = ifelse(p.adj<0.05, "*", ""))%>%
  #filter(condition %in% c("Fib"), source %in% pways$source)%>%
  filter(coef== "pre_r_v_pre_nr", net!= "naba") %>%
  ggplot(aes(x= reorder(source,score), y= condition, fill = score))+
  geom_tile(color="black", size= 0.5)+
  geom_text(aes(label = label), size= 7)+
  scale_fill_gradient2(low= col_list$general_up_down_colors[[1]], 
                       high = col_list$general_up_down_colors[[2]])+
  theme_cowplot()+
  theme(axis.line = element_blank(), 
        axis.text.x = element_text(angle= 90,vjust= 0.5, hjust= 1))+
  labs(x= "", y="", fill ="Enrichment\nscore")
p.fetal.sets<- df %>% 
  filter(condition %in% c("CM", "Fib"), 
         source %in% c(CM_fetal, Fib_fetal), 
         Factor_agreement)%>%
  mutate(label = ifelse(p.adj<0.05, "*", ""))%>%
  #filter(condition %in% c("Fib"), source %in% pways$source)%>%
  filter(coef== "FET_v_CT", net!= "naba") %>%
  ggplot(aes(x= reorder(source,score), y= condition, fill = score))+
  geom_tile(color="black", size= 0.5)+
  geom_text(aes(label = label), size= 7)+
  scale_fill_gradient2(low= col_list$general_up_down_colors[[1]], 
                       high = col_list$general_up_down_colors[[2]])+
  theme_cowplot()+
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle= 90, vjust= 0.5,hjust= 1))+
  labs(x= "", y="", fill ="Enrichment\nscore")

p.comp<- cowplot::plot_grid(p.lvad.sets+coord_equal(), 
                   p.fetal.sets+coord_equal(), ncol = 1
                  )

pdf("output/figures/project_interpret_shared_signal.pdf", 
    width = 7, height= 10)
p.lvad.sets+coord_equal()
p.fetal.sets+coord_equal()
dev.off()

p.fetal_CM <- df %>% 
  filter(condition %in% c("CM"), source %in% CM_fetal)%>%
  filter(coef== "FET_v_CT") %>%
  ggplot(aes(y= reorder(source,-score), x= score, fill = Factor_agreement))+
  geom_col(color= "black", width= 0.8, position = "dodge")+
  scale_fill_manual(values = c("grey", "darkgreen"))+
  #scale_fill_gradient2(low= "blue", high = "red")+
  ggtitle("CM")+
  theme_cowplot()+
  labs(y= "", x="ES score")+
  coord_flip()+
  theme(axis.text.x = element_text(angle= 90, hjust= 1),
        legend.position = "none")
p.fetal_CM

p.fetal_Fib <- 
  df %>% 
  filter(condition %in% c("Fib"), source %in% Fib_fetal)%>%
  filter(coef== "FET_v_CT", net!= "naba") %>%
  ggplot(aes(y= reorder(source,-score), x= score, fill = Factor_agreement))+
  geom_col(color= "black", width= 0.8, position = "dodge")+
  #scale_fill_gradient2(low= "blue", high = "red")+
  scale_fill_manual(values = c("grey", "darkgreen"))+
  ggtitle("Fib")+
  theme_cowplot()+
  labs(y= "", x="ES score")+
  coord_flip()+
  theme(axis.text.x = element_text(angle= 90, hjust= 1),
        legend.position = "none")
p.fetal_Fib


p.lvad_CM <- 
  df %>% 
  filter(condition %in% c("CM"), source %in% CM_lvad)%>%
  filter(coef== "pre_r_v_pre_nr") %>%
  ggplot(aes(y= reorder(source,-score), x= score, fill = Factor_agreement))+
  geom_col(color= "black", width= 0.8, position = "dodge")+
  scale_fill_manual(values = c("grey", "darkgreen"))+
  #scale_fill_gradient2(low= "blue", high = "red")+
  ggtitle("CM")+
  theme_cowplot()+
  labs(y= "", x="ES score")+
  coord_flip()+
  theme(axis.text.x = element_text(angle= 90, hjust= 1),
        legend.position = "none")
p.lvad_CM

p.lvad_Fib <- 
  df %>% 
  filter(condition %in% c("Fib"), source %in% Fib_lvad)%>%
  #filter(condition %in% c("Fib"), source %in% pways$source)%>%
  filter(coef== "pre_r_v_pre_nr", net!= "naba") %>%
  ggplot(aes(y= reorder(source,-score), x= score, fill = Factor_agreement))+
  geom_col(color= "black", width= 0.8, position = "dodge")+
  #scale_fill_gradient2(low= "blue", high = "red")+
  scale_fill_manual(values = c("grey", "darkgreen"))+
  ggtitle("Fib")+
  theme_cowplot()+
  labs(y= "", x="ES score")+
  coord_flip()+
  theme(axis.text.x = element_text(angle= 90, hjust= 1),
        legend.position = "none")
p.lvad_Fib

p.lvad_Endo <- 
  df %>% 
  filter(condition %in% c("Endo"), source %in% Endo_lvad)%>%
  #filter(condition %in% c("Fib"), source %in% pways$source)%>%
  filter(coef== "pre_r_v_pre_nr", net!= "naba") %>%
  ggplot(aes(y= reorder(source,-score), x= score, fill = Factor_agreement))+
  geom_col(color= "black", width= 0.8, position = "dodge")+
  #scale_fill_gradient2(low= "blue", high = "red")+
  scale_fill_manual(values = c("grey", "darkgreen"))+
  ggtitle("Endo")+
  theme_cowplot()+
  labs(y= "", x="ES score")+
  coord_flip()+
  theme(axis.text.x = element_text(angle= 90, hjust= 1),
        legend.position = "none")
p.lvad_Endo

p.lvad_myeloid <- 
  df %>% 
  filter(condition %in% c("Myeloid"), source %in% Myeloid_lvad)%>%
  #filter(condition %in% c("Fib"), source %in% pways$source)%>%
  filter(coef== "pre_r_v_pre_nr", net!= "naba") %>%
  ggplot(aes(y= reorder(source,-score), x= score, fill = Factor_agreement))+
  geom_col(color= "black", width= 0.8, position = "dodge")+
  #scale_fill_gradient2(low= "blue", high = "red")+
  scale_fill_manual(values = c("grey", "darkgreen"))+
  ggtitle("Myeloid")+
  theme_cowplot()+
  labs(y= "", x="ES score")+
  coord_flip()+
  theme(axis.text.x = element_text(angle= 90, hjust= 1),
        legend.position = "none")
p.lvad_myeloid

    
  
  ggtitle("Myeloid")+
  labs(y= "", x="ES score")+
  coord_flip()+
  theme(axis.text.x = element_text(angle= 90, hjust= 1),
        legend.position = "none")

pdf("output/figures/fetal_fib_test.pdf", 
    width= 2, height= 6)
p.fetal_CM
p.fetal_Fib
dev.off()

pdf("output/figures/lvad_fib_test.pdf", 
    width= 4, height= 7)
plot_grid(p.lvad_CM,p.lvad_Fib+labs(x= ""), ncol = 2, align = "h")
dev.off()

?plot_grid()
pdf("output/figures/fetal_lvad_deg_correation.pdf",
    width= 12, 
    height= 3)
p.FET
p.LVAD
dev.off()

pdf("output/figures/fetal_lvad_deg_correation_reduced.pdf",
    width= 4, 
    height= 2.5)
p.FET_red
p.LVAD_red
dev.off()

  res.M<- df %>% ungroup()%>%
  mutate(p.adj= p.adjust(p_value))%>%
  filter(source %in% pways,
         Factor_agreement==TRUE)%>%
  select(-p_value, -p.adj,-statistic, -net, -Factor_agreement)%>%
  arrange(source)%>%
  pivot_wider(names_from = condition, values_from= score, values_fn = mean)%>%
  as.data.frame()%>% 
  column_to_rownames("source")%>%
  as.matrix()
ComplexHeatmap::Heatmap(res.M)  

res.M<- df %>% ungroup()%>%
  mutate(p.adj= p.adjust(p_value))%>%
  filter(source %in% pways,
         Factor_agreement==FALSE)%>%
  select(-p_value, -p.adj,-statistic, -net, -Factor_agreement)%>%
  arrange(source)%>%
  pivot_wider(names_from = condition, values_from= score, values_fn = mean)%>%
  as.data.frame()%>% 
  column_to_rownames("source")%>%
  as.matrix()
ComplexHeatmap::Heatmap(res.M)  


## compare the fibroblast enrichemnts of fetal and dcm

names(net_list)
metabo<-c(#   "G2M_CHECKPOINT",
  "HALLMARK_GLYCOLYSIS", 
  "HALLMARK_FATTY_ACID_METABOLISM",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
)
