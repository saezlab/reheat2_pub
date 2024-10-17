# HEADER --------------------------------------------
#
# Author:     Jan D. Lanzer
# Copyright   Copyright 2024 - Jan D. Lanzer
# Email:      lanzerjan@gmail.com
#
# Date:     2024-04-15
#
# Script Name:    ~/R-projects/reheat2_pilot/deconvolute_compare_bulk_sc_comps.R
#
# Script Description:
# Try to see consistency of composition values (not composition changes)


library(tidyverse)
library(cowplot)

#load sc comps
sc<- read_csv("data/sample_comps.csv")
cm<- read_csv("data/comps_meta.csv")

# load reheat comps
bulk <- read_csv("output/deconvo/reheat_devonvoluted_composition.csv")
bulk_meta<- read_csv("output/deconvo/reheat_devonvoluted_meta.csv")


#integrate info 


composition_df = rbind(sc, bulk) %>% 
  as.data.frame()%>% 
  column_to_rownames("sample_id") %>%
  compositions::acomp() %>%
  compositions::clr()

composition_meta = rbind(bulk_meta%>% mutate(tech ="bulk"), 
                         cm %>% 
                           select(sample_id, heart_failure, study)%>%
                           rename(HeartFailure = heart_failure)%>%
                           mutate(HeartFailure = ifelse(HeartFailure == "HF", "yes", "no"))%>%
                           mutate(tech = "sc"))


saveRDS(list("composition_matrix"= composition_df, 
             "sample_meta"= composition_meta), 
        "output/composition_df.rds")

# PCA ---------------------------------------------------------------------

PCA  = prcomp(composition_df, scale= T, center = T)

pca.df <- PCA$x %>%
  as.data.frame() %>% 
  rownames_to_column("sample_id")%>%
  left_join(composition_meta, by= "sample_id")

pca.all.samps<- ggplot(pca.df, aes(x= PC1, y= PC2, shape = tech, color= HeartFailure))+
geom_point(size= 2)+
  labs(x= paste0("PC1 (",as.character(round(PCA$sdev[1]^2/sum(PCA$sdev^2)*100)),"%)"),
       y= paste("PC2 (",as.character(round(PCA$sdev[2]^2/sum(PCA$sdev^2)*100)),"%)"),
       shape= "Data modality",
       color = "Heart failure")+
  theme_cowplot()+
  scale_color_manual(values=rev(c("black", "darkgrey")))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), 
        axis.line = element_blank())+
  coord_equal()
pca.all.samps

pdf("output/figures/deconv_sc_vs_bulk_PCA.pdf", 
    height= 5, width= 5)
  pca.all.samps
dev.off()

pca.all.samps<- ggplot(pca.df, aes(x= PC5, y= PC6, shape = tech, color= HeartFailure))+
  geom_point(size= 2)+
  labs(x= paste0("PC1 (",as.character(round(PCA$sdev[1]^2/sum(PCA$sdev^2)*100)),"%)"),
       y= paste("PC2 (",as.character(round(PCA$sdev[2]^2/sum(PCA$sdev^2)*100)),"%)"))+
  theme_cowplot()+
  scale_color_manual(values=rev(c("black", "darkgrey")))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), 
        axis.line = element_blank())+
  coord_equal()
pca.all.samps
#test for associated variance of hf, tech and study

pca.df2= pca.df %>% 
  as_tibble() %>% 
  pivot_longer(cols= c(starts_with("PC")), names_to = "PC_name", values_to = "PC_value")%>%
  group_by(PC_name )%>%
  nest()

#add expl var
expl_var<- unlist(map(1:7, function(x){PCA$sdev[x]^2/sum(PCA$sdev^2)*100}))
pca.df2$expl_var<- expl_var

## hf asscoiation
hf<-pca.df2%>%
  mutate(p_val = map(data, function(dat){
    fit= lm(PC_value ~ HeartFailure, data= dat)
    # t.test(dat$PC_value, alternative = "greater") %>%
    #   broom::tidy()
    # 
    summary(fit)$coefficients[2,4]  %>% unlist()
  })) %>% 
  unnest(p_val)
  
hf_sum <- hf %>% filter(p_val< 0.001) %>% pull(expl_var) %>% sum()

##tech association
tech<-pca.df2%>%
  mutate(p_val = map(data, function(dat){
    fit= lm(PC_value ~ tech, data= dat)
    # t.test(dat$PC_value, alternative = "greater") %>%
    #   broom::tidy()
    # 
    summary(fit)$coefficients[2,4]  %>% unlist()
  })) %>% 
  unnest(p_val)

tech_sum <- tech %>% filter(p_val< 0.001) %>% pull(expl_var) %>% sum()

## study 

study<-pca.df2%>%
  mutate(p_val = map(data, function(dat){
    fit= lm(PC_value ~ study, data= dat)
    # t.test(dat$PC_value, alternative = "greater") %>%
    #   broom::tidy()
    # 
    summary(fit)$coefficients[2,4]  %>% unlist()
  })) %>% 
  unnest(p_val)

study_sum <- study %>% filter(p_val< 0.001) %>% pull(expl_var) %>% sum()
hf_sum
study_sum
tech_sum
  
test.res= pca.df %>% 
  as_tibble() %>% 
  pivot_longer(cols= c(starts_with("PC")), names_to = "PC_name", values_to = "PC_value")%>%
  group_by(PC_name )%>%
  nest()%>%
  mutate(test_res = map(data, function(dat){
    fit= lm(PC_value ~ tech, data= dat)
    # t.test(dat$PC_value, alternative = "greater") %>%
    #   broom::tidy()
    # 
    summary(fit)$coefficients[2,4]  %>% unlist()
  })) %>% 
  unnest(test_res)


test.res= pca.df %>% 
  as_tibble() %>% 
  pivot_longer(cols= c(starts_with("PC")), names_to = "PC_name", values_to = "PC_value")%>%
  group_by(PC_name )%>%
  nest()%>%
  mutate(test_res = map(data, function(dat){
    fit= lm(PC_value ~ study, data= dat)
    # t.test(dat$PC_value, alternative = "greater") %>%
    #   broom::tidy()
    # 
    summary(fit)$coefficients[2,4]  %>% unlist()
  })) %>% 
  unnest(test_res)


## only NF samples
samps<- composition_meta%>% filter(HeartFailure=="no")%>% pull(sample_id)
samps = samps[samps != "hs_lv_079_C"]
PCA  = prcomp(composition_df[samps, ], scale= F, center = T)
pca.df <- PCA$x %>%as.data.frame() %>% rownames_to_column("sample_id")%>%
  left_join(composition_meta, by= "sample_id")

ggplot(pca.df, aes(x= PC1, y= PC2, shape = tech, color= study))+
  geom_point()+
  labs(x= paste0("PC1 (",as.character(round(PCA$sdev[1]^2/sum(PCA$sdev^2)*100)),"%)"),
       y= paste("PC2 (",as.character(round(PCA$sdev[2]^2/sum(PCA$sdev^2)*100)),"%)"))

"hs_lv_079_C"
pca.df%>% filter(PC2< -0.5)



# now we see if the tech effect gets mitigated when comparing comp --------

study_diff_stats <- composition_df %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(-sample_id, names_to = "cell_type") %>%
  left_join(composition_meta, by = "sample_id") %>%
  group_by( cell_type, study) %>%
  mutate(HeartFailure = factor(HeartFailure, levels=c("yes", "no")))%>%
  nest() %>%
  dplyr::mutate(HF_diff = map(data, function(dat) {
    
    t.test(value ~ HeartFailure, data = dat) %>%
      broom::tidy()
    
  })) %>%
  dplyr::select(-data) %>%
  unnest(HF_diff)


diff_stats_df <-study_diff_stats%>% select(study, cell_type, estimate) %>%
  pivot_wider(names_from = study, values_from= estimate)%>%
  column_to_rownames("cell_type")

PCA  = prcomp(t(diff_stats_df), scale= T, center = T)

pca.df <- PCA$x %>% as.data.frame() %>% 
  rownames_to_column("sample_id")%>%
  mutate(tech= ifelse(grepl("CM", sample_id), "sc", "bulk"))

p.pca.diff<- ggplot(pca.df, aes(x= PC1, y= PC2, color = tech))+
  geom_point(size= 4)+
  labs(x= paste0("PC1 (",as.character(round(PCA$sdev[1]^2/sum(PCA$sdev^2)*100)),"%)"),
       y= paste("PC2 (",as.character(round(PCA$sdev[2]^2/sum(PCA$sdev^2)*100)),"%)"))+
  ggrepel::geom_label_repel(aes(label = sample_id))+
  theme_cowplot()

cowplot::plot_grid(pca.all.samps,  p.pca.diff)

study_diff_stats %>%
  left_join(composition_meta%>% distinct(study, tech))%>%
  ggplot(aes(x= reorder(cell_type,estimate), y= estimate , fill = tech))+
  #geom_point(alpha= 0.5)+
  geom_boxplot(alpha=1, width = 0.5)+
    geom_hline(yintercept = 0)+
    theme_cowplot()+
  scale_fill_manual(values = unname(color_list$fib_state_colors))+
  labs(x= "",
       fill ="data\nmodality",
       y= "compositional change")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), 
        axis.line = element_blank(), 
        axis.text.x= element_text(angle= 45, hjust= 1))

study_diff_stats %>%
    left_join(composition_meta%>% distinct(study, tech))%>%
    ggplot(aes(x= reorder(cell_type,estimate),  y= estimate ))+
    geom_boxplot()+
    geom_hline(yintercept = 0)+
    theme_cowplot()
