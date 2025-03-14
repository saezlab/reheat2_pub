# HEADER --------------------------------------------
#
# Author:     Jan D. Lanzer
# Copyright   Copyright 2023 - Jan D. Lanzer
# Email:      lanzerjan@gmail.com
#
# Date:     2023-10-24
#
# Script Name:    ~/R-projects/reheat2_pilot/downstream_analysis/state_marker_MCP_ORA.R
#
# Script Description:
# Perform use MCPs and perform ORA with state marker

library(tidyverse)
library(decoupleR)
library(cowplot)

source(save_source_data())

g.load= read.csv("output/mofa/gene_loadings.csv")
#meta= read.csv("output/mofa/metamodel_meta.csv")
states <- read.csv("../reheat2_cleaned_jan/data/fibroblasts/cluster_markers_processed.csv")
col_list <- readRDS("color_list_figures.rds")

# cell_state_marker_ORA ---------------------------------------------------
col_list$fib_state_colors
#prepare the state marker, select n top genes
n_top_genes= 200

net= states%>%
  group_by(cluster)%>% 
  filter(logfc>0)%>%
  arrange(pval, desc(logfc))%>%
  slice_head(n= n_top_genes)%>%
  dplyr::rename(target= "gene", 
         source= "cluster", 
         )%>%
  mutate(mor= 1)
  
# prepare the gene loading matrix:

mat= g.load %>% 
  as_tibble()%>% 
  pivot_wider(names_from = Factor, values_from = value)%>%
  group_by(ctype)%>% nest()
  
mat <- mat %>%
  dplyr::mutate(data = map(data, function(x){
    #print(head(x))
    x %>% 
      as.data.frame()%>%
      column_to_rownames("feature")%>%
      as.matrix()
  })) 

# run ulm 
res= mat %>% 
  mutate(decopler_res= map(data, function(x){
    decoupleR::run_ulm(mat = x,network = net)#,n_up = 200 )
  }))

# plot  result --------------------------------------------

# we used fibroblast states, thus we will only plot results of Fib factors:

df= res %>% 
  filter(ctype =="Fib")%>%
  unnest(decopler_res)%>% 
  mutate(p_adj=p.adjust(p_value))%>% 
  mutate(significant = ifelse(p_adj < 0.01, "*",""))
df <- df %>%
  mutate(condition = factor(condition, levels = unique(df$condition)))
df<-df %>% group_by(condition)%>% mutate(sc_score = scale(score))
p_hmap= df %>% 
  ggplot(., aes(x= source, y= condition, fill = score))+
  geom_tile(color="black") +
  scale_fill_gradient2(low = "blue", mid= "white", high = "red") +  # Adjust the color gradient as needed
  #labs(title = "State marker enriched in factor loadings") +
  geom_text(aes(label= significant))+
  theme_minimal()+
  theme_cowplot()+
  labs(x= "Fibro states", 
       y= "Multicellular Factors")

p_hmap


# compare to t-stats from composition -------------------------------------

# to see whether the enrichment score associates with composition changes we will 
# compare them to t-values from the composition test:

study_diff_stats= readRDS( file= "output/fib_sub_analysis/fib_state_composition_tvals.rds")%>%
  as_tibble()

df 
#join with ora results
df.j= study_diff_stats %>% 
  select(cell_type, batch,  statistic)%>%
  left_join(df%>% select(-statistic)%>% dplyr::rename(cell_type =source), by= c("cell_type")  )


p1= df.j %>% #dplyr::filter(condition ==  "Factor1")%>%
  group_by(cell_type, condition)%>% 
  mutate(m.composition_t= mean(statistic))%>%
  ggplot(., aes(x= statistic, y= score, shape= batch))+
  facet_grid(~condition)+
  geom_point(size= 2)+
  geom_smooth(aes(colour=batch),method = "lm")+
  theme_cowplot()+
  labs(x= "t-statistic for compositional change",
       y= "Fib marker enrichment")
p1

# plot facot
p_mean_t_score_asso= df.j %>% #dplyr::filter(condition ==  "Factor1")%>%
  group_by(cell_type, condition)%>% 
  mutate(m.composition_t= mean(statistic))%>%
  ggplot(., aes(x= m.composition_t, y= score, color = cell_type))+
  geom_point(size= 2)+
  #geom_smooth(method = "lm")+
  facet_grid(~ condition) +
  theme_cowplot()

df.j%>% 
  group_by(cell_type, condition)%>% 
  mutate(m.composition_t= mean(statistic))%>%
  group_by(condition)%>% 
  mutate(corr= cor(m.composition_t, score))%>%
  distinct(corr, condition)%>%
  ggplot(., aes(x= reorder(condition,corr),  y= corr))+
  geom_col(width= 0.05)+
  geom_point()+
  theme_cowplot()+
  geom_hline(yintercept = 0)+
  coord_flip()+
  labs(y= "Pearson correlation",
      x= "")

p2= df.j %>% #dplyr::filter(condition ==  "Factor1")%>%
  group_by(cell_type, condition)%>% 
  mutate(m.composition_t= mean(statistic))%>%
  ggplot(., aes(x= m.composition_t, y= score))+
  #geom_point(aes(shape= cell_type), size= 2)+
  geom_smooth(aes(color = condition), method = "lm", alpha= 0.2)+
  scale_color_brewer(palette = "Dark2")+
  geom_text(label= condition)+
  #facet_grid(~ condition) +
  theme_cowplot()
  
p2
df2= df.j %>% #dplyr::filter(condition ==  "Factor1")%>%
  group_by(source, condition)%>% 
  mutate(m.composition_t= mean(composition_t, na.rm =T))%>%
  distinct(score, source, condition, m.composition_t)%>%
  ungroup()%>%
  group_by(condition)%>%
  mutate(cor = cor(score, m.composition_t))

df2%>% distinct( condition, cor)%>%
  ggplot(., aes(x= cor, y= condition,fill= cor))+
  geom_col()


pdf("output/fib_sub_analysis/figures/dol_comp_and_state_ora.pdf", 
    width = 10)
  p_hmap+theme_cowplot()
  p_mean_t_score_asso
dev.off()

library(RColorBrewer)
# plot only factor1 
p_mean_t_score_asso_fact1= df.j %>%
  dplyr::filter(condition %in% c("Factor1", "Factor2"))%>%
  mutate(condition= str_replace_all(condition, "Factor", "Fib_MCP"))%>%
  group_by(cell_type, condition)%>% 
  mutate(m.composition_t= mean(statistic))%>%
  ggplot(., aes(x= m.composition_t, y= sc_score))+#, color = cell_type))+
  geom_smooth(method = "lm", color="black")+
  scale_color_manual(values = c("#FFB200", "#E4003A", "#17BECF" ,"#B60071" ,
                                "#1F77B4", "black"))+
  #scale_color_brewer(type="qual", palette = 6)+
  geom_point(size= 3, aes(color=cell_type))+
  facet_grid(~ condition, scales = "free") +
  theme_cowplot()+
  theme(panel.grid.major = element_line(color = "darkgrey",
                                        size = 0.5,
                                        linetype = 2), 
        strip.background =element_rect(fill="lavender"), 
    
        )+
  labs(y= "State marker\nenrichment", 
       x= "composition change (t-statistic)", 
       color ="")

p_mean_t_score_asso_fact1

pdf("output/figures/fib_ES_vs_comp.pdf", 
    height= 2.5,
    width = 5.5)
p_mean_t_score_asso_fact1
dev.off()

##save source data

df.j %>%
  dplyr::filter(condition %in% c("Factor1"))%>%
  mutate(condition= str_replace_all(condition, "Factor", "Fib_MCP"))%>%
  group_by(cell_type, condition)%>% 
  mutate(mean_tval_compositionchange= mean(statistic))%>%
  distinct(cell_type, statistic, score, mean_tval_compositionchange)%>%
  save_source_data(data= ., T, 4 , "D")

#get correlations for paper 

df.j%>% 
  dplyr::filter(condition %in% c("Factor1", "Factor2"))%>%
  group_by(cell_type, condition)%>% 
  mutate(m.composition_t= mean(statistic))%>%
  group_by(condition)%>% 
  mutate(corr= cor.test(m.composition_t, score)$estimate,
         corr_p= cor.test(m.composition_t, score)$p.val)%>%
  distinct(condition,corr,corr_p )
  View()

df.j %>%
  mutate(condition= str_replace_all(condition, "Factor", "MCP"))%>%
  group_by(cell_type, condition)%>% 
  mutate(m.composition_t= mean(statistic))

