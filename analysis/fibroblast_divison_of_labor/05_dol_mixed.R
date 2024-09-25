# HEADER --------------------------------------------
#
# Author:     Jan D. Lanzer
# Copyright   Copyright 2023 - Jan D. Lanzer
# Email:      lanzerjan@gmail.com
#
# Date:     2023-10-30
#
# Script Name:    
#
# Script Description:
# divison of labor via linear mixed model

library(tidyverse)
library(cowplot)
library(edgeR)
library(SummarizedExperiment)
library(lme4)
library(lmerTest)
library("r2glmm")
library(tune)
library(ggpubr)
library(ggrepel)
library(decoupleR)

col_list <- readRDS("color_list_figures.rds")
g.load= read.csv("output/mofa/gene_loadings.csv") %>% as_tibble()
pb = readRDS( "output/fib_sub_analysis/all_state_pb_processed_sub_reg.rds")

# prep data  -------------------------------------------------------------

## select features from the factor loadings of fibro view
features_fac1= g.load %>%
  filter(Factor %in% c("Factor1", "Factor2"),
         ctype == "Fib")%>% 
  filter( abs(value)>0.1)%>%
  pull(feature)%>% 
  unique()

length(features_fac1)

pb_long = pb %>% do.call(rbind, .)%>%
  mutate(cell_state= as.factor(cell_state),
         heart_failure= as.factor(heart_failure),
         sample_id =as.factor(sample_id))
#filter for low cell counts
pb_long<- pb_long%>% group_by(sample_state)%>%
  filter(cell_counts>40)

features_fac1= features_fac1[features_fac1 %in% pb_long$feature]
length(features_fac1)

# test to run lmer for one gene -------------------------------------------

# function to plot normalized expression value to demonstrate
# individual gene's expression pattern:

plot_exp= function(genes, ...){
  pls<- lapply(genes, function(gene){
    
    mod.df= pb_long %>% filter(feature ==gene) %>% 
      mutate(cell_state= paste0("Fib_", cell_state))
    
    p.exp= ggplot(mod.df, aes(x= cell_state, 
                              y= (value),
                             
                              fill = heart_failure))+
      geom_boxplot(outlier.shape = NA,
                   aes(color = cell_state),
                   show.legend = T)+
      scale_fill_manual(values= c("white", "darkgrey"))+
      scale_color_manual(values= c(unname(col_list$fib_state_colors), "black"))+
      labs(y= paste0(gene), 
           x= "", fill ="")+
      theme_cowplot()+
      scale_y_continuous(limits = quantile(mod.df$value, c(0.05, 0.95)))+
      theme(axis.text.x = element_text(angle =90, hjust= 1, vjust= 0.5))+
      guides(color = "none" )
  }) 
  
  legend <- get_legend(pls[[1]])
  pls<-lapply(pls, function(x){x+theme(legend.position="none")})
  pls$leg= legend
  plot_grid(plotlist =pls, ...)
  
}

p.showcase_exp <- plot_exp(genes = c("CDK8", "FGF14", "COL24A1"), ncol = 1, align = "v")
p.showcase_exp

plot_qq= function(x){
  plot(x) ## TA-plot
  par(mfrow = c(1, 2))
  qqnorm(ranef(x)$batch[,"(Intercept)"], 
         main = "Random effects")
  qqnorm(resid(x), main = "Residuals")
  
}


# lmm to assign gene programs --------------------------------------------------------------

#fit the same linear mixed model modelling cell labels
excluded_feats<- character(0)

res= lapply(features_fac1, function(gene){
  print(gene)
  
  mod.df= pb_long %>% filter(feature == gene)
  
  #check whether there are samples across all 4 batches
  x= table(mod.df$cell_state, mod.df$heart_failure, mod.df$batch)
  
  length(attr(x, "dimnames")[[3]])
  
  if(length(attr(x, "dimnames")[[3]]) <4 ){
    excluded_feats<<-c(excluded_feats,gene)
    #print(paste0(x, " skipped"))
    return(NULL)
  }
  
  fit.clust = lmer(value ~ 1 + cell_state + heart_failure + (1| batch:sample_id) , data = mod.df)

  #extract random effect variance 
  a =summary(fit.clust)$varcor$batch[1] # random effect batch variance 
  b = attr(summary(fit.clust)$varcor, "sc")^2  # random effect residual variance 
  #calc explained var by random effect (batch)
  r.eff.var= a/(a+b)
  
  # use this to get partial R2 values
  df= r2beta(fit.clust)%>% as.data.frame()

  df= df %>% mutate(r.eff.var= r.eff.var, 
                feature = gene)
  
  return(df)
})%>% do.call(rbind, .)

res  

# plot results ------------------------------------------------------------

## extract R squared
df.plot= res %>% select(Effect, Rsq, r.eff.var, feature)%>%
  pivot_wider(names_from = Effect, values_from= Rsq)%>%
  mutate(ratio= heart_failure/cell_state, 
         logratio= log10(ratio))%>%
  left_join(g.load %>% filter(ctype=="Fib", Factor== "Factor1"), by= "feature")

#define threshold for R2 
thresh=0.1

df.plot<- df.plot %>%
  mutate(sign = ifelse(value>0, "NF", "HF"),
         #lab = ifelse(feature %in% label.genes, feature, ""),
         group= ifelse(heart_failure <thresh & cell_state> thresh, "Compositional","Undefined"),
         group= ifelse(heart_failure >thresh & cell_state< thresh, "Molecular",group),
         group= ifelse(heart_failure >thresh & cell_state> thresh, "Comp/Mol",group)) %>%
  dplyr::rename("Random_eff_var"= r.eff.var,
                "R2_model"= Model,
                "R2_HF"= heart_failure, 
                "R2_CS"= cell_state, 
                "HF_CS_ratio"= ratio
  )


saveRDS(df.plot,"output/fib_sub_analysis/fib_state_dol_mixed_results_reg.rds")

df.plot= readRDS("output/fib_sub_analysis/fib_state_dol_mixed_results_reg.rds")

label.genes= c("CDK8", "SCN7A", #comp hf
               "GSN", "FBN1", # comp nf
               "COL24A1", "PLCE1", "FKBP5", #trans hf
               "SPOCK1", "MALAT1",  #trans nf
               "POSTN","FAP", "CTHRC1", "THBS4","FGF14",  #comp/trans hf
               "DCN", "CFD") # comp mol nf


df.plot%>% arrange(desc(R2_HF))%>% filter(R2_CS<0.1)
df.plot%>% arrange(desc(R2_CS))%>% filter(R2_HF<0.1)
df.plot%>% arrange(desc(R2_CS))%>% filter(R2_HF>0.4, R2_CS > 0.4)

p.r2= df.plot %>%
  mutate(sign = ifelse(value>0, "Nonfailing", "Heart failure"),
         lab = ifelse(feature %in% label.genes, feature, ""))%>%
  filter(abs(value)>0.1)%>%
  ggplot(., aes(x= R2_HF, y= R2_CS, color = group))+
  geom_abline(slope = 1)+
  facet_wrap(~sign)+
  geom_point()+
  scale_color_manual(values=c( "#E41A1C", "#377EB8", "#984EA3","#4DAF4A"))+
  theme_cowplot()+
  theme(strip.background =element_rect(fill="lavender"))+
  ggrepel::geom_label_repel(aes(label = lab), 
                            #color="black",
                            alpha= 0.7,
                             show.legend = F,
                             max.overlaps = 100, 
                             size= 4)+
  coord_obs_pred()+
  labs(x= "partial R² heart failure", 
       y= "parial R² cell state", 
       col= "Expression\ngroup")
p.r2

pdf("output/figures/dol_genes_overview.pdf", 
    width= 8, 
    height= 5)
p.r2
dev.off()

## add the expression plot to show case the labeled genes

p.showcase_exp <- plot_exp(genes = c("CDK8", "FGF14", "COL24A1"), ncol = 2)
p.showcase_exp

pdf("output/figures/fib_expression_showcase_dol.pdf", 
    width = 5, height = 4)
p.showcase_exp
dev.off()

genes = c("CDK8", "FGF14", "COL24A1")
mod.df= pb_long %>% filter(feature %in% genes) %>% 
  mutate(cell_state= paste0("Fib_", cell_state))

p.showcase<- ggplot(mod.df, aes(x= cell_state, 
                          y= (value),
                          
                          fill = heart_failure))+
  geom_boxplot(outlier.shape = NA,
               aes(color = cell_state),
               show.legend = T,
               width= 0.8, size= 0.7)+
  scale_fill_manual(values= c("white", "darkgrey"))+
  facet_wrap(~feature, ncol =1, scales="free_y")+
  scale_color_manual(values= c(unname(col_list$fib_state_colors), "black"))+
  labs(y= "", 
       x= "", fill ="")+
  theme_cowplot()+
  scale_y_continuous(limits = quantile(mod.df$value, c(0.05, 0.95)))+
  theme(axis.text.x = element_text(angle =90, hjust= 1, vjust= 0.5),
        strip.background = element_blank())+
  guides(color = "none" )
pdf("output/figures/fib_expression_showcase_dol2.pdf", 
    width = 3.2, height = 5)
p.showcase
dev.off()


  p.mof= df.plot %>%
   ggplot(., aes(x= R2_model, y= abs(value)))+
   geom_abline(slope = 1)+
   geom_point()+
   scale_color_gradient2(low="blue", mid= "grey", high = "red" )+
   theme_cowplot()+
   #geom_smooth()+
  # stat_cor()+
   labs(x= "R2, linear mixed model", y= "loading of MetaMofacell model")
   #coord_obs_pred()
 
 
p.mof= df.plot %>%
   ggplot(., aes(x= Model, y= abs(value), color = r.eff.var))+
   geom_abline(slope = 1)+
   geom_point()+
   scale_color_gradient2(low="blue", mid= "grey", high = "red" )+
   theme_cowplot()+
   #geom_smooth()+
   stat_cor()+
   labs(x= "R2, linear mixed model", y= "loading of MetaMofacell model")+
   coord_obs_pred()

p.overview= plot_grid(p.r2, plot_grid(p.mof, p.r.eff), ncol = 1)

pdf("output/fib_sub_analysis/figures/dol_exp_of_top_mixed_model_overview.pdf",
    width= 20, 
    height= 12)
p.overview
dev.off()

top_hf_genes= df.plot%>% arrange(desc(R2_HF))%>% 
  filter(value>0) %>%
  slice_head(n= 10)%>% pull(feature)
top_state_genes= df.plot%>% arrange(desc(R2_CS))%>% 
  filter(value>0) %>% 
  slice_head(n= 10)%>% pull(feature)


pdf("output/fib_sub_analysis/figures/dol_exp_of_top_mixed_model_genes.pdf",
    width= 20, 
    height= 12)
  plot_exp(top_hf_genes)+theme_cowplot()
  plot_exp(top_state_genes)
dev.off()


# select marker -----------------------------------------------------------


df.plot%>% arrange(desc(R2_CS))%>% filter(R2_HF<0.1, sign=="HF")
df.plot%>% arrange((value))%>% filter(group=="Molecular",  sign=="HF")
df.plot%>% arrange((value))%>% filter(group=="Molecular",  sign=="HF")%>% print(n=50)
df.plot%>% arrange(desc(R2_HF))%>% filter(group=="Molecular",  sign=="HF")%>% pull(feature)
df.plot%>% filter(grepl("NAB1", feature))


plot_exp(c("NAB2", "NAB1"))
plot_exp(genes = c("MXRA5"))


# compare mofa loadings between groups ------------------------------------

p.box.mcp.dol<- df.plot %>% 
  filter(ctype=="Fib",
         Factor=="Factor1")%>%
  drop_na()%>%
  mutate(sign = ifelse(sign=="HF", "Heart failure", "Nonfailing"))%>%
  #filter(abs(value)< 0.1)%>%
  ggplot(aes(x=group , y= abs(value), color= group))+
  geom_hline(yintercept = 0)+
  facet_grid(~sign)+
  geom_boxplot()+
  stat_compare_means(comparisons = list(c("Comp/Mol", "Molecular"), 
                                        c("Molecular", "Compositional")))+
  theme_cowplot()+
  scale_color_manual(values=c( "#E41A1C", "#377EB8", "#984EA3","#4DAF4A"))+
  theme(axis.text.x = element_text(angle=  90,vjust=0.5,  hjust= 1), 
        legend.position = "none")+
  theme(strip.background =element_rect(fill="lavender"))+
  labs(y= "absolute MCP1 loading", x="")+
  ylim(c(0,2))
  
p.box.mcp.dol

pdf("output/figures/fib_dol_mcp_loadings.pdf", 
    width = 3.5, 
    height= 4.1)
p.box.mcp.dol
dev.off()


# compare reheat scores ---------------------------------------------------

directed_signature = readRDS("~/R-projects/Collaborations/cheerio/app_data/signature.rds")
new_sig<- read_delim("output/reheat1/meta_analysis_summary.txt")
new_sig<- new_sig %>% mutate(weight= -log10(fisher_pvalue) * sign(mean_t))%>%arrange(desc(weight))
#new_sig %>% filter(gene=="C2orf83")
#MAT= as.matrix(directed_signature %>% as.data.frame()%>% column_to_rownames("gene"))
MAT= as.matrix(new_sig %>% distinct(gene, weight) %>% as.data.frame() %>% column_to_rownames("gene"))
df.plot = df.plot %>% mutate(gex_group= paste0(group, "_", sign))

groups= unique(df.plot$gex_group)

df =lapply(groups, function(x){
  
  net= df.plot %>%
    filter(gex_group ==x)%>%
    select(gex_group, feature)%>%
    mutate(mor= 1)
  
  run_ulm(MAT, net, .source = "gex_group", 
          .target = "feature")
}) %>% do.call(rbind,.)

p.reheat_enrich_dol <- df %>% 
  separate( source, into = c("group", "heart_failure"), sep = "_")%>%
  mutate(p_value_adj = p.adjust(p_value)) %>%
  mutate(sig= p_value_adj<0.01)%>%
  mutate(heart_failure = ifelse(heart_failure=="HF", "Heart failure", "Nonfailing"))%>%
  ggplot(aes(x= group, y= score , fill= -log10(p_value_adj)))+
  geom_col(aes(col=group),col ="black",  show.legend = T)+
  geom_vline(xintercept = 0)+
  facet_wrap(~heart_failure)+
  scale_fill_gradient(low = "grey", high = col_list$general_up_down_colors[[2]])+
  #scale_color_manual(values = c("white", "black"))+
  scale_color_manual(values=c( "#E41A1C", "#377EB8", "#984EA3","#4DAF4A"))+
  theme_cowplot()+
  geom_hline(yintercept = 0)+
  theme(strip.background =element_rect(fill="lavender"))+
  theme(axis.text.x = element_text(angle= 90,vjust=0.5,  hjust= 1))+
  labs(x="", 
       fill="-log(p)",
       y= "enrichment score\nconsensus bulk signature")



pdf("output/figures/reheat_enrich_dol.pdf", 
    width = 3, 
    height= 3)
p.reheat_enrich_dol
dev.off()

p.comb <- cowplot::plot_grid(p.box.mcp.dol, 
                             p.reheat_enrich_dol,
                   rel_widths = c(1,1.3))

pdf("output/figures/fib_dol_composite.pdf", 
    width = 6.5, 
    height= 4)
p.comb
dev.off()




# compare expression of state marker with dol -----------------------------

state.m <- read_csv("output/fib_sub_analysis/cluster_markers_leiden_0.7_processed_reg.csv")

state.m$cluster %>% unique()

top.marker<- state.m %>% 
  group_by(cluster)%>% 
  filter(logfc>.5, pval<0.01)%>%
  slice_min(order_by = logfc, n = 200)

df.jacc<- sapply(unique(top.marker$cluster), function(x){
  
  sapply(unique(df.plot$gex_group), function(y){

    genes_statemarker<-top.marker %>% filter(cluster ==x)%>% pull(gene)
    genes_dol_group<- df.plot %>% filter(gex_group ==y)%>% pull(feature)
    jacc = length(intersect(genes_dol_group, genes_statemarker)) / length(union(genes_dol_group, genes_statemarker))
    })
  })

print(ComplexHeatmap::Heatmap(scale(df.jacc), name = "scaled\nJaccard's\nIndex"))
dev.off()

x<- top.marker %>% left_join(df.plot%>% distinct(feature, gex_group), join_by("gene"=="feature"))

x %>% ggplot(aes(x= cluster, fill = gex_group))+
  geom_bar(stat="count")
