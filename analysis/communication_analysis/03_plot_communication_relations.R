# HEADER --------------------------------------------
#
# Author:     Jan D. Lanzer
# Copyright   Copyright 2024 - Jan D. Lanzer
# Email:      lanzerjan@gmail.com
#
# Date:     2024-07-17
#
# Script Name:    
#
# Script Description:
# plot 

library(tidyverse)
library(cowplot)
library(igraph)
library(ggrepel)
library(RColorBrewer)
library(ggpubr)
library(igraph)
library(circlize)
library(ComplexHeatmap)

#results from running liana lr pairs w nn 
df.nn <- readRDS("output/communication/nn_meta_mofa_fact_liana.rds")
df.lr= readRDS("output/communication/lr_meta_mofa.rds")
g.loads= read_csv("output/mofa/gene_loadings.csv")
df.depend <- read.csv("output/communication/info_comms_net_MCP1.csv")

# calculate sums of interactions of LR pairs
df.graph <- df.lr%>% 
  select( receiver, sender,everything())%>%
  group_by(sender, receiver, group, Factor)%>%
  summarise(n= n(), 
            w= sum(interact, na.rm = T))%>%
  group_by(group,)%>% 
  mutate(scaled_weight= as.numeric(scale(w)))

df.graph<- df.graph %>%
  rename(predictor = sender, 
         target = receiver)%>% 
  filter(Factor=="Factor1")%>% 
  select(-w, -scaled_weight, -Factor)%>%
  pivot_wider(names_from = group, values_from= n)%>%
  mutate(NF= ifelse(is.na(NF), 0, NF),
         HF= ifelse(is.na(HF), 0, HF))%>%
  rename(lr_HF= HF, lr_NF = NF)

# join dependendencies with LR sums
df.depend <- df.depend%>% 
  left_join(df.graph, by= c("predictor", "target"))%>%
  select(-com_HF, -com_NF)%>%# remove old numbers
  mutate(lr_NF = ifelse(is.na(lr_NF),0, lr_NF))%>%
  mutate(lr_HF = ifelse(is.na(lr_HF),0, lr_HF))

# plottiing of info net ---------------------------------------------------

# main figure
#hmaps of dependencies
df.depend %>% 
  ggplot(aes(x= predictor, y= target, fill  =HF ))+
  geom_tile(color="black", size= 0.5)+
  scale_fill_gradient2(low="blue", mid="white", high ="red", na.value = "darkgrey")+
  theme_cowplot()+
  coord_equal()+
  theme(axis.text.x = element_text(angle= 90, hjust= 1), 
        axis.line = element_blank(), 
        axis.ticks = element_blank())
  
#hmapso of LR
df.depend %>% 
  ggplot(aes(x= predictor, y= target, fill  = lr_HF ))+
  geom_tile(color="black", size= 0.5)+
  scale_fill_gradient2(low="blue", mid="white", high ="red", na.value = "darkgrey")+
  theme_cowplot()+
  coord_equal()+
  theme(axis.text.x = element_text(angle= 90, hjust= 1), 
        axis.line = element_blank(), 
        axis.ticks = element_blank())+
  labs(x= "sender", y= "receiver", fill = "n Lig-Rec\ninteractions")

# MAIN Figure
# plot correlation between dependency and LR interactions
col_list<- readRDS("color_list_figures.rds")

p.corr.HF <- df.depend%>% 
  #mutate(label = ifelse(target =="CM", predictor, ""))%>%
  mutate(label = ifelse(predictor =="Fib", "Fib", ""))%>%
  #mutate(label = ifelse(target =="Fib", predictor, ""))%>%
  ggplot(aes(x= HF, y= lr_HF))+
  geom_smooth(method = "lm", color="darkgrey",se = F)+
  geom_point(aes(color= target), size= 2.4)+
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, method = "spearman")+
  geom_text_repel(aes(label = label), arrow = arrow(length=unit(0.01, "npc")), max.overlaps = 40, min.segment.length = 0.01,force = 20, force_pull = 20, nudge_x = 0.1)+
  scale_color_manual(values = col_list$ctype_colors)+ 
  #scale_color_brewer(type = "qual", palette = 3)+
  theme_cowplot()+
  labs(x= "predictive importance", 
       y= "n(L-R pairs)", 
       color="Target cell type")
p.corr.HF
p.corr.NF <- 
  df.depend%>% 
  mutate(label = ifelse(target =="Myeloid", predictor, ""))%>%
  ggplot(aes(x= NF, y= lr_NF))+
  geom_smooth(method = "lm", color="darkgrey",se = F)+
  geom_point(aes(color= target), size= 2.4)+
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, method = "spearman")+
  geom_text_repel(aes(label = label), arrow = arrow(length=unit(0.01, "npc")), max.overlaps = 40, min.segment.length = 0.01,force = 20, force_pull = 20, nudge_x = 0.1)+
  scale_color_manual(values = col_list$ctype_colors)+ 
  #scale_color_brewer(type = "qual", palette = 3)+
  theme_cowplot()+
  labs(x= "predictive importance", 
       y= "n(L-R pairs)", 
       color="Target cell type")

p.corr.predictive.importances<- df.depend%>% 
  mutate(label = ifelse(target =="CM", predictor, ""))%>%
  ggplot(aes(x= HF, y= NF))+
  geom_smooth(method = "lm", color="darkgrey",se = F)+
  geom_point(aes(color= target), size= 2.4)+
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, method = "spearman")+
  geom_text_repel(aes(label = label), arrow = arrow(length=unit(0.01, "npc")), max.overlaps = 40, min.segment.length = 0.01,force = 20, force_pull = 20, nudge_x = 0.1)+
  scale_color_manual(values = col_list$ctype_colors)+ 
  #scale_color_brewer(type = "qual", palette = 3)+
  theme_cowplot()+
  labs(x= "predictive importance in HF", 
       y= "predictive importance in NF", 
       color="Target cell type")

pdf("output/figures/comm_correlations.pdf", 
    width= 4.5, height=3)
  p.corr.HF#+coord_equal()
  p.corr.NF#+coord_equal()
  p.corr.predictive.importances+coord_equal()
dev.off()

hist(df.depend$HF, breaks = 20)

g1= graph_from_data_frame(df.depend%>% 
                            filter((HF)>0.2)%>%
                            select(predictor, 
                                   target, 
                                   HF),
                          directed = T)
#l = layout_in_circle(g1)
l= layout_nicely(g1)

col_fun = colorRamp2(c(min(E(g1)$HF), max(E(g1)$HF)),
                     c(col_list$general_up_down_colors[[1]],
                       col_list$general_up_down_colors[[2]] )
                     )

lgd = Legend(col_fun = col_fun, title = "predictive\nimportance")
E(g1)$color<- col_fun(E(g1)$HF)

pdf("output/figures/comm_network.pdf",
    width = 6, 
    height= 6)
 plot(g1, 
     layout= l,
     main = "",
     # === vertex
     vertex.color = "lavender",#rgb(0.8,0.4,0.3,0.8),          # Node color
     vertex.frame.color = "white",                 # Node border color
     vertex.shape="circle",                        # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
     vertex.size=40,                               # Size of the node (default is 15)
     vertex.size2=NA,                              # The second size of the node (e.g. for a rectangle)
     
     # === vertex label
     #vertex.label= V(n)$name,                   # Character vector used to label the nodes
     vertex.label.color="black",
     vertex.label.family="Helvetica",                  # Font family of the label (e.g.“Times”, “Helvetica”)
     vertex.label.font=2,                          # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
     vertex.label.cex=1.3,                           # Font size (multiplication factor, device-dependent)
     vertex.label.dist=0,                          # Distance between the label and the vertex
     vertex.label.degree=0 ,                       # The position of the label in relation to the vertex (use pi)
     
     # === Edge
     #edge.color = ifelse(E(g1)$n2 > 0, "red", "darkblue"),                         # Edge color
     edge.width=10* abs(E(g1)$HF),                                 # Edge width, defaults to 1
     #edge.width=E(g1)$w *10,                                 # Edge width, defaults to 1
     edge.arrow.size=1.4* abs(E(g1)$HF),                            # Arrow size, defaults to 1
     edge.arrow.width=1.15,
     #edge.arrow.mode=1,
     #edge.arrow.color= "black", # Arrow width, defaults to 1
     edge.lty="solid",                             # Line type, could be 0 or “blank”, 1 or “solid”, 2 or “dashed”, 3 or “dotted”, 4 or “dotdash”, 5 or “longdash”, 6 or “twodash”
     edge.curved=0.3                            # Edge curvature, range 0-1 (FALSE sets it to 0, TRUE to 0.5)
)
draw(lgd, x = unit(11.4, "cm"), y = unit(12, "cm"), just = c("right", "top"))
dev.off()

# join  NN results with LR and identify top interactions ------------------

# we disregard the target weights here and only join on the ligand 
# to get a df that for NF and HF and for each cell type pair 
# has all LR interactions, a median interactions score for each ligand
# and the nichenet metrics for the ligand 

df.comm <- df.nn %>% 
  mutate(group= toupper(group), 
         group= str_replace_all(group, "CT", "NF"))%>%
  select(-target, -weight, -expl_genes)%>%
  rename(source_genesymbol = test_ligand) %>%
  left_join(df.lr %>% filter(Factor=="Factor1"))%>%
  group_by(source_genesymbol, sender,receiver, Factor, group)%>% 
  mutate(med.interact=median(interact, na.rm = T))%>%
  distinct()

df.comm<- df.comm%>% select(pair, 
                  sender, 
                  receiver, 
                  source_genesymbol, 
                  target_genesymbol,
                  value_ligand,
                  value_receptor,
                  interact,med.interact,
                  everything())%>% 
  filter(!is.na(interact))

## explore which ligands are expressed in multiple cell types

x<- df.comm %>% 
  ungroup()%>%
  select(sender, source_genesymbol)%>%
  distinct()%>%
  group_by(source_genesymbol)%>%
  summarise(n.cell = n())%>%
  arrange(desc(n.cell))


df.comm<- df.comm %>% left_join(x)

# save table for a collaborator

df.comm%>%ungroup()%>%
  select("sender","receiver","source_genesymbol","target_genesymbol","aupr_corrected","group","n.cell")%>% 
  write.csv("output/all_communication_events_reheat2.csv")

df.comm%>%ungroup()%>%
  select("sender","receiver","source_genesymbol","target_genesymbol","aupr_corrected","group","n.cell")%>% 
  filter(n.cell== 1, 
         sender =="Fib", receiver =="CM")%>%
  write.csv("output/fib_to_cm_communication_events_reheat2.csv")

df.comm%>% ungroup() %>%
  select("sender","receiver","source_genesymbol","target_genesymbol","aupr_corrected","group","n.cell")%>% 
  filter(n.cell== 1, 
         sender =="Fib", receiver =="Myeloid")%>%
  write.csv("output/fib_to_myeloid_communication_events_reheat2.csv")


# explore specific cell type pairs ----------------------------------------

sender_oi= "Fib"
receiver_oi = "Myeloid"

# function to plot return a filtered data frame for a cell pair of interest
# and plot 3 metrics for their importance
get_top_hits <- function(sender_oi ,receiver_oi){
  df<-df.comm %>% 
    filter(receiver == receiver_oi,
           Factor=="Factor1",
           sender==sender_oi,
           group =="HF")%>% 
    arrange(aupr_corrected)
  p<- df%>%
    distinct(aupr_corrected, med.interact, pearson, source_genesymbol, n.cell)%>%
    ggplot(aes(x= aupr_corrected, y= -med.interact, label = source_genesymbol))+
    ggrepel::geom_label_repel(col= "black")+
    geom_point(size =5.2, color="black" )+
    geom_point(aes(color=factor(n.cell)), size =5)+
    theme_half_open()+
    scale_color_brewer(type= "seq", palette = 1)+
    labs(x="Nichenet regulatory potential",
         y= "L-R score", 
         color= "n(cell types)\nexpressing\nligand")
  print(p)
  
  return(df)
  
}

# plot target genes for selected ligands and cell type pairs --------------

# add plots with target genes 
x<-get_top_hits("Fib", "Myeloid")

y <- get_top_hits("Fib", "CM")

plot_targets <- function(sender_oi ,receiver_oi, ligands_oi,
                         target_gene_cutoff =-0.1){
 
  f_c = df.nn%>% 
    filter(receiver == receiver_oi,
           Factor=="Factor1",
           sender==sender_oi,
           group =="hf", 
           test_ligand %in% ligands_oi)
  
  ## plot the target matrix
  p1 <- f_c%>%
    #rename(weight= pearson)%>%
    ggplot(aes(y= reorder(target,weight), 
               x= reorder(test_ligand,weight),  fill = weight))+
    geom_tile()+
    #facet_grid(~sender, scales = "free")+
    scale_fill_gradient(low = "blue", high= "red")+
    theme_cowplot()+
    theme(axis.text.x = element_text(size= 9),
          axis.text.y= element_blank(), 
          axis.ticks.y =element_blank(),
          axis.line = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)
          )+
    labs(y= paste0(receiver_oi, " target genes"), 
         x= paste0(sender_oi, " ligands"))+
    scale_x_discrete(position = "top") 
  p1
  
  target.genes.oi <- g.loads%>% 
    filter(Factor=="Factor1",
           ctype == receiver_oi,
           (value) < target_gene_cutoff)%>% 
    pull(feature)
  ## plot with target gene names
  p2 = 
    f_c %>% 
    filter(target %in% target.genes.oi)%>%
    ggplot(aes(x= reorder(target,weight), y= reorder(test_ligand,weight),  fill = weight))+
    geom_tile()+
    #facet_grid(~sender, scales = "free")+
    scale_fill_gradient(low = "grey", high= col_list$general_up_down_colors[[2]])+
    theme_cowplot()+
    theme(axis.text.x = element_text(size= 10,  angle=90, hjust= 1, vjust= 0.5),
          axis.text.y= element_text(size=10),
          axis.line= element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)
          )+
    labs(x= paste0(receiver_oi, " target genes"), 
         y= paste0(sender_oi, " ligands"), 
         fill ="Regulatory\npotential")+
    coord_equal()
  p2
  list(p1, p2)
}

plot_targets("Fib", "Myeloid", ligands_oi = x %>% 
               ungroup()%>%
               filter(n.cell==1)%>%
               slice_max(order_by = aupr_corrected, n = 10) %>%
               pull(source_genesymbol)
             )


p.fib.cm.targets<- plot_targets("Fib", "CM", y %>% ungroup()%>%
               filter(n.cell==1)%>%
               slice_max(order_by = aupr_corrected, n = 10) %>%
               pull(source_genesymbol), -0.25)


p.fib.cm.targets

ligs <- c("MXRA5", "NRG1", "BMP4", "LAMA4", "SLIT1", "SLIT2", "COL1A1")

p.fib.cm.targets<- plot_targets(sender_oi = "Fib",receiver_oi = "CM", ligands_oi= ligs, 
             target_gene_cutoff = -0.3)

pdf("output/figures/comm_nn_targets_cm_fib.pdf", 
    width=5, height=2.5)
p.fib.cm.targets
dev.off()

