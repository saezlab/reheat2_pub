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

library(grid)  # Needed for text customization
library(circlize)  # Needed for colorRamp2

source("make_source_data.R")

#results from running liana lr pairs w nn 
df.nn <- readRDS("output/communication/nn_meta_mofa_fact_liana.rds")

df.lr= readRDS("output/communication/lr_meta_mofa.rds")

df.lr %>% 
  filter(Factor =="Factor1")%>% select(sender, receiver, group, source_genesymbol)
#save filtered version for shiny app
# df.nn%>% 
#   filter( Factor=="Factor1")%>%
#   select(sender, receiver, aupr_corrected, test_ligand,target, weight, group)%>%
#   left_join(g.loads%>% filter(Factor=="Factor1"), 
#             by = c("sender" = "ctype", "target" = "feature"))%>%
#   saveRDS("output/communication/nn_target_weights_mcp1.rds")
#         


g.loads= read_csv("output/mofa/gene_loadings.csv")
df.depend <- read.csv("output/communication/info_comms_net_MCP1.csv")

 df.nn %>%
  filter(Factor == "Factor1") %>%
  mutate(group =ifelse(group =="hf", "HF", "NF")) %>%
  select(sender, receiver, aupr_corrected, test_ligand, target, weight, group) %>%
  left_join(g.loads %>% filter(Factor == "Factor1"), 
            by = c("sender" = "ctype", "target" = "feature")) %>%
  semi_join(df.lr %>% 
              filter(Factor == "Factor1") %>% 
              select(sender, receiver, group, source_genesymbol),
            by = c("sender", "receiver", "group", "test_ligand" = "source_genesymbol"))%>%
  saveRDS("output/communication/nn_target_weights_mcp1.rds")
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
df.depend %>% saveRDS(., "output/communication/df_communications_mcp1_dependencies.rds")
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

#save source data
df.depend %>%
  select(-NF, -com_NF)%>%
  save_source_data(T, fig_number = 3, 
                   panel_letter = "F",
                   data = ., 
                   bottom_description = "HF, predictive importance in HF; lr_HF, number of lr pairs in HF" )

df.depend %>%
  select(-HF, -lr_HF)%>%
  save_source_data(F, fig_number = 7, 
                   panel_letter = "D",
                   data = ., 
                   bottom_description = "NF, predictive importance in NF; lr_NF, number of lr pairs in NF" )
df.depend %>%
  select(-lr_NF, -lr_HF)%>%
  save_source_data(F, fig_number = 7, 
                   panel_letter = "B",
                   data = ., 
                   bottom_description = "NF, predictive importance in NF; HF, predictive importance in HF" )


######
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

df.depend%>% 
  filter((HF)>0.2)%>%
  select(predictor, 
         target, 
         HF)%>%
  save_source_data(T, 3, "E", data= .)

# revision update
# add nf network 


g1= graph_from_data_frame(df.depend%>% 
                            filter((NF)>0.2)%>%
                            select(predictor, 
                                   target, 
                                   NF),
                          directed = T)
#l = layout_in_circle(g1)
l= layout_nicely(g1)

col_fun = colorRamp2(c(min(E(g1)$NF), max(E(g1)$NF)),
                     c(col_list$general_up_down_colors[[1]],
                       col_list$general_up_down_colors[[2]] )
)

lgd = Legend(col_fun = col_fun, title = "predictive\nimportance")
E(g1)$color<- col_fun(E(g1)$NF)

pdf("output/figures/comm_network_NF.pdf",
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
     edge.width=10* abs(E(g1)$NF),                                 # Edge width, defaults to 1
     #edge.width=E(g1)$w *10,                                 # Edge width, defaults to 1
     edge.arrow.size=2.5* abs(E(g1)$NF),                            # Arrow size, defaults to 1
     edge.arrow.width=1.15,
     #edge.arrow.mode=1,
     #edge.arrow.color= "black", # Arrow width, defaults to 1
     edge.lty="solid",                             # Line type, could be 0 or “blank”, 1 or “solid”, 2 or “dashed”, 3 or “dotted”, 4 or “dotdash”, 5 or “longdash”, 6 or “twodash”
     edge.curved=0.3                            # Edge curvature, range 0-1 (FALSE sets it to 0, TRUE to 0.5)
)
draw(lgd, x = unit(15, "cm"), y = unit(13, "cm"), just = c("right", "top"))
dev.off()

df.depend%>% 
  filter((NF)>0.2)%>%
  select(predictor, 
         target, 
         HF)%>%
  save_source_data(F, 7, "A", data= .)

# join  NN results with LR and identify top interactions ------------------

# we disregard the target weights here and only join on the ligand 
# to get a df that for NF and HF and for each cell type pair 
# has all LR interactions, a median interactions score for each ligand
# and the nichenet metrics for the ligand 

df.comm <- df.nn %>% 
  mutate(group= toupper(group), 
         group= str_replace_all(group, "CT", "NF"))%>%
  select(-target, -weight, -expl_genes)%>%
  dplyr::rename(source_genesymbol = test_ligand) %>%
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

saveRDS(df.comm, "output/communication/df_communications_mcp1.rds")
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
  
  return(list("df" =df, "p"=p))
  
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

plot_targets("Fib", "Myeloid", ligands_oi = x$df %>% 
               ungroup()%>%
               filter(n.cell==1)%>%
               slice_max(order_by = aupr_corrected, n = 10) %>%
               pull(source_genesymbol)
             )


p.fib.cm.targets<- plot_targets("Fib", "CM", y$df %>% ungroup()%>%
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


#save source data for fib- cm target gene heatmap
f_c = df.nn%>% 
  filter(receiver == "CM",
         Factor=="Factor1",
         sender=="Fib",
         group =="hf", 
         test_ligand %in% ligs)
target.genes.oi <- g.loads%>% 
  filter(Factor=="Factor1",
         ctype == "CM",
         (value) < 0.1)%>% 
  pull(feature)
f_c %>% 
  filter(target %in% target.genes.oi)%>%
  select(test_ligand, aupr_corrected, pearson, sender, receiver, weight, group)%>%
  save_source_data(T, 3, "G", .)

# revision update, address cell types targeting fibs -----------------------

df.depend%>% 
  filter(target =="Fib")%>%
  select(predictor, 
         target, 
         NF, HF)%>%
  pivot_longer(cols= c(NF, HF), names_to = "HF_status", values_to = "predictive_importance")%>%
  ggplot(aes(x= reorder(predictor, predictive_importance), y= predictive_importance))+
  geom_col(fill ="lavender", color="black", width = 0.4)+
  theme_cowplot()+
  #coord_flip()+
  labs(y= "Predictive Importance", 
       x= "") +
  theme(plot.margin = margin(20, 20, 20, 20, "pt"), 
        axis.text.x = element_text(angle= 90, hjust= 1, vjust= 0.5)) +
  facet_grid(~HF_status)

# we will look at possible ligands
p1= get_top_hits("CM", "Fib")
p1$p+ggtitle("CM to Fib ligands") 
p2= get_top_hits("Myeloid", "Fib")


p3= plot_grid(  p2$p+
              ggtitle("Myeloid")+
              theme(legend.position = "none"), 
            p1$p+
              ggtitle("CM"), 
            rel_widths= c(1,1.4)
          
)

pdf("output/figures/comm_ligands_to_fibs.pdf",
    width= 8, height= 3.5)
p3
dev.off()
# LR network MCP1 as a heatmap --------------------------------------------


hmap1%>%
  filter(LR=="POSTN_PTK7") %>%
select(sender, receiver, LR, n_receiver, n_sender, group, n_group)

LRs= df.comm %>%
  #filter(abs(interact)>0.2) %>%
  filter(Factor == "Factor1")%>%
  #group_by(receiver) %>%
  ungroup()%>%
  slice_max(order_by = abs(interact), prop= 0.1)%>%
  mutate(LR = paste(source_genesymbol, target_genesymbol, sep = "_"))%>%
  pull(LR)
length(LRs)
sort(LRs)
df.comm %>%
  filter(abs(interact)>0.2) %>%
  filter(Factor == "Factor1", 
         sender=="Fib", receiver =="CM", source_genesymbol=="NRG1")
hmap1 = df.comm %>%
  filter(Factor == "Factor1") %>%
  mutate(LR = paste(source_genesymbol, target_genesymbol, sep = "_")) %>%
  filter(LR %in% LRs) %>%
  distinct(sender, receiver, LR, group) %>%
  group_by(sender, LR) %>%
  mutate(n_sender = n()) %>%
  group_by(LR, receiver) %>%
  mutate(n_receiver = n()) %>%
  group_by(LR, group) %>%
  mutate(n_group = n()) %>%
  ungroup()

  # Prepare the data for heatmaps
  hmap_data <- hmap1 %>%
    select(sender, receiver, LR, n_receiver, n_sender) %>%
    distinct() %>% arrange(LR, sender, receiver)
  
  # Convert to wide format for heatmaps
  # Convert to wide format for heatmaps
  matrix_sender <- hmap_data %>%
    distinct(LR, sender, n_sender) %>%
    pivot_wider(names_from = sender, values_from = n_sender, values_fill = 0) %>%
    column_to_rownames("LR") %>%
    as.matrix()
  
  matrix_receiver <- hmap_data %>%
    distinct(LR, receiver, n_receiver) %>%
    pivot_wider(names_from = receiver, values_from = n_receiver, values_fill = 0) %>%
    column_to_rownames("LR") %>%
    as.matrix()
  
  # Prepare stacked barplot data for n_group
  barplot_data <- hmap1 %>%
    distinct(LR, group, n_group) %>% arrange(LR)%>%
    pivot_wider(names_from = group, values_from = n_group, values_fill = list(n_group = 0)) %>%
    column_to_rownames("LR") %>%
    #select(any_of(c("NF", "HF"))) %>%  # Ensure order
    as.matrix()
  
  group_colors <- c("HF" = "black", "NF" = "grey")  # Adjust colors as needed
  group_colors <- group_colors[colnames(barplot_data)]
  
  bar_legend <- Legend(
    labels = names(group_colors),  # "Group1", "Group2"
    legend_gp = gpar(fill = group_colors),  # Apply colors
    title = "Group"
  )
  # Define barplot annotation (stacked)
  bar_anno <- rowAnnotation(
    "n° cell type\npairs" = anno_barplot(
      barplot_data, 
      beside = FALSE,  # Stacked bars
      gp = gpar(fill = group_colors),
      border = TRUE,
      axis_param = list(side = "top",
                        labels_rot = 0)  # Place axis on top for easier reading
    )
  )
  
  
  
  matrix_receiver = matrix_receiver[,sort(colnames(matrix_receiver))]
  matrix_sender = matrix_sender[,sort(colnames(matrix_sender))]
  # Custom color function ensuring 0 is grey
  color_function_sender <- colorRamp2(c(0, 1, max(matrix_sender, na.rm = TRUE)), c("grey", "white", "blue"))
  color_function_receiver <- colorRamp2(c(0, 1, max(matrix_receiver, na.rm = TRUE)), c("grey", "white", "red"))
  
  # Function to display values inside cells, skipping 0s
  cell_text_fun_sender <- function(j, i, x, y, width, height, fill) {
    value <- matrix_sender[i, j]
    if (!is.na(value) && value != 0) {  # Do not print zero values
      grid.text(value, x, y, gp = gpar(fontsize = 8, col = "black"))
    }
  }
  
  cell_text_fun_receiver <- function(j, i, x, y, width, height, fill) {
    value <- matrix_receiver[i, j]
    if (!is.na(value) && value != 0) {  # Do not print zero values
      grid.text(value, x, y, gp = gpar(fontsize = 8, col = "black"))
    }
  }
  
  # Define heatmaps
  heatmap_sender <- Heatmap(matrix_sender, 
                            name = "n_sender",
                            col = color_function_sender,
                            cell_fun = cell_text_fun_sender,
                            #col = colorRampPalette(c("white", "blue"))(50),
                            cluster_rows = FALSE, 
                            cluster_columns = FALSE,
                            show_row_names = TRUE, 
                            show_column_names = TRUE,
                            row_names_side = "left",
                            row_names_gp = gpar(fontsize = 8),
                            show_row_dend = F, 
                            show_column_dend = F,
                             show_heatmap_legend = FALSE )
  heatmap_sender
  heatmap_receiver <- Heatmap(matrix_receiver, 
                              name = "n_receiver",
                              #col = colorRampPalette(c("white", "red"))(50),
                              col = color_function_receiver,
                              cell_fun = cell_text_fun_receiver,
                              cluster_rows = FALSE, 
                              cluster_columns = FALSE,
                              show_row_names = FALSE,  # Avoid duplication
                              show_column_names = TRUE,
                              show_column_dend = F,
                              show_heatmap_legend = FALSE )
  
  # Draw heatmaps next to each other
  draw(heatmap_sender + heatmap_receiver + bar_anno , 
       heatmap_legend_side = "bottom",
       annotation_legend_list = list(bar_legend))
  
pdf("output/figures/comm_lr_network.pdf", height= 10.4, 
    width= 5)
  draw(heatmap_sender + heatmap_receiver + bar_anno , 
     heatmap_legend_side = "bottom",
     annotation_legend_list = list(bar_legend))
dev.off()

## target validation
target.genes.oi <- g.loads%>% 
  filter(Factor=="Factor1",
         ctype == "CM",
         (value) < -0.1)%>% 
  pull(feature)
df.nn%>%
  filter(sender == "Fib" & receiver =="CM",
         test_ligand %in% c("MXRA5", "BMP4", "NRG1"))%>%
  View()
df.nn%>%
  filter(sender == "Fib" & receiver =="CM",
         test_ligand %in% c("MXRA5", "BMP4", "NRG1"))%>%
  filter(target %in% target.genes.oi)%>%
  distinct(test_ligand, target, weight)%>%
  group_by(test_ligand)%>% arrange(test_ligand, -weight)%>%
  write_csv("output/communication/identifying_ligand_targets.csv")
  