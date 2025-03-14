# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de
# modified [2024] Jan D. Lanzer

#' In this script we generate the replicability 
#' main figure from the manuscript.
#' Tile plots are unified from result objects

library(tidyverse)
library(cowplot)

source("make_source_data.R")

METAheart = readRDS(file = "output/reheat1/METAheart2023.rds") #main object
# names(METAheart)<- str_replace_all(names(METAheart), "2019", "19")
# saveRDS(METAheart, "output/reheat1/METAheart2023.rds")
experiments = names(METAheart)
names(experiments) = experiments

new_study_ids= c("Forte22", "Flam19", "Wang22", "Rao21", "Hua19")
old_study_ids <- as.character(experiments[!experiments %in% new_study_ids])
study_order <- c(old_study_ids, new_study_ids)
a <- c(rep("black", length(old_study_ids)), rep("orange4", length(new_study_ids)))

# For labeling
experiment_size = sort(unlist(lapply(METAheart,
                                     function(x) ncol(x$GEX))),
                       decreasing = T)

plot(experiment_size)
sum(experiment_size)

# create a sample size overview:

df= lapply(names(METAheart), function(x){
  METAheart[[x]]$TARGETS%>% select(HeartFailure)%>% mutate(study= x)
})%>% do.call(rbind, .)

df <- df %>% group_by(study, HeartFailure)%>% summarise(x= n())%>% 
  mutate(new= ifelse(study %in% new_study_ids, "y", "n"),
         HeartFailure= factor(HeartFailure, levels= c("yes", "no")))

p.samplesize<- df%>% ggplot(., aes(x= reorder(study, x),y= x, fill= HeartFailure, color= new))+
  geom_col()+
  scale_color_manual(values= rev(c("black", "white")))+
  scale_fill_manual(values = c("#eb4034", "#0b4f32"))+
  coord_flip()+
  labs(y= "", x="")+
  theme_cowplot()

pdf("output/figures/reheat1_samplesize.pdf",
    width= 7, height= 4)
p.samplesize
dev.off()

#1. Jaccard plot
jaccard_df = readRDS(file = "output/reheat1/figure_objects/jaccard_df.rds")

jaccard_tile = ggplot(jaccard_df, aes(x = Var1, 
                                      y = Var2,
                                      fill = value)) +
  geom_tile(color= "darkgrey") +
  theme_minimal() + xlab("Predicted") +
  ylab("Predictor") +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1),
        axis.text = element_text(size=10),
        legend.key.size = unit(.6, "cm"),
        legend.text = element_text(size = 9)) +
  scale_fill_gradientn(colours=RColorBrewer::brewer.pal(9, 'RdPu'), limits=c(0, .35)) +
  labs(fill = "Jaccard \nIndex")

jaccard_tile

#2. AUC of Disease scores

pairwise_200 = readRDS(file = "output/reheat1/figure_objects/pairwise_200.rds")

pcolors = RColorBrewer::brewer.pal(9, 'Purples')[1:6]

pairwise_plot_200 = pairwise_200 %>% ggplot(aes(x=factor(PredictedExperiment,
                                                         levels = rev(study_order)),
                                                y=factor(PredictorExperiment,
                                                         levels = study_order)
                                                ,fill = single)) + 
  geom_tile(color="black") +
  theme_minimal() + xlab("Predicted") +
  ylab("Predictor") +
  theme(
        axis.text = element_text(size=10),
        axis.title = element_text(size=12),
        legend.text = element_text(size = 9),
        axis.text.x = element_text(angle = 90, hjust = 1, colour = a),
        axis.text.y = element_text(colour = rev(a))) +
  #scale_fill_gradientn(colours= pcolors, limits=c(0, 1)) + 
  #scale_fill_gradient(low="white", high= "purple4")+
  scale_fill_gradient2(low = "white", midpoint = 0.5)+
  coord_flip() + labs(fill = "AUROC")

pairwise_plot_200
pdf("output/reheat1/pariwise_AUC_hmap.pdf",
    width = 5.08,
    height = 4)
pairwise_plot_200

dev.off()
#save sourcedata

pairwise_200 %>% 
  rename(AUROC = single) %>% 
  save_source_data(fig_type = T, 
                   fig_number = 2, 
                    panel_letter = "B", 
                   data =  ., description = NULL)
#3. ES plots

pcolors_up = RColorBrewer::brewer.pal(9, 'Reds')[1:7]
pcolors_down = rev(RColorBrewer::brewer.pal(9, 'Blues')[1:7])
pcolors = c(pcolors_down,pcolors_up)

#upregulation
up_ES = readRDS(file = "output/reheat1//figure_objects/up_ES.rds")

up_ES_plot = up_ES %>%
  mutate(sig= ifelse(padj < 0.05, "*", ""), 
         sig = ifelse(Reference==DEG, "" , sig))%>%
  ggplot(aes(x=factor(Reference, levels = rev(study_order)),
             y=factor(DEG, levels = study_order),
             fill = ES,
             label= sig)) + 
  geom_tile(color="white") +
  theme_minimal() + 
  xlab("Reference") +
  ylab("Individual DEG") +
  geom_text(hjust = 0.5, vjust = 0.7, size = 4) +
  theme(#axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text = element_text(size=10),
        axis.title = element_text(size=12),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 9),
        axis.text.x = element_text(angle = 90, hjust = 1, colour = a),
        axis.text.y = element_text(colour = rev(a))) +
  scale_fill_gradientn(colours= pcolors, limits=c(-1, 1)) + 
  coord_flip() + 
  ggtitle("Upregulated genes")

up_ES_plot = up_ES_plot + theme(legend.position = "none")

#downregulation

down_ES = readRDS(file = "output/reheat1//figure_objects/down_ES.rds")

down_ES_plot = down_ES %>% mutate(sig= ifelse(padj < 0.05, "*", ""), 
                                  sig = ifelse(Reference==DEG, "" , sig))%>%
  ggplot(aes(x=factor(Reference, levels = rev(study_order)),
             y=factor(DEG, levels = study_order),
             fill = ES,
             label= sig)) + 
  geom_tile(color="white") +
  theme_minimal() + 
  xlab("Reference") +
  ylab("Individual DEG") +
  geom_text(hjust = 0.5, vjust = 0.7, size = 4) +
  theme(#axis.text.x = element_text(angle = 90, hjust = 1),
    axis.text = element_text(size=10),
    axis.title = element_text(size=12),
    axis.title.x = element_blank(),
    legend.text = element_text(size = 9),
    axis.text.x = element_text(angle = 90, hjust = 1, colour = a),
    axis.text.y = element_text(colour = rev(a))) +
  #scale_fill_gradientn(colours= pcolors, limits=c(-1, 1)) + 
  coord_flip() + 
  coord_flip() + ggtitle("Downregulated genes")


es_legend = get_legend(down_ES_plot)

down_ES_plot = down_ES_plot + theme(legend.position = "none")
down_ES_plot

# Align all plots

left_panel = plot_grid(jaccard_tile, "",
                       align = "v",
                       ncol = 1, rel_heights = c(1,1))
                      #ä labels = c("A","B"))
right_panel = plot_grid(up_ES_plot,down_ES_plot,
                        align = "v",ncol = 1, 
                        rel_heights = c(1,1))
                        #labels = c("C",""))

right_panel = plot_grid(right_panel, es_legend,
                        nrow = 1, rel_widths = c(1,.25))

pdf("output/reheat1/supp_figure_jaccard_ES.pdf",
    width = 10,
    height = 8)

plot(plot_grid(left_panel, right_panel,
          nrow = 1, rel_widths = c(1,1),
          align = "h"))

dev.off()


pdf("output/reheat1/pariwise_AUC_hmap.pdf",
    width = 5,
    height = 4.3)
pairwise_plot_200

dev.off()

# get stats for paper -----------------------------------------------------

add_info<- function(df, col1, col2){
  col1 <- sym(col1)
  col2 <- sym(col2)
  df%>%
  mutate(version= ifelse({{col1}} %in% new_study_ids, "2024", "2021"))%>%
    mutate( {{ col1 }} := factor( {{ col1 }} ,levels= c(new_study_ids, old_study_ids))
           )%>%
    group_by(version)%>%
    mutate(m_score= mean({{col2}}))%>%
    ungroup()
}


pairwise_200<- add_info(pairwise_200, "PredictorExperiment", "single")
down_ES<- add_info(down_ES, "Reference", "ES")
up_ES <- add_info(up_ES, "Reference", "ES")
jaccard_df<- add_info(as_tibble(jaccard_df), "Var1", "value")

jaccard_df%>% distinct(version, m_score)
pairwise_200%>% distinct(version, m_score)
down_ES%>% distinct(version, m_score)
up_ES%>% distinct(version, m_score)

plot_df<- pairwise_200 %>% 
  mutate(new= ifelse(PredictorExperiment %in% new_study_ids, "y", "no"))

# test new vs old 
wilcox.test( plot_df %>% filter(version=="2021" )%>%
          pull(single), 
        plot_df %>% 
          filter(version == "2024")%>% 
          pull(single))

lapply(new_study_ids, function(x){
  wilcox.test( plot_df %>% filter(version=="2021" & !PredictedExperiment %in% new_study_ids )%>%
                 pull(single), 
               plot_df %>% 
                 filter(PredictorExperiment == x & !PredictedExperiment %in% new_study_ids)%>% 
                 pull(single))
  
})

# get mean old vs new
plot_df %>% filter(version=="2021" )%>%
  pull(single)%>% mean() 
plot_df %>% 
  filter(version == "2024")%>% 
  pull(single)%>% mean()
plot_df %>% 
  #filter(version == "2024")%>% 
  pull(single)%>% mean()





p.boxplot_auc <- plot_df%>%
  ggplot(aes(x= PredictorExperiment, y= single, col = version))+
  geom_boxplot()+
  geom_hline(yintercept = 0.5,  color= "darkgrey", linetype=3)+
  coord_flip()+
  theme_cowplot()+
  scale_color_brewer(type = "qual",palette = 2)+
  labs(y="pairwise AUROC", x= "Predictor Study")
p.boxplot_auc

pdf("output/reheat1/pariwise_AUC_box.pdf",
    width = 5,
    height = 5)
  p.boxplot_auc
dev.off()

plot_df%>%
  ggplot(aes(x= new, y= single, col = new))+
  geom_boxplot()+
  geom_hline(yintercept = 0.5,  color= "darkgrey", linetype=3)+
  coord_flip()+
  theme_cowplot()+
  scale_color_brewer(type = "qual",palette = 2)+
  labs(y="AUC", x= "Predictor Study")

