# HEADER --------------------------------------------
#
# Author:     Jan D. Lanzer
# Copyright   Copyright 2024 - Jan D. Lanzer
# Email:      lanzerjan@gmail.com
#
# Date:     2024-08-29
#
# Script Name:    ~/R-projects/reheat2_pilot/meta_data_plots.R
#
# Script Description:
# plot various meta data plots:

library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

col_list<- readRDS("color_list_figures.rds")

meta_data <- readRDS( "output/meta_data_for_plotting.rds")

# etiologies --------------------------------------------------------------

p.etio_summed<- meta_data$etiologies %>% 
  mutate(modality= ifelse(modality=="sc", "Single Nuc","Bulk"))%>%
  filter(heart_failure %in% c("HF", "yes"))%>%
  group_by(disease_code, modality)%>%
  count()%>% 
  group_by(modality)%>%
  mutate(prop= n/sum(n))%>%
  ggplot(aes(x=modality, y= prop , fill =disease_code))+
  geom_col(color="black" )+
  # geom_text(aes(label =scales::percent(prop, accuracy = 1)), 
  #           position = position_stack(vjust = 0.5), # Center the labels within the bars
  #           color = "black")+
  scale_fill_manual(values = col_list$etiologies_colors)+
  theme_cowplot()+
  scale_y_continuous(labels = scales::percent)+
  theme(axis.text.x = element_text(angle= 90, hjust= 1, vjust= 0.5), 
        axis.line = element_blank(), 
        axis.ticks.x = element_blank())+
  labs(x="", y="", fill ="HF\netiology")

p.etio_summed
ggsave(p.etio_summed,
       device = "pdf", 
       filename = "output/figures/fig1_etio_summed.pdf",
       width=4, height=5)

p.etio_per_study<-meta_data$etiologies %>% 
  filter(heart_failure %in% c("HF", "yes"))%>%
  group_by(disease_code, study)%>%
  count()%>% 
  ungroup()%>%
  mutate(prop= n/sum(n))%>% 
  mutate(study= factor(study, levels= rev(sort(unique(meta_data$etiologies$study)))))%>%
  ggplot(aes(y= study,x= disease_code, size = n))+
  geom_point()+
  scale_color_manual(values= c("black", "grey"))+
  #facet_grid(~disease_code)+
  theme_cowplot()+
  theme(axis.text.x= element_text(angle= 90, hjust= 1, vjust=0.5))+
  #guides(color = guide_legend(override.aes = list(size = 4)),)+
  labs(x= "", y="", color="Heart failure")+
  coord_flip()

ggsave(p.etio_per_study,
       device = "pdf", 
       filename = "output/figures/fig1_etio_study_pdf",
       width=5, height=5)


# sex ---------------------------------------------------------------------
fem_hf<-mean(subset(x= meta_data$demo$proportion_female,
            subset = meta_data$demo$heart_failure=="HF"), na.rm = T)

fem_nf<- mean(subset(x= meta_data$demo$proportion_female,
            subset = meta_data$demo$heart_failure=="NF"), na.rm = T)

p.sex_summed<- meta_data$demo %>% 
  mutate(modality= ifelse(modality=="sc", "Single Nuc","Bulk"))%>%
  mutate(group= paste0(modality, "_", heart_failure))%>%
    ggplot(aes(x=group, y= proportion_female/100 ))+
  geom_hline(yintercept= 0.5)+
  geom_boxplot(aes( fill =heart_failure), show.legend = F)+
  scale_y_continuous(labels = scales::percent)+
  geom_jitter(width=0.1)+
  theme(axis.text.x = element_text(angle= 90, hjust= 1, vjust= 0.5))+
  #scale_fill_manual(values = unname(col_list$etiologies_colors))+
  scale_fill_manual(values = c("#387780", "darkgrey"))+
  labs(x="", y="Female patients", fill ="HF\netiology")
p.sex_summed

p.sex_per_study<- meta_data$demo %>% 
  arrange(modality)%>%
  mutate(study= factor(study, levels= rev(sort(unique(meta_data$demo$study)))))%>%
  #mutate(study= factor(study, levels= c( y,x)))%>%
  ggplot(aes(y=study,  x= proportion_female, fill = heart_failure))+
  geom_col(position = "dodge", width = 0.7, color="black")+
  coord_flip()+
  geom_hline(yintercept = 50)+
  theme_cowplot()+
  scale_fill_manual(values= c("black", "grey"))+
  labs(y="", x= "female sex (%)")+
  theme(axis.text.x= element_text(angle= 90, hjust= 1, vjust=0.5))

p.sex_per_study

# age-------------------------------------------------------------------------

p.age_box<- 
  meta_data$demo %>% 
  mutate(modality= ifelse(modality=="sc", "Single Nuc","Bulk"))%>%
  mutate(group= paste0(modality, "_", heart_failure))%>%
  ggplot(aes(x=group, y= age_mean))+
  geom_boxplot(aes( fill =heart_failure), show.legend = F, outlier.color = NA)+
  geom_jitter(width=0.1)+
  theme(axis.text.x = element_text(angle= 90, hjust= 1, vjust= 0.5))+
  scale_fill_manual(values = c("#387780", "darkgrey"))+
  labs(x="", y="Age (y)\nstudy mean", fill ="HF\netiology")

p.age_box

p.age<- meta_data$demo%>%
  mutate(label= ifelse(heart_failure=="NF","",study ))%>%
  ggplot(aes(x=year, y=proportion_female, color=heart_failure))+
  geom_hline(yintercept = mean(dfcomb%>% filter(heart_failure=="HF")%>% pull(proportion_female), na.rm = T),
             linetype=3, 
             linewidth=0.8,
             color="black")+
  geom_hline(yintercept = mean(dfcomb%>% filter(heart_failure!="HF")%>% pull(proportion_female), na.rm = T),
             linetype=3,
             linewidth=0.8,
             color="#b83b54")+
  geom_col(position="dodge", fill = "grey", alpha= 0.5)+
  ggrepel::geom_label_repel(aes(label=label, fill= modality),color="white",
                            max.overlaps = 2000,
                            show.legend = T, 
                            alpha= 0.8
  )+
  geom_smooth(formula = y~x, na.rm = T, alpha= 0.2, se = F)+
  scale_color_manual(values= c("#E83151", "#757780"))+
  scale_fill_manual(values= c("darkblue", "black"))+
  theme_cowplot()+
  ylim(c(0,100))+
  xlim(c(2005,2023))+
  labs(x="Year", y= "Female Patients\nProportion per study (%)", 
       fill="Modality",color= "Cohort")+
  theme(panel.background = element_rect(fill="lavender"))
p.age

p.age<- dfcomb%>%
  mutate(label= ifelse(heart_failure=="NF","",study ))%>%
  ggplot(aes(x=year, y=age_mean,  color=heart_failure))+
  geom_hline(yintercept = mean(dfcomb%>% filter(heart_failure=="HF")%>% pull(age_mean), na.rm = T),
             linetype=3,linewidth=0.8, color="black")+
  geom_hline(yintercept = mean(dfcomb%>% filter(heart_failure!="HF")%>% pull(age_mean), na.rm = T), 
             linetype=3, linewidth=0.8, color="darkgreen")+
  geom_errorbar(aes(ymin=age_mean-age_sd, ymax=age_mean+age_sd),alpha= .4, 
               linetype=1,
              linewidth=.4)+
  ggrepel::geom_label_repel(aes(label=label, fill= modality),color="white",
                            max.overlaps = 2000,
                            show.legend = T, 
                            alpha= 0.8  )+
  geom_smooth(formula = y~x, na.rm = T, alpha= 0.2, se= F)+
  scale_color_manual(values= c("#E83151", "#757780"))+
  scale_fill_manual(values= c("darkblue", "black"))+
  geom_point()+
  theme_cowplot()+
  ylim(c(25,80))+
  xlim(c(2005,2023))+
  labs(y="Age\n(study mean ± sd)", x="Year",  fill="Modality",color= "Cohort")+
  theme(panel.background = element_rect(fill="lavender"))
p.age

p.age<- meta_data$demo%>%
  mutate(label= ifelse(heart_failure=="NF","",study ))%>%
  ggplot(aes(x=heart_failure, y=age_mean,  fill=heart_failure))+
  geom_hline(yintercept = mean(meta_data$demo%>% filter(heart_failure=="HF")%>% pull(age_mean), na.rm = T),
             linetype=3,linewidth=0.8, color="black")+
  geom_hline(yintercept = mean(meta_data$demo%>% filter(heart_failure!="HF")%>% pull(age_mean), na.rm = T), 
             linetype=3, linewidth=0.8, color="grey")+
  geom_errorbar(aes(ymin=age_mean-age_sd, ymax=age_mean+age_sd),alpha= .4, 
                linetype=1,
                linewidth=.8)+
  facet_wrap(~study, ncol = 13)+
  # ggrepel::geom_label_repel(aes(label=label, fill= modality),color="white",
  #                           max.overlaps = 2000,
  #                           show.legend = T, 
  #                           alpha= 0.8  )+
  # #geom_smooth(formula = y~x, na.rm = T, alpha= 0.2, se= F)+
  scale_color_manual(values= c("#E83151", "#757780"))+
  scale_fill_manual(values= c("grey", "black"))+
  geom_col(position= "dodge")+
  theme_cowplot()+
  #ylim(c(25,80))+
  #xlim(c(2005,2023))+
  labs(y="Age\n(study mean ± sd)", x="",  fill="",color= "Cohort")+
  theme(panel.background = element_rect(fill="white"),
        axis.text.x= element_blank(),
        axis.ticks.x= element_blank(), 
        strip.text = element_text(size= 7))
  #coord_flip()
p.age
# NA matrix ---------------------------------------------------------------




# sample_location and accquisition---------------------------------------------------------

df_location_sum<-meta_data$location%>%
  group_by(modality, sample_site, heart_failure)%>% summarise(n= sum(n))

p_biopsy_location_summed<- meta_data$location %>% 
  mutate(modality= ifelse(modality=="sc", "Single Nuc","Bulk"))%>%
  mutate(sample_site= str_replace_all(sample_site, "LV_", ""))%>%
  filter(heart_failure %in% c("HF", "yes"))%>%
  group_by(sample_site, modality)%>%
  count()%>% 
  group_by(modality)%>%
  mutate(prop= n/sum(n))%>%
  ggplot(aes(x= modality,y= prop, fill = sample_site))+
  geom_col(position = "fill", color= "black")+
  # geom_text(aes(label =scales::percent(prop, accuracy = 0.1)), 
  #           position = position_stack(vjust = 0.5), # Center the labels within the bars
  #           color = "black")+
  #scale_fill_manual(values= c("#428af5","#e08816", "grey"))+
  scale_fill_manual(values=unname(col_list$etiologies_colors[-1]))+
  theme_cowplot()+
  scale_y_continuous(labels = scales::percent)+
  theme(axis.text.x = element_text(angle= 90, hjust= 1, vjust= 0.5), 
        axis.line = element_blank(), 
        axis.ticks.x = element_blank())+
  #theme(axis.text.x = element_text(angle= 90, hjust= 1))+
  labs(x= "", y= "", fill = "LV biopsy\nlocation")
p_biopsy_location_summed


# get proportions of sample loc 
# and sample acc
rbind(bulk_sample_acc_df, sc_sample_acc_df)%>%
  filter(heart_failure=="HF")%>%
  group_by(modality, sample_accquisition)%>%
  summarise(n2= sum(n, na.rm = T))%>%
  mutate(prop= n2/sum(n2))
filter(!is.na(n2))
rbind(bulk_sample_acc_df, sc_sample_acc_df)%>%
  filter(heart_failure=="HF")%>%
  group_by( sample_accquisition)%>%
  summarise(n2= sum(n, na.rm = T))%>%
  mutate(prop= n2/sum(n2))
filter(!is.na(n2))


df_location %>%
  filter(heart_failure=="HF")%>%
  group_by(modality, sample_site)%>%
  summarise(n2= sum(n, na.rm = T))%>%
  mutate(prop= n2/sum(n2))

df_location %>%
  filter(heart_failure=="HF")%>%
  group_by(sample_site)%>%
  summarise(n2= sum(n, na.rm = T))%>%
  mutate(prop= n2/sum(n2))

df_location %>%
  group_by(modality, sample_site)%>%
  summarise(n2= sum(n, na.rm = T))%>% 
  filter(!is.na(n2))

# Create test data.

data <- df_location %>%
  filter(modality=="bulk")%>%
  group_by(modality, sample_site)%>%
  summarise(n2= sum(n, na.rm = T))%>% 
  filter(!is.na(n2))

# Compute percentages
data$fraction = data$n2 / sum(data$n2)
# Compute the cumulative percentages (top of each rectangle)
data$ymax <- cumsum(data$fraction)

# Compute the bottom of each rectangle
data$ymin <- c(0, head(data$ymax, n=-1))

# Compute label position
data$labelPosition <- (data$ymax + data$ymin) / 2

# Compute a good label
data$label <- paste0(data$sample_site, "\n value: ", round(data$fraction, 2))
bulk.p <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=sample_site)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=2) +
  scale_fill_brewer(palette=3) +
  coord_polar(theta="y") +
  #xlim(c(2, 4)) +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")


data <- df_location %>%
  filter(modality=="sc")%>%
  group_by(modality, sample_site)%>%
  summarise(n2= sum(n, na.rm = T))%>% 
  filter(!is.na(n2))

# Compute percentages
data$fraction = data$n2 / sum(data$n2)
# Compute the cumulative percentages (top of each rectangle)
data$ymax <- cumsum(data$fraction)

# Compute the bottom of each rectangle
data$ymin <- c(0, head(data$ymax, n=-1))

# Compute label position
data$labelPosition <- (data$ymax + data$ymin) / 2

# Compute a good label
data$label <- paste0(data$sample_site, "\n value: ", round(data$fraction, 2))

sc.p <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=sample_site)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=2) +
  scale_fill_brewer(palette=3) +
  coord_polar(theta="y") +
  #xlim(c(2, 4)) +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")

data <- df_location %>%
  #filter(modality=="sc")%>%
  group_by( sample_site)%>%
  summarise(n2= sum(n, na.rm = T))%>% 
  filter(!is.na(n2))%>%
  arrange(desc(sample_site))
data
# Compute percentages
data$fraction = data$n2 / sum(data$n2)
# Compute the cumulative percentages (top of each rectangle)
data$ymax <- cumsum(data$fraction)

# Compute the bottom of each rectangle
data$ymin <- c(0, head(data$ymax, n=-1))

# Compute label position
data$labelPosition <- (data$ymax + data$ymin) / 2

# Compute a good label
data$label <- paste0(data$sample_site, "\n value: ", round(data$fraction, 2))

sc.p <-data%>% arrange(desc(label))%>% 
  mutate(sample_site = ifelse(is.na(sample_site), "unkown", sample_site))%>%
  
  ggplot( aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=sample_site)) +
  geom_rect(color="black") +
  geom_label( x=3.5, aes(y=labelPosition, label=label),alpha=1,  size=3) +
  scale_fill_brewer(palette=3) +
  coord_polar(theta="y") +
  #xlim(c(2, 4)) +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")
sc.p


# accquisition ---------------------------------------------------------------------

p.sample_acc_summed<- meta_data$accquisition%>%
  mutate(modality= ifelse(modality=="sc", "Single Nuc","Bulk"))%>%
  filter(heart_failure %in% c("HF", "yes"))%>%
  group_by(sample_accquisition, modality)%>%
  count()%>% 
  group_by(modality)%>%
  mutate(prop= n/sum(n))%>%
  ggplot(aes(x= modality,y= prop, fill = sample_accquisition))+
  geom_col(position = "fill", color= "black")+
  # geom_label(aes(label =scales::percent(prop, accuracy = 1)), 
  #            label.padding = unit(0.25, "lines"),
  #            fill="white", 
  #            position = position_stack(vjust = 0.5), # Center the labels within the bars
  #            color = "black")+
  #scale_fill_manual(values= c("#428af5","#e08816", "grey"))+
  scale_fill_manual(values=c(  "#B16D91",   "#D2CCA1"))+
                      #unname(col_list$etiologies_colors))+
  theme_cowplot()+
  scale_y_continuous(labels = scales::percent)+
  theme(axis.line = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.x= element_text(angle= 90, hjust= 1, vjust= 0.5))+
  labs(x= "", y= "", fill = "Reason for\nHF biopsy")
p.sample_acc_summed

p.sample_acc_study<- 
  meta_data$accquisition%>%
  mutate(study= factor(study, levels= rev(sort(unique(meta_data$etiologiess$study)))))%>%
  ggplot(aes(y= study,x= heart_failure, size = n, color = heart_failure))+
  geom_point()+
  scale_color_manual(values= c("black", "grey"))+
  facet_grid(~sample_accquisition)+
  theme_cowplot()+
  coord_flip()+
  theme(axis.text.x = element_text(angle= 90, hjust= 1,vjust= 0.5))+
  guides(color = guide_legend(override.aes = list(size = 4)),
  )+
  labs(x= "", y="", color="Heart failure")
p.sample_acc_study

p.sample_loc<- plot_grid(p_biopsy_location,
                         p.sample_acc, align= "v", ncol = 1)

pdf("output/figures/meta_data_presentation/biopsy_location.pdf",
    width= 3, height= 5)
p.sample_loc
dev.off()




# genetic -----------------------------------------------------------------

p.genetic_summed<- meta_data$target_sc%>%
  #distinct(genetic_HF, study)%>%
  filter(heart_failure=="HF")%>%
  group_by(genetic_HF, study)%>%count()%>%
  group_by(study)%>%
  mutate(prop=  n / sum(n))%>%
  ggplot(aes(x= study, y = prop, fill = genetic_HF))+
  geom_col(color="black")+
  scale_y_continuous(labels = scales::percent)+
  theme(axis.line = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.x= element_text(angle= 90, hjust= 1, vjust= 0.5))+
  labs(x="", y="", fill ="familial/\ngenetic\nHF")+
  # geom_text(aes(label =scales::percent(prop, accuracy = 1)), 
  #            angle= 90,
  #            position = position_stack(vjust = .5), # Center the labels within the bars
  #            color = "black")+
  # geom_label(aes(label =scales::percent(prop, accuracy = 1)), 
  #           label.padding = unit(0.25, "lines"),
  #           fill="white", 
  #           position = position_stack(vjust = .1), # Center the labels within the bars
  #           color = "black")+
  scale_fill_manual(values =unname(col_list$etiologies_colors)[2:3])
p.genetic_summed

p.genetic_per_study<- meta_data$target_sc%>%
  filter(genetic_HF == "yes")%>%
  #distinct(genetic_HF, study)%>%
  group_by(pv, study)%>%count()%>%
  group_by(study)%>%
  mutate(prop=  n / sum(n))%>%
  ggplot(aes(x= study, y = prop, fill = pv))+
  geom_col(color="white")+
  geom_text(aes(label =pv), 
            position = position_stack(vjust = 0.5), # Center the labels within the bars
            color = "black")+
  scale_y_continuous(labels = scales::percent)+
  theme(axis.line = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.x= element_text(angle= 90, hjust= 1), 
        legend.position = "none")+
  labs(x="", y="", fill ="pathogenic variant")
p.genetic_per_study


# lvef, bmi and race ------------------------------------------------------
p.lvef<- meta_data$target_sc %>%
  ggplot(aes(x = study, y = LVEF, fill = heart_failure)) +
  geom_quasirandom(
    #aes(color = heart_failure),  # Color points by heart_failure to match the fill of the violins
    size = 1, 
    width = 0.1, 
    varwidth = FALSE,  # Ensure uniform width for the beeswarm plot
    dodge.width = 0.8,  # Align points with violin plot using dodge.width
    position = position_dodge(0.8)  # Position the points to match the dodge position
  ) +
  geom_boxplot(position = position_dodge(0.8), alpha = 0.6) +  # Make violins slightly transparent
  scale_fill_manual(values = c("#387780", "darkgrey"))+  # Custom fill colors for the violins
  scale_color_manual(values = c("grey", "black")) +  # Matching custom colors for the points
  #theme_minimal() +
  labs(x = "", y = "LVEF (%)", fill = "Heart Failure", color = "Heart Failure")+
  theme(#axis.line = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.x= element_text(angle= 90, hjust= 1, vjust= 0.5))
p.lvef
# Load necessary libraries
library(ggplot2)
library(ggbeeswarm)
p.bmi<- meta_data$target_sc %>%
    ggplot(aes(x = study, y = BMI, fill = heart_failure)) +
    geom_quasirandom(
      #aes(color = heart_failure),  # Color points by heart_failure to match the fill of the violins
      size = 1, 
      width = 0.1, 
      varwidth = FALSE,  # Ensure uniform width for the beeswarm plot
      dodge.width = 0.8,  # Align points with violin plot using dodge.width
      position = position_dodge(0.8)  # Position the points to match the dodge position
    ) +
    geom_boxplot(position = position_dodge(0.8), alpha = 0.6) +  # Make violins slightly transparent
  scale_fill_manual(values = c("#387780", "darkgrey"))+ # Custom fill colors for the violins
    scale_color_manual(values = c("grey", "black")) +  # Matching custom colors for the points
    labs(x = "", y = "BMI (kg/m²)", fill = "Heart Failure", color = "Heart Failure")+
  theme(#axis.line = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.text.x= element_text(angle= 90, hjust= 1, vjust= 0.5))
p.bmi  
# sample size  ------------------------------------------------------------

p.sample_size<- meta_data$sample_size %>%
  mutate(new = ifelse(new == "y", "2024", "2021"),
         modality = ifelse(modality=="bulk", "Bulk", "Sn"),
         label = ifelse(new == "2024" & heart_failure == "HF", "*", "")) %>%
  ggplot(aes(x = reorder(study, n), y = (n), fill = heart_failure)) +
  geom_hline(yintercept = c(0, 10,20, 30,40,  50, 100, 200, 300), color = "darkgrey") +
  geom_col(width = 0.8, aes(color = new), show.legend = F) +
  scale_color_manual(values = rev(c("black", "#F9F7F3"))) +
  scale_fill_manual(values = c("#387780", "darkgrey"))+
  #scale_fill_manual(values = rev(c("grey", "black"))) + 
  coord_flip() +
  #geom_text(aes(label = label), hjust = -0.2, vjust = 1, size = 10, color="red") +  # Adds stars at the top
  facet_grid(rows = vars(modality), scales = "free", space = 'free') +
  labs(y = "", x = "") +
  theme_cowplot() +
  scale_y_continuous(breaks = c(10, 30,  50, 100, 200, 300)) +
  theme(legend.position = "right",
        strip.background =element_rect(fill="lavender"),
        axis.text.x = element_text(angle = 90, hjust = 1,
                                   size = 12, vjust= 0.5)) +
  labs(fill = "Cohort", color = "Version")

 
p.sample_size
# combine -----------------------------------------------------------------

## main figure
legend <- get_legend(p.lvef)

p.compo2<-plot_grid(
  plot_grid(p.sample_size, 
          plot_grid(p.age_box,
                    p.sex_summed,
                    ncol = 1),
          rel_widths = c(1, 0.25)
          ),
  plot_grid(ncol = 3,
            p.etio_summed,
            p_biopsy_location_summed,
            p.sample_acc_summed, axis = "tb", align ="hv"
            ),
  plot_grid(p.genetic_summed, 
            plot_grid(ncol=3,
                      p.lvef+theme(legend.position = "none"),
                      p.bmi+theme(legend.position = "none"),
                      legend,
                      rel_widths = c(1,1,0.4)),
            nrow= 1, 
            rel_widths = c(1,1.6)),
  ncol = 1, 
  rel_heights = c(1,0.6, 0.6)
  )
#p.compo2<- plot_grid(p.sample_size, p.compo, ncol= 1, rel_heights = c(1,1.7))
p.compo2
pdf("output/figures/meta_data_presentation/meta_composite.pdf", 
    width= 9.4, height= 11.2)
p.compo2
dev.off()


  ## supp figure
p.sample_acc_study
p.genetic_per_study
p.age
p.sex_per_study
p.etiogen<- plot_grid(p.etio_per_study, p.genetic_per_study,ncol = 2, rel_widths = c(1,0.35))
p.supp<- plot_grid(ncol =1,
          align ="hv",
  p.age,
  p.sex_per_study,
  p.etiogen,
  p.sample_acc_study,
  rel_heights = c(1,0.7, 0.8)
)
p.etiogen
pdf("output/figures/meta_data_presentation/meta_composite_supp.pdf", 
    width= 10, height= 15)
p.supp
dev.off()


# get numbers for paper ---------------------------------------------------
meta_data$sample_size%>%
  group_by(modality)%>% 
  summarise(n2= sum(n))%>% mutate(total = sum(n2))
meta_data$sample_size %>% pull(study)%>% unique()%>% length()
meta_data$sample_size %>% View()
meta_data$target_sc%>% 
  group_by(genetic_HF)%>%
  filter(heart_failure =="HF", disease_code != "ICM")%>%
count()%>%ungroup() %>%  mutate(prop= n/sum(n))

meta_data$target_sc %>% 
  filter(heart_failure =="NF")%>%
  summarise(ef= median(LVEF, na.rm= T),
            bmi= median(BMI, na.rm=T))

meta_data$target_sc %>% 
  #filter(heart_failure =="NF")%>%
  count(race)%>% 
  ungroup()%>% 
  mutate(prop= n/sum(n))
