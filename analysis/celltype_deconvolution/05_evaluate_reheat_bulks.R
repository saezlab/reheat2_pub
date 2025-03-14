## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2023-01-24
##
## Copyright (c) Jan D. Lanzer, 2023
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
## Evaluate the resulsts of the deconvolution
## For different mixtures, different strategies are used for evaluation:
## 1. for artifical mixtures from sc data, we will calculate pearson and RMSE on this 
## ground truths
## 2. for real bulks, we will evaluate mean composition changes per study and compare separability
## of HF and CT bulks

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(cowplot)
library(ggrepel)
library(cluster)
library(lmerTest)
source("make_source_data.R")
source("celltype_deconvolution/utils_decon.R")
source("aesthetics.R")

color_list<-readRDS("color_list_figures.rds")

#deconvo results:
res <- readRDS("output/deconvo/cibersort_results_reheat23.rds")
names(res) =str_replace_all(names(res), "2019", "19")
meta_heart<-readRDS("output/reheat1/Metaheart_counts2023.rds")

meta_heart2 <-readRDS("output/reheat1/METAheart2023.rds")
meta_heart$Hua2019$target<-meta_heart2$Hua19$TARGETS
meta_heart$Hua19 <- meta_heart$Hua2019
meta_heart$Hua2019<- c()
sig.matrices_names<- c("comp")
## first, we take the mean if there are multiple runs (subsample runs per matrix)
res= lapply(res, function(study){
 res=  map(sig.matrices_names, function(sig){
    #names(study)[grepl(sig, names(study))]
    require(purrr)
    x= study[grepl(sig, names(study))]
    v = reduce(x, `+`) / length(x)
  })
 names(res)= sig.matrices_names
 return(res)
})

# creating a large tidy df -------------------------------------

target.df <- lapply(names(meta_heart), function(study){
  meta_heart[[study]]$target%>% dplyr::select(Sample, HeartFailure)%>%mutate(study= study)
})%>%  do.call(rbind, .)
target.df = target.df %>% mutate(study= str_replace_all(study, "2019", "19"))

M <-lapply(sig.matrices_names, function(sig){
  x <- unlist(res, recursive = F)
  x <- x[grepl(x = names(x), pattern = sig)]
  x <- do.call(rbind, x)
})
names(M)<- sig.matrices_names

# we only use the comp matrix 
M<- M$comp[, wanted_cells]%>% 
  as.matrix()
#translate 
colnames(M)<- cell_dic[match(colnames(M), cell_dic$pred),"real"]

write_csv(as.data.frame(M)%>% rownames_to_column("sample_id"), "output/deconvo/reheat_devonvoluted_composition.csv")
write_csv(as.data.frame(target.df)%>% rename(sample_id = "Sample"), "output/deconvo/reheat_devonvoluted_meta.csv")


# calculate whether compositions associate with hf or study label ---------

# Here we can calculate the sillouhette score based on distances in the compositions

patient_dists <- stats::dist(M[target.df$Sample,],
                      method = "euclidean")
target.df <- target.df%>% mutate(HF_dummy = as.factor(HeartFailure) %>%
                                   as.integer(),
                                 study_dummy = as.factor(study) %>%
                                   as.integer())

si_hf <- (cluster::silhouette(x = target.df$HF_dummy, dist=  patient_dists))%>%
  as.data.frame() %>%
  dplyr::select(cluster, sil_width) %>%
  cbind(target.df[, c("Sample", "HeartFailure")])%>%
  dplyr::mutate(covar = "Disease") %>%
  dplyr::rename("var_name" = HeartFailure)%>%
  as_tibble()
si_hf
si_study <- (cluster::silhouette(x = target.df$study_dummy, dist=  patient_dists))%>%
  as.data.frame() %>%
  dplyr::select(cluster, sil_width) %>%
  cbind(target.df[, c("Sample", "study")])%>%
  dplyr::mutate(covar = "Study") %>%
  dplyr::rename("var_name" = study)%>%
  as_tibble()

si <- bind_rows(si_hf, si_study)%>% as_tibble()
si
# Plot

si_plt <- si %>%
  ggplot(aes(x = var_name, y = sil_width)) +
  geom_boxplot() +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10)) +
  ylab("silhouette score") +
  facet_grid(.~covar,scales = "free", space='free') +
  xlab("")
si_plt
pdf("output/figures/deconv_reheat_si.pdf", height = 3, width = 5)

plot(si_plt)

dev.off()

# Finding which groups are influenced?

si %>%
  group_by(var_name) %>%
  nest() %>%
  mutate(zero_dif = map(data, function(dat){
    
    t.test(dat$sil_width, mu = 0, alternative = "greater") %>%
      broom::tidy()
    
  })) %>%
  dplyr::select(-data) %>%
  unnest() %>%
  ungroup() %>%
  dplyr::mutate(adj_pval = p.adjust(p.value))

# Numbers for the paper

si %>%
  group_by(var_name) %>%
  summarize(median(sil_width))



# run t-test per study of composition changes in HF -----------------------

## for each study we want to calculate a statistic for separating hf from nf
# study="Hua19"
# sig= "comp"

t.test.res <-lapply(names(res), function(study){
  print(study)
  table(rownames(res[[study]]$comp) %in% meta_heart[[study]]$target$Sample)
  
  lapply(sig.matrices_names, function(sig){
    prop_data_clr <-res[[study]][[sig]][,wanted_cells]%>%
      as.matrix()%>%
      compositions::acomp() %>%
      compositions::clr()
    
    study_diff_stats <- prop_data_clr %>%
      as.data.frame() %>%
      rownames_to_column("Sample") %>%
      pivot_longer(-Sample, names_to = "cell_type") %>%
      left_join(meta_heart[[study]]$target, by = "Sample") %>%
      mutate(HeartFailure= factor(HeartFailure, levels =c("yes", "no")))%>%
      group_by( cell_type) %>%
      nest() %>%
      dplyr::mutate(HF_diff = map(data, function(dat) {
        
        t.test(value ~ HeartFailure, data = dat) %>%
          broom::tidy()
        
      })) %>%
      dplyr::select(-data) %>%
      unnest(HF_diff) %>%
      mutate(signature_matrix_category= sig, 
             mix_name= study)
    
      
  })%>% do.call(rbind,.)
  

  
})%>% do.call(rbind,.)

# plot those test statistics
stats= t.test.res %>% 
  dplyr::select(cell_type, statistic, signature_matrix_category, mix_name, p.value)
stats$cell_type<- cell_dic[match(stats$cell_type, cell_dic$pred),"real"]
#boxplot
stats %>%
  ggplot(aes(x= cell_type, y= statistic, fill = cell_type))+
  geom_hline(yintercept = 0)+
  geom_boxplot()+
  geom_point(aes(color= mix_name, alpha= 0.5))+
  facet_grid(~signature_matrix_category)+
  theme_cowplot()+

  theme(axis.text.x= element_blank())

#same but as hmap:
hmap_study_wise <- stats %>%
  mutate(star= ifelse(p.value<0.05, "*", ""))%>%
  ggplot(aes(x= cell_type, y= mix_name, fill = statistic))+
  geom_tile(color="black", size= 0.5)+
  geom_text(aes(label=star))+
  theme_cowplot()+
  scale_fill_manual(values= unname(color_list$general_up_down_colors))+
  coord_flip()+
  scale_fill_gradient2(low = color_list$general_up_down_colors[1],  mid= "white", high =color_list$general_up_down_colors[2])+
  theme(axis.text.x= element_text(angle= 45, hjust= 1),
        axis.line = element_blank())+
  labs(x="", y="", fill = "t-value\nHF - NF")
pdf("output/figures/deconv_reheat_tstats.pdf", 
    width= 4.5, height= 3)
hmap_study_wise
dev.off()

map(unique(stats$signature_matrix_category), function(sig){
  stats2<-stats %>%
    filter(signature_matrix_category==sig)%>%
    pivot_wider(id_cols= -signature_matrix_category, 
                names_from = cell_type, values_from= estimate)%>%
    as.data.frame() %>%
    column_to_rownameslibrary(ComplexHeatmap)("mix_name")%>%
    as.matrix()%>%
    Heatmap( name = sig)
})



# perform a linear mixed model to get a single change statistic --------------------------------------------

concat_M <-lapply(names(res), function(study){
  print(study)
  table(rownames(res[[study]]$allGenes) %in% meta_heart[[study]]$target$Sample)
  prop_data_clr <-res[[study]][[sig]][,wanted_cells]%>%
    as.matrix()%>%
    compositions::acomp() %>%
    compositions::clr()%>%
    as.data.frame() %>%
    rownames_to_column("Sample")%>%
    mutate(study =study)%>%
    left_join(target.df, by= c("Sample", "study"))
})%>%
  do.call(rbind, .)

study_diff_stats_lmer <- map(wanted_cells, function(cell){
  model_df <- concat_M[,c(cell, "study", "HeartFailure")]
  fit= lmerTest::lmer(formula = as.formula(paste0(cell, " ~ HeartFailure + ( 1 | study)")),
              data= model_df)
  model_res <- fit%>% summary()
  
  fixed_summ <- model_res$coefficients %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    dplyr::filter(grepl(pattern = "HeartFailureyes", term)) %>%
    dplyr::rename("p_val" = "Pr(>|t|)")
  
  random_summ <- model_res$varcor %>%
    as.data.frame() %>%
    dplyr::select(grp, vcov) %>%
    pivot_wider(names_from = grp, values_from = vcov) %>%
    dplyr::mutate(perc_studyvar = study/(study + Residual)) %>%
    pull(perc_studyvar)
  
  all_res <- fixed_summ %>%
    dplyr::mutate(perc_studyvar = random_summ)%>%
    dplyr::mutate(cell_type = cell)
  
  return(all_res)
}) %>%
  do.call(rbind,.)%>%
  ungroup() %>%
  dplyr::mutate(adj_pval = p.adjust(p_val))

study_diff_stats_lmer$cell_type
#translate 
study_diff_stats_lmer$cell_type<- cell_dic[match(study_diff_stats_lmer$cell_type, cell_dic$pred),"real"]

write_csv(df, "output/lmer_bulk_devonvo_res.csv")

ct_order <- study_diff_stats_lmer %>%
  arrange(desc(abs(Estimate))) %>%
  pull(cell_type)

lmer_plot <- study_diff_stats_lmer %>%
  dplyr::mutate(direction = ifelse(Estimate<0, "decrease", "increase"),
                cell_type = factor(cell_type, levels = ct_order)) %>%
  ggplot(., aes(x = cell_type, y = abs(Estimate), fill = direction,label = round(adj_pval, 3))) +
  geom_bar(stat = "identity") +
  geom_text(vjust= 0)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_manual(values = unname(color_list$general_up_down_colors)) +
  xlab("")
lmer_plot

pdf("output/figures/lmer_test_bulk_deconv.pdf", 
    width= 4, height= 3)
lmer_plot
dev.off()


# plot heatmap ------------------------------------------------------------

stats2<-stats %>%
  select(-signature_matrix_category, -p.value)%>%
  #filter(signature_matrix_category==sig)%>%
  mutate(statistic= ifelse(is.na(statistic), 0, statistic))%>%
  pivot_wider(names_from = cell_type, values_from= statistic)%>%
  as.data.frame() %>%
  column_to_rownames("mix_name")%>%
  as.matrix()

p.vals<-stats %>%
  select(-signature_matrix_category, -statistic)%>%
  #filter(signature_matrix_category==sig)%>%
  mutate(p.value= ifelse(is.na(p.value), 1, p.value))%>%
  pivot_wider(names_from = cell_type, values_from= p.value)%>%
  as.data.frame() %>%
  column_to_rownames("mix_name")%>%
  as.matrix()


study_diff_stats_lmer

save_source_data(T, 5, "G2", data= study_diff_stats_lmer, 
                 description = "results from the linear mixed model",
                 row_names =F)

color_fun <- function(estimates, p_vals) {
  ifelse(
    p_vals < 0.05,  # Check if p-value is less than 0.05
    ifelse(estimates > 0, 
           color_list$general_up_down_colors["up"], 
           color_list$general_up_down_colors["down"]
    ),  # Positive estimates are blue, negative are red
    "grey"  # Use grey if p-value is not significant
  )
}

heatmap_color_fun <- colorRamp2(
  c(min(stats2),-5, 0,5, max(stats2)),  # Define the range of the data
  #c(color_list$general_up_down_colors["down"],"blue",  "white","red",  color_list$general_up_down_colors["up"])
  c("darkblue","blue",  "white","red", "darkred")
)

# Create the bar plot annotation with color coding
column_ha <- HeatmapAnnotation(
  estimate = anno_barplot(
    study_diff_stats_lmer$Estimate,
    gp = gpar(fill = color_fun(study_diff_stats_lmer$Estimate, study_diff_stats_lmer$p_val)),
    ylim = c(min(study_diff_stats_lmer$Estimate, 0),
             max(study_diff_stats_lmer$Estimate, 0)),
    height = unit(1.5, "cm") 
    
  )
)

star_matrix <- (ifelse(p.vals < 0.05, "*", ""))

#stats2[(stats2) > 5] = 5
#stats2[(stats2) < -5] = -5

# Create the heatmap with the left annotation and stars
reheat1_comp_changes<- ComplexHeatmap::Heatmap(
  (stats2),
  name = "t-value\nHF-NF",
  rect_gp = gpar(col = "black", lwd = 1),
  top_annotation = column_ha,
  show_column_dend = F, 
  show_row_dend = F, 
  col = heatmap_color_fun,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(star_matrix[i, j], x, y, gp = gpar(fontsize = 12, col = "black"))
  }
)
reheat1_comp_changes  

#save source data
stats2%>%
  as.data.frame()%>%
  rownames_to_column("study")%>%
  save_source_data(T, 5, "G1", data= .,
                 description = "Cell type compositon changes", row_names = F)

pdf("output/figures/deconv_reheat_hmap_clusterd.pdf",
    width = 3.5, height= 4)
print(reheat1_comp_changes)
dev.off()
# compare study var -------------------------------------------------------

# Plot distributions
significant_cells <- study_diff_stats_lmer %>%
  dplyr::filter(adj_pval <= 0.05) %>%
  pull(cell_type)

dens_plt <- M %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "cell_type") %>%
  left_join(target.df, by = "Sample") %>%
  dplyr::filter(cell_type %in% significant_cells) %>%
  ggplot(aes(x = value, y = HeartFailure, color = HeartFailure)) +
  geom_violin(width=.8, alpha = 0.1, aes(fill = HeartFailure)) +
  geom_boxplot(width=0.15) +
  facet_wrap(. ~ cell_type,
             ncol = 2,
             nrow = 4,
             scales = "free") +
  scale_fill_manual(values = c("black", "darkgrey")) +
  scale_color_manual(values = c("black", "darkgrey")) +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("") +
  xlab("proportions")
dens_plt

study_var <- study_diff_stats_lmer%>%
  ggplot(aes(x = -log10(adj_pval),
             y = perc_studyvar,
             label = cell_type)) +
  geom_vline(xintercept = -log10(0.05)) +
  geom_point() +
  geom_label_repel(size =3) +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12)) +
  ylab("% of var by study")
study_var

###
pdf("output/figures/deconv_lmer_random.pdf", height = 2, width = 5)
 plot(study_var)
dev.off()


# plot compositions as stacked vis -------------------------------------------------------------

res %>% select(-value, -Correlation,-RMSE)
  column_to_rownames("sample_id") %>%
  as.matrix()
  
sample_compositions<-  as.data.frame(M)%>% 
  #select(-P.value, -Correlation,-RMSE)%>%
  as.matrix()

# CT order
ct_order <- hclust(dist(t(sample_compositions)))
ct_order <- colnames(sample_compositions)[ct_order$order]
# Patient order
pat_order <- hclust(dist((sample_compositions)))
pat_order <- rownames(sample_compositions)[pat_order$order]

# Make heatmap of compositions
sample_compositions2<- sample_compositions%>% as.data.frame() %>%
  rownames_to_column("sample_id")%>%# <- read_csv("./results/compositions/sample_comps.csv") %>%
  pivot_longer(-sample_id, names_to = "cell_type", values_to = "props")  %>%
  dplyr::mutate(sample_id = factor(sample_id, levels = pat_order),
                cell_type = factor(cell_type, levels = ct_order))

comps_plot <- ggplot(sample_compositions2, aes(fill=cell_type, y=props, x=sample_id)) +
  geom_bar(position="stack", stat="identity", width= 1.4,colour=NA) +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent)+
  theme(axis.text.y = element_blank(),
        axis.text.x =  element_text(size= 13, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), 'cm'),
        legend.position = "right") +
  scale_fill_manual(values= unname(color_list$ctype_colors[levels(sample_compositions2$cell_type)]))+
  xlab("") +
  ylab("")+
  theme(plot.margin = unit(c(1, 1, 1, 0), "cm"))+
  coord_flip()
comps_plot

target.df
all_meta <- target.df %>%
  filter(Sample %in% sample_compositions2$sample_id)%>%
  mutate(sample_id = factor(Sample, levels = pat_order)) %>%
  dplyr::select(sample_id, HeartFailure,  study)

study_plt <-ggplot(all_meta, aes(fill=study, y=1, x=sample_id)) +
  geom_tile(colour=NA) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), 'cm'),
        legend.position = "none") +
  xlab("") +
  ylab("")+
  scale_fill_manual(values= unname(color_list$study_colors_bulk))
  #scale_fill_manual(values= unname(color_list$study_colors_bulk[all_meta$study%>% unique()]))
  study_plt
hf_plt <- ggplot(all_meta, aes(fill=HeartFailure, x=1, y=sample_id)) +
  geom_tile(colour=NA) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), 'cm'),
        legend.position = "none") +
  xlab("") +
  ylab("") +
  scale_fill_manual(values = c("darkgrey","black"))
hf_plt
plot_grid(study_plt+theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
          hf_plt+theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
          ncol= 1)
stacked_plt <- cowplot::plot_grid(comps_plot,
                                  #study_plt+theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), 
                                  hf_plt+theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
                                  rel_widths = c(1,.5),
                                  ncol = 2, 
                                  align = "hv")
                                  
stacked_plt

pdf("output/figures/comps_reheat1.pdf", 
    height= 5, 
    width = 6)
stacked_plt
dev.off()


#save source data

sample_compositions2 %>%
  left_join( target.df%>%
               select(Sample, HeartFailure, study), 
             by= join_by(sample_id == Sample)) %>%
  save_source_data(data=., fig_type =T, 5, "F")


# calculate silhouettes ---------------------------------------------------
## evaluate via silhouette :
calculate_sw <- function(scores, meta, test_label) {
  meta_info <- meta %>%
    dplyr::select_at(c("Sample", test_label, "study")) %>%
    column_to_rownames("Sample") %>%
    dplyr::mutate(label_ix = .data[[test_label]] %>% 
                    as.factor() %>% 
                    as.integer())
  
  scores <- scores[rownames(meta_info),wanted_cells]
  
  patient_dists <- dist(scores,
                        method = "euclidean")

  si <- silhouette(x = meta_info$label_ix, patient_dists) %>%
    as.data.frame() %>%
    dplyr::select(cluster, sil_width) 
  si<- cbind(si, meta_info)
  #left_join(unique(meta_info), by = c("cluster" = "label_ix"))
  
  return(si)
}

df <- map2(M,sig.matrices_names, function(X, sig){
  sil<- calculate_sw(X, meta, "HeartFailure")%>%
    mutate(sig.matrix= sig)
  })%>% do.call(rbind, .)%>% as_tibble()

ggplot(df, aes(x=sig.matrix, y= sil_width, fill = HeartFailure ))+
  geom_boxplot()

ggplot(df, aes(x=sig.matrix, y= sil_width ))+
    geom_boxplot()
ggplot(df, aes(x=study, y= sil_width, fill = sig.matrix ))+
    geom_boxplot()


# evaluate by comparing distances to real comps -------------------------------------------------------------------------
#
sc<- read_csv("data/sample_comps.csv")
cm<- read_csv("data/comps_meta.csv")
cell_dic <- read.csv("output/celltype_dictionary_HCA_ReHeaT2.csv")

sc.matrix<- sc %>% as.data.frame()%>% column_to_rownames("sample_id")%>%as.matrix()
sig= M
sample.comp.distances<- lapply(M, function(sig){
  A<- sig[, cell_dic$pred]%>% 
    as.matrix()
  #translate 
  colnames(A)<-cell_dic[match(colnames(A), cell_dic$pred),"real"]
  
  A <- rbind(sc.matrix[,colnames(A)], A)
  
  tidydist<-as.matrix(dist(A), method="eucledian")%>%
    as.data.frame()%>% 
    rownames_to_column("sample1")%>% 
    pivot_longer(-sample1, values_to = "dist", names_to = "sample2")%>%
    left_join(cm, by= c("sample1"= "sample_id"))%>% # join for reheat2
    left_join(target.df, by=c(sample2 = "Sample"))  
  
})
names(tidydist)<-names(M)

sample.comp.distances<- map(names(tidydist), function(x){
    sample.comp.distances[[x]] %>% 
    mutate(signature_matrix= x)%>%
    drop_na()
})%>% do.call(rbind,.)

#plot
tidydist%>%
#sample.comp.distances %>%
  #filter(heart_failure== "HF" & HeartFailure=="yes")%>%
  ggplot(aes(x= HeartFailure, y= dist, fill = heart_failure))+
  geom_boxplot()+
  #facet_wrap(~signature_matrix)+
  theme_cowplot()+
  ggtitle("Distance ReHeaT1 & 2")+
  labs(x= "ReHeaT1 bulk", 
       fill = "ReHeaT2 sc")

sample.comp.distances %>%
  #filter(heart_failure== "HF" & HeartFailure=="yes")%>%
  group_by(sample2, heart_failure)%>%
  summarise(m.dist= mean(dist))%>% arrange(sample2)%>%
  ggplot(aes(x= signature_matrix, y= dist, fill = signature_matrix))+
  geom_boxplot()+
  theme_cowplot()+
  ggtitle("Distance ReHeaT1 & 2 HF")

sample.comp.distances %>%
  filter(heart_failure== "NF" & HeartFailure=="no")%>%
  ggplot(aes(x= signature_matrix, y= dist, fill = signature_matrix))+
  geom_boxplot()+
  theme_cowplot()+
  ggtitle("Distance ReHeaT1 & 2 NF")

sample.comp.distances %>%
  #filter(heart_failure== "NF" & HeartFailure=="no")%>%
  filter(heart_failure== "HF" & HeartFailure=="yes")%>%
  group_by(sample2, heart_failure, signature_matrix)%>%
  summarise(m.dist= mean(dist))%>% arrange(sample2)%>%
  ggplot(aes(x= signature_matrix, y= m.dist, fill = heart_failure))+
  geom_boxplot()+
  theme_cowplot()+
  ggtitle("Distance ReHeaT1 & 2 HF")

