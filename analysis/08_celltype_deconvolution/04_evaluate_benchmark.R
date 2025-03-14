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
library(cowplot)
library(rstatix)
library(ComplexHeatmap)
library(circlize)

source("celltype_deconvolution/utils_decon.R")
source("aesthetics.R")
source("make_source_data.R")
color_list <- readRDS("color_list_figures.rds")


#mixtures and proportions
pb= readRDS("data/HF_studiespb.rds")
# get proportions  
sc<- read_csv("data/sample_comps.csv")
cm<- read_csv("data/comps_meta.csv")#

pb.summed<-readRDS("data/val_studies_samplebp.rds")
sc_combined<- read_csv("data/val_studies_composition.csv")
comb_meta<- read_csv("data/val_studies_meta.csv")%>%
  select(study, orig.ident, disease_code, heart_failure)%>%
  rename(sample_id = orig.ident)

#filter hca out 
amrute_nf <- comb_meta %>% filter(study== "Amrute" & heart_failure=="NF")%>% pull(sample_id)%>% unique()
amrute_hf <- comb_meta %>% filter(study== "Amrute" & heart_failure=="HF")%>% pull(sample_id)%>% unique()
cm%>% filter(study=="Koenig2022_DCM")%>% pull(sample_id)
hca_samples <- comb_meta %>% filter(study== "HCA_pb" & heart_failure=="NF")%>% pull(sample_id)%>% unique()

sc_combined<-sc_combined %>% 
  filter(!sample_id %in% c(hca_samples, amrute_nf))

comb_meta<-comb_meta %>% filter(!sample_id %in% c(hca_samples, amrute_nf), study!= "HCA_pb")

cm<-rbind(cm%>% mutate(core="yes"), comb_meta%>% mutate(core ="no"))
sc <- rbind(sc, sc_combined)

sc <- sc %>%
  mutate_all(~ ifelse(is.na(.), 0, .))

mixtures_names <- str_replace_all(list.files("output/deconvo/mixtures_reheat", full.names = F),".txt", "")


#meta data of mixtures
#obs <- read.csv("data/metadata/combined_meta.csv")

#deconvo results:
#res <- readRDS("output/deconvo/cibersort_results.rds")
res <- readRDS("output/deconvo/cibersort_results_markers2.rds")

# align cell type labels --------------------------------------------------

#the labels are slightly different between the mixtures and the reference and 
# will be harmonized via a dictionary_

#reich_meta= read.csv("../MOFAcelldata/processed/reichart/Reichart2022_meta_data.csv")[,-1]%>% mutate(study= "Reichart2022")
#chaff_meta= read.csv("data/processed/scell_chaffin/Chaffin2022_meta_data.csv")[,-1]%>% mutate(study= "Chaffin2022")

#obs %>% as_tibble()

#chaff= pb$Chaffin2022$props %>% column_to_rownames("donor_id")%>% t()
#reich= pb$Reichart2022$props %>% column_to_rownames("donor_id")%>% t()



# evaluation of the real mixtures -------------------------------------------------------------------------

# We have n different signature matrix variations per feature space
# we will first plot the predicted vs real proportions for each of the matrices
# to see if there is a lot of variation due to single cell sampling 
# then we take the mean of the proportions from the n signature matrices which
# we then proceed to calculate two metrics, RMSE and pearson, to compare the matrix 
# performances


# part 1 - plot the variation in signature matrices due to sampling

#res$mix_Chaffin2022

#sig.matrices_names <-str_replace_all(sig.matrices <- list.files("output/deconvo/signature_Ms", full.names = F), ".txt", "")
res$mix_Amrute
sig.matrices_names= c("allGenes", 
                      "comp_", 
                      "mol", "not_reg")
                      #"unknown",
                      #"markers")
# mix= "mix_Chaffin2022"
# sig = "comp_"
# sig.name= names(sig.list)[[1]]
# 
#loop over mixtures
tidy_results <-map(names(res), function(mix){
  
  # loop over signature matrix names
  sig_df <- map(sig.matrices_names, function(sig){
    #extract those 
    print(sig)
    sig.list <- res[[mix]][grepl(names(res[[mix]]), pattern = sig)]
    print(paste0(sig, " and -> ", names(sig.list)))
    
    #bring in tidy format
    sig_df = map(names(sig.list) , function(sig_name){
      as.data.frame(sig.list[[sig_name]])%>%
        rownames_to_column("sample_id")%>%
        select(sample_id, cell_dic$pred)%>%
        mutate(signature_matrix= sig_name)%>%
        as_tibble()
    })%>% do.call(rbind, .)%>%
      mutate(signature_matrix_category = sig)
  })%>%
    do.call(rbind, .) %>%
      pivot_longer(cols = (cell_dic$pred), names_to = "cell_type", values_to = "predicted_prop")%>%
    mutate(mix_name= mix)
  
  #compare the different signature matrices based on the subsampling:
  ggplot(sig_df, aes(x= cell_type, y= predicted_prop, fill= signature_matrix))+
    geom_boxplot()+
    facet_grid(~signature_matrix_category)+
    theme(legend.position = "none", 
          axis.text.x = element_text(angle= 90, hjust= 1))
  
  #take the mean:
  sig_df_summarized= sig_df %>% 
    group_by(signature_matrix_category, sample_id, cell_type, mix_name)%>%
    summarise(mean_predicted_prop= mean(predicted_prop))
  
  return(sig_df_summarized)
})%>% do.call(rbind,.)

tidy_results

update_sig_names <- function(df){
  df$signature_matrix_category<- str_replace_all(df$signature_matrix_category,"allGenes", "dereg")
  df$signature_matrix_category<- str_replace_all(df$signature_matrix_category,"comp_", "dereg_comp")
  df$signature_matrix_category<- str_replace_all(df$signature_matrix_category,"mol", "dereg_mol")
  df$signature_matrix_category <- str_replace_all(df$signature_matrix_category,
                                                  "not_reg","unreg")
  df<- df %>% 
    mutate(signature_matrix_category= factor(signature_matrix_category, 
                                             levels= c("unreg", 
                                                       "dereg", 
                                                       "dereg_mol", 
                                                       "dereg_comp")))
  return(df)
}


tidy_results<- update_sig_names(tidy_results)

#translate cell types via dictionary
tidy_results$cell_type <- 
  cell_dic[match(tidy_results$cell_type, cell_dic$pred),"real"]

#this was before adding more studiesadd real proportions
# real_props= map(names(pb), function(x){
#   pb[[x]][["props"]]%>% 
#     pivot_longer(cols = (cell_dic$real), names_to = "cell_type", values_to = "real_prop")%>%
#     mutate(mix_name= paste0("mix_", x))
# })%>% do.call(rbind,.)%>% 
#   rename(sample_id = donor_id)
sc.matrix <- sc%>%
  distinct()%>%
  as.data.frame()%>% 
  column_to_rownames("sample_id")%>%
  as.matrix()
real_props = sc.matrix %>% 
  as.data.frame()%>% 
  rownames_to_column("sample_id")%>%
  pivot_longer(col= -sample_id, names_to = "cell_type", values_to = "real_prop")

tidy_results<- tidy_results %>% left_join(real_props, by= c("cell_type", "sample_id"))

tidy_results <- tidy_results %>%
  #group_by(signature_matrix_category, sample_id, cell_type)%>%
  rowwise()%>%
  mutate(delta = mean_predicted_prop- real_prop)

saveRDS(tidy_results,"output/deconvo/tidy_results.rds")

tidy_results<- tidy_results %>% filter(!mix_name %in% c("mix_Reichart2022", "mix_HCA_pb"))

samps_amrute_koenig<- tidy_results%>% filter(mix_name %in% c("mix_Amrute", "mix_Koenig2022"))%>% group_by(sample_id)%>%count()%>%
  filter(n==56)%>% pull(sample_id)

tidy_results<- tidy_results %>% filter(!(mix_name == "mix_Amrute" & sample_id %in% samps_amrute_koenig))


# Calculate errors --------------------------------------------------------

## Part 2- calculate errors and do statisical tests
#calculate pearson and rmse, for both patients and cell types:

cell_wise_error <- tidy_results%>%
  group_by(cell_type,mix_name, signature_matrix_category)%>%
  mutate(cell_type_rmse= sqrt(mean(delta^2)))
cell_wise_error

patient_error <- tidy_results%>%
  left_join(cm%>% select(-study)%>% distinct(), by= "sample_id")%>%
  group_by(signature_matrix_category, mix_name)%>%
  mutate(sample_rmse= sqrt(mean(delta^2)))

patient_error_hf <- tidy_results%>%
  left_join(cm%>% select(-study)%>% distinct(), by= "sample_id")%>%
  group_by(signature_matrix_category, mix_name, heart_failure)%>%
  mutate(sample_rmse= sqrt(mean(delta^2)))

cell_corr<- tidy_results%>%
  group_by(cell_type, signature_matrix_category, mix_name)%>%
  summarise(cell_type_corr= cor(mean_predicted_prop, real_prop), .groups = "keep")

pat_cor<- tidy_results%>%
  group_by(signature_matrix_category, mix_name)%>%
  summarise(corr= cor(mean_predicted_prop, real_prop), .groups = "keep")

pat_cor_hf<- tidy_results%>%
  left_join(cm%>% select(-study)%>% distinct(), by= "sample_id")%>%
  group_by(signature_matrix_category, mix_name,heart_failure)%>%
  summarise(corr= cor(mean_predicted_prop, real_prop), .groups = "keep")

#numbers for paper
patient_error_hf %>% 
  group_by(signature_matrix_category, heart_failure)%>%
  summarise(m = median(sample_rmse, na.rm= T))
pat_cor_hf %>% 
  group_by(signature_matrix_category, heart_failure)%>%
  summarise(m = median(corr, na.rm= T))


# Plot errors -------------------------------------------------------------

## for the manuscript we will show that errors differe between sig matrices
## and between hf cohorts 
library(ggpubr)
# correlation plot 

my_comparisons <- list(c("unreg", "dereg_comp"),
                       c("dereg", "dereg_comp"),
                       c("dereg_mol", "dereg_comp")
                       )

#pat_cor_hf<- update_names(pat_cor_hf)
p.study_corr<- pat_cor_hf %>%
  ggplot(.,aes(x= signature_matrix_category, y= corr))+
  facet_grid(~heart_failure)+
  geom_boxplot(alpha= 0.8, outlier.colour = NA)+
  geom_jitter(width= 0.1)+
  theme_cowplot()+
  labs(x= "", 
       y= "Pearson's Correlation")+
  stat_compare_means(comparisons = my_comparisons, 
                     paired = TRUE, method = "wilcox.test"
  )+ # Add pairwise comparisons p-value
  theme(axis.text.x = element_text(angle = 90,vjust= 0.5, hjust= 1))+
#  theme(axis.text.x = element_blank())+#(angle = 90,vjust= 0.5, hjust= 1))+
  geom_hline(yintercept = 1)+
  scale_y_continuous(breaks= c(0.4, 0.6, 0.8, 1))

p.study_corr
# save the soruce data
save_source_data(T, 5, "E", data=pat_cor_hf)
# get Ns from this
pat_cor_hf%>%
  ungroup()%>%
  group_by(signature_matrix_category, heart_failure) %>%
  count()

# rmse plot 


p.study_rmse<- patient_error_hf %>% 
  select(-sample_id,-disease_code)%>%
  distinct(signature_matrix_category, sample_rmse, heart_failure, core)%>%
  ggplot(.,aes(x= signature_matrix_category, y= sample_rmse))+
  facet_grid(~heart_failure)+
  geom_boxplot(alpha= 0.8, outlier.color = NA)+
  geom_jitter(width= 0.1)+
  theme_cowplot()+
  geom_hline(yintercept = 0)+
  labs(x= "", 
       y= "RMSE")+
  stat_compare_means(comparisons = my_comparisons, 
                     paired = TRUE,method = "wilcox.test"
                    )+ # Add pairwise comparisons p-value
  theme(axis.text.x = element_text(angle = 90, vjust= 0.5, hjust= 1))
p.study_rmse
#save source data
patient_error_hf %>% 
  select(-sample_id,-disease_code)%>%
  distinct(signature_matrix_category, sample_rmse, heart_failure, core)%>%
  save_source_data(T, 5, "E2", data=.)

p1= plot_grid(p.study_rmse,
          p.study_corr)

pdf("output/figures/deconv_rmse_corr.pdf", 
    width = 5, height= 5)
  print(p1)
dev.off()

## supp figure
#plot: 

p.cell_rmse <- 
  cell_wise_error %>% 
  distinct(mix_name, signature_matrix_category, cell_type_rmse)%>%
  #filter(mix_name != "mix_Reichart2022")%>% 
  #ggplot(.,aes(x= signature_matrix_category, y= cell_type_rmse))+
  ggplot(.,aes(x= signature_matrix_category, col = signature_matrix_category, 
               y= cell_type_rmse))+
  #geom_col(position = "dodge", col ="black")+
  geom_boxplot(outlier.colour = NA)+
  facet_grid(~cell_type)+
  scale_color_brewer(type="qual",palette = 6)+
  geom_jitter(size= 1, width= 0.1)+
  theme_minimal()+
  labs(col= "Signature matrix", 
       y= "RMSE", 
       x="")+
  theme(axis.text.x = element_blank())

p.cell_rmse

## correlations
p.cell_corr<- cell_corr %>% 
  distinct(mix_name, signature_matrix_category, cell_type_corr)%>%
  ggplot(.,aes(x= signature_matrix_category, y= cell_type_corr, 
               col= signature_matrix_category))+
  geom_boxplot(outlier.colour = NA)+
  facet_grid(~cell_type)+
  scale_color_brewer(type="qual",palette = 6)+
  geom_jitter(size= 1, width= 0.1)+
  theme_minimal()+
  labs(col= "Signature matrix", 
       y= "Correlation", 
       x="")+
  theme(axis.text.x = element_blank())
p.cell_corr

p.cell_performance<- plot_grid(p.cell_rmse,
                               p.cell_corr+theme(legend.position = "none"), ncol = 1, 
                               align= "v")

pdf("output/figures/deconvo_performance_cells.pdf",
    width = 7, height= 4)
p.cell_performance
dev.off()

# plot compositions! 

p.compositions<-tidy_results%>%
  select(signature_matrix_category,  sample_id, cell_type, mix_name, mean_predicted_prop, real_prop)%>%
  drop_na()%>%
  filter(mix_name != "mix_Reichart2022")%>% 
  ggplot(., aes( x= real_prop, y= mean_predicted_prop, color= cell_type))+
  geom_point(alpha= 0.9)+
  facet_wrap(~signature_matrix_category, ncol=1)+
  #geom_smooth(method= "lm")+
  geom_abline(slope= 1)+
  scale_colour_manual(values= unname(color_list$ctype_colors))+
  xlim(c(0, 1))+
  ylim(c(0,1))+
  theme_cowplot()+
  coord_equal()+
  labs(x= "true composition", 
       y= "predicted composition")
p.compositions

# instead of points plot linear models 
p.compositions2<-
  tidy_results%>%
  drop_na()%>%
  ggplot(., aes( x= real_prop, y= mean_predicted_prop, color= cell_type, shape= mix_name))+
  #geom_point(alpha= 0.4)+
  facet_wrap(~signature_matrix_category, ncol=1)+
  geom_smooth(method= "lm")+
  scale_colour_manual(values= unname(color_list$ctype_colors))+
  geom_abline(slope= 1)+
  xlim(c(0, 1))+
  ylim(c(0,1))+
  theme_cowplot()+
  coord_equal()+
  labs(x= "true composition", 
       y= "predicted composition")


pdf("output/figures/deconvo_sc_bench_compositions.pdf", 
    width= 6, height= 8)
p.compositions
p.compositions2
dev.off()



# do stats: ---------------------------------------------------------------

# do stats test 
pat_cor_hf%>% 
  select(-mix_name)%>% 
  ungroup()%>%
  group_by(heart_failure)%>%
  pairwise_wilcox_test(as.formula("corr ~ signature_matrix_category"), paired = TRUE,
                  p.adjust.method = "BH"
  )

patient_error_hf%>% 
  select(-mix_name)%>% 
  ungroup()%>%
  pairwise_t_test(as.formula("sample_rmse ~ signature_matrix_category"), paired = TRUE,
                  p.adjust.method = "BH"
  )

pat_cor%>% select(-mix_name)%>% ungroup()%>%
  pairwise_t_test(as.formula("corr ~ signature_matrix_category"), paired = TRUE,
    p.adjust.method = "BH"
  )


patient_error%>% ungroup()%>%drop_na()%>%
  #group_by(heart_failure)%>%
  pairwise_t_test(as.formula("sample_rmse ~ signature_matrix_category"), paired = TRUE,
                       p.adjust.method = "BH"
  )

patient_error_hf%>% ungroup()%>%drop_na()%>%
  group_by(signature_matrix_category)%>%
  pairwise_t_test(as.formula("sample_rmse ~ heart_failure"), paired = TRUE,
                  p.adjust.method = "BH"
  )

df.cell_corr <- cell_corr %>% 
  distinct(cell_type, signature_matrix_category, cell_type_corr)%>% 
  group_by(cell_type)%>%
  pairwise_t_test(as.formula("cell_type_corr ~ signature_matrix_category"), paired = TRUE,
                  p.adjust.method = "BH"
  )

df.cell_corr%>% 
  filter(group1=="dereg_mol", group2=="dereg_comp")
  filter(group1=="unreg", group2=="dereg_comp")
  
df.cell_rmse <- cell_wise_error %>% 
    distinct(cell_type, signature_matrix_category, cell_type_rmse)%>% 
    group_by(cell_type)%>%
    pairwise_t_test(as.formula("cell_type_rmse ~ signature_matrix_category"), paired = TRUE,
                    p.adjust.method = "BH"
    )

df.cell_rmse%>% 
  filter(group1=="dereg_mol", group2=="dereg_comp")

# evaluate on NF-HF contrast level ----------------------------------------
# now we will calculate a composition change and see if the different signature matrices can recover this
#change 

#join with meta data:
tidy_results= tidy_results %>% 
  left_join(cm%>% select(-study)%>% distinct(), by= "sample_id")
  #left_join(obs, by= "sample_id")


prop_data <- tidy_results %>%
  dplyr::select(sample_id, cell_type, real_prop,mean_predicted_prop, signature_matrix_category, mix_name)

prop_data_pred <- prop_data %>%
  pivot_wider(id_cols = -real_prop, values_from = mean_predicted_prop, names_from = cell_type, values_fill = 0) 

prop_data_real <- prop_data %>%
  select(-signature_matrix_category, -mean_predicted_prop)%>%
  distinct()%>%
  pivot_wider(values_from = real_prop, names_from = cell_type, values_fill = 0)
  #column_to_rownames("sample_id") %>%
  #as.matrix()
sig= "allGenes"
mix= "mix_CHD_HILL_pb"
#studies with only one condition need to be excluded 
exclude<- c("mix_Amrute", "mix_HCM_dimmeler_pb", "mix_Sarcoidosis_pb")

prop_m<-  lapply(unique(prop_data_pred$signature_matrix_category), function(sig){
  print(sig)
  lapply(unique(prop_data_pred$mix_name), function(mix){
    print(mix)
    if(mix %in% exclude){return(NULL)}
    prop_data_clr<- prop_data_pred%>% 
      filter(signature_matrix_category== sig, 
             mix_name== mix)%>%
      select(-signature_matrix_category,-mix_name)%>%
      column_to_rownames("sample_id") %>%
      as.matrix()%>%
      compositions::acomp() %>%
      compositions::clr()
    
    study_diff_stats <- prop_data_clr %>%
      as.data.frame() %>%
      rownames_to_column("sample_id") %>%
      pivot_longer(-sample_id, names_to = "cell_type") %>%
      left_join(cm%>% select(-study)%>% distinct(), by= "sample_id") %>%
      group_by( cell_type) %>%
      nest() %>%
      dplyr::mutate(HF_diff = map(data, function(dat) {
        
        t.test(value ~ heart_failure, data = dat) %>%
          broom::tidy()
        
      })) %>%
      dplyr::select(-data) %>%
      unnest(HF_diff) %>%
      mutate(signature_matrix_category= sig, 
             mix_name= mix)
  })%>% do.call(rbind,.)
})%>% do.call(rbind,.)

real_m <-  lapply(unique(prop_data_real$mix_name), function(mix){
  print(mix)
  if(mix %in% exclude){return(NULL)}
    prop_data_clr<- prop_data_real%>% 
      filter( mix_name== mix)%>%
      select(-mix_name, -signature_matrix_category)%>%
      distinct()%>%
      column_to_rownames("sample_id") %>%
      as.matrix()%>%
      compositions::acomp() %>%
      compositions::clr()
    
    study_diff_stats <- prop_data_clr %>%
      as.data.frame() %>%
      rownames_to_column("sample_id") %>%
      pivot_longer(-sample_id, names_to = "cell_type") %>%
      left_join(cm%>% select(-study)%>% distinct(), by = "sample_id") %>%
      group_by( cell_type) %>%
      nest() %>%
      dplyr::mutate(HF_diff = map(data, function(dat) {
        
        t.test(value ~ heart_failure, data = dat) %>%
          broom::tidy()
        
      })) %>%
      dplyr::select(-data) %>%
      unnest(HF_diff) %>%
      mutate(signature_matrix_category= "real", 
             mix_name= mix)
  })%>% do.call(rbind,.)

# combine :
prop.df <- rbind(real_m,prop_m) %>% 
  filter(mix_name != "mix_Reichart2022")

prop.df%>% 
  ggplot(., aes(x= cell_type, y= statistic, color = signature_matrix_category ))+
  geom_point()+
  facet_wrap(~mix_name)

map(unique(prop.df$mix_name), function(mix){
  print(mix)
  prop.df%>%
    filter(mix_name== mix)%>%
    select(cell_type, statistic, signature_matrix_category)%>%
    pivot_wider(values_from = statistic, names_from=cell_type)%>%
    as.data.frame()%>%
    column_to_rownames("signature_matrix_category")%>%
    Heatmap(name= mix)
})
mix= "mix_Chaffin2022"
## calculate distance of the change vectors between real and the sig matrices per study
# and plot them together
dists_change<- map(unique(prop.df$mix_name), function(mix){
  print(mix)
  M= prop.df%>%
    filter(mix_name== mix)%>%
    select(cell_type, estimate, signature_matrix_category)%>%
    pivot_wider(values_from = estimate, names_from=cell_type)%>%
    as.data.frame()%>%
    column_to_rownames("signature_matrix_category")%>%
    as.matrix()
  
  dists<- dist(M)%>% as.matrix()%>% as.data.frame()%>% select(real)
  print(dists)
  return(t(dists))
})%>% do.call(rbind,.)

#rmse instead of eucledian
dists_change<- map(unique(prop.df$mix_name), function(mix){
  print(mix)
  M= prop.df%>%
    filter(mix_name== mix)%>%
    select(cell_type, estimate, signature_matrix_category)%>%
    pivot_wider(values_from = estimate, names_from=cell_type)%>%
    as.data.frame()%>%
    column_to_rownames("signature_matrix_category")%>%
    as.matrix()
  
  M2<- M - M[rep(1, nrow(M) ), ]
  rmse_M<- sqrt(rowMeans(M2^2))
  
  print(dists)
  return(rmse_M)
})%>% do.call(rbind,.)

rownames(dists_change)<- unique(prop.df$mix_name)
dists_change%>% 
  Heatmap()

dists_change%>% colMeans()%>%
  Heatmap()

dists_change_tidy<- 
  dists_change%>% as.data.frame()%>% as_tibble()%>%
  mutate(mix_name= unique(prop.df$mix_name))%>%
  pivot_longer(-mix_name, names_to= "signature_matrix", values_to = "change_distance")%>%
  filter(signature_matrix != "real")

my_comparisons <- list(c("allGenes", "mol"),
                       c("compMol", "mol"), 
                       c("allGenes", "compMol"),
                       c("comp_", "mol"))

p.dists<-
dists_change_tidy%>%
  ggplot(aes(x= signature_matrix, y= change_distance))+
  geom_boxplot()+
  geom_point(aes(color= mix_name))+
  theme(axis.text.x =element_text(angle= 45, hjust= 1))+
  stat_compare_means(comparisons = my_comparisons, paired = T)
p.dists
#stats test
dists_change_tidy%>% 
  pairwise_t_test(change_distance ~ signature_matrix, paired = TRUE )
dists_change_tidy%>% 
  pairwise_wilcox_test(change_distance ~ signature_matrix, paired = TRUE )


# calculate F1 score  -----------------------------------------------------
prop.df$mix_name

df.f1 = prop.df%>% #filter(signature_matrix_category=="real")%>% print(n=100)%>%
  #filter(!mix_name %in% c("mix_Kuppe"))%>%
  mutate(prediction= ifelse(p.value < 0.05 & estimate >0 , 
                            "increase", 
                            ifelse(p.value < 0.05 & estimate <0 , "decrease", 
                                   "no_change")))

df.f1 <- df.f1%>% filter(signature_matrix_category=="real")%>% 
  select(cell_type, mix_name, prediction)%>%
  left_join(df.f1 %>%rename(prediction2 = prediction)%>%
              select(cell_type, mix_name, signature_matrix_category, prediction2),
            by= c("cell_type", "mix_name"))%>%
  ungroup()

#  filter(signature_matrix_category != "real")

  df.f1 %>%
  print(n=100)

df.f1 <- df.f1 %>% group_by(signature_matrix_category)%>%
  mutate(x= (prediction==prediction2))
table(df.f1$prediction, df.f1$prediction2)

# Example contingency matrix

# Function to calculate F1 Score for each class
calculate_f1 <- function(conf_matrix) {
  f1_scores <- numeric(nrow(conf_matrix))
  for (i in 1:nrow(conf_matrix)) {
    tp <- conf_matrix[i, i]  # True Positives for class i
    fp <- sum(conf_matrix[-i, i])  # False Positives for class i
    fn <- sum(conf_matrix[i, -i])  # False Negatives for class i
    
    precision <- tp / (tp + fp)
    recall <- tp / (tp + fn)
    
    f1_scores[i] <- ifelse((precision + recall) == 0, 0, 2 * (precision * recall) / (precision + recall))
  }
  names(f1_scores) <- rownames(conf_matrix)
  return(f1_scores)
}
library(colorspace)
f1= sapply(unique(df.f1$signature_matrix_category), function(y){
         print(y)
         df<-df.f1%>% 
           ungroup()%>%
           filter(signature_matrix_category== y)
         print(df)
         confusion_m <- table(df$prediction, df$prediction2)
         print(confusion_m)
         f1_scores <- calculate_f1(confusion_m)
       })
f1

heatmap_color_fun <- colorRamp2(
  c( min(f1),max(f1[,-1])),  # Define the range of the data
  c( "white", color_list$general_up_down_colors["up"])
)
library(circlize)
column_ha = HeatmapAnnotation(mean_F1 = anno_barplot(colMeans(f1[,-1])),
                              col= list(bar1= "black"))
p.f1.hmap = Heatmap(f1[,-1], name ="F1-score",
        rect_gp = gpar(col = "black", lwd = 1),
        top_annotation = column_ha,
        col = heatmap_color_fun,
        show_row_dend = F, 
        show_column_dend = F
       
        )
p.f1.hmap
colMeans(f1)
pdf("output/figures/deconvolute_change_f1.pdf", 
    width = 3, 
    height= 2.5)
p.f1.hmap
dev.off()
# print plots -------------------------------------------------------------



theme_set(theme_cowplot())
theme_update(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
pdf("output/figures/deconvolute_change_dists.pdf", 
    width = 4, 
    height= 4)
p.dists
dev.off()



# sandbox -----------------------------------------------------------------
real_m%>% rename(real_statistic= statistic)%>%
  select(cell_type, real_statistic, mix_name)%>%
  right_join(prop_m, by= c("cell_type", "mix_name"))%>%
  group_by(mix_name, signature_matrix_category, cell_type)%>%
  mutate(delta_t = abs(abs(statistic)- abs(real_statistic)))%>%
  ggplot(. , aes(x= signature_matrix_category ,y= delta_t))+
  geom_boxplot()+
  geom_point(aes(color= cell_type))+
  facet_grid(~mix_name)
# Testing differences in HF
prop_data_clr <- prop_data %>%
  compositions::acomp() %>%
  compositions::clr()

# Compare the vector of comparisons

study_diff_stats <- prop_data_clr %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(-sample_id, names_to = "cell_type") %>%
  left_join(meta_data, by = "sample_id") %>%
  #dplyr::filter(!disease_code %in% test_exclusion) %>%
  group_by(batch, cell_type) %>%
  nest() %>%
  dplyr::mutate(HF_diff = map(data, function(dat) {
    
    t.test(value ~ heart_failure, data = dat) %>%
      broom::tidy()
    
  })) %>%
  dplyr::select(-data) %>%
  unnest(HF_diff) 

t_mat= study_diff_stats%>%
  dplyr::select(cell_type, batch, statistic) %>%
  pivot_wider(names_from = batch, values_from = statistic) %>%
  column_to_rownames("cell_type") %>%
  as.matrix()



