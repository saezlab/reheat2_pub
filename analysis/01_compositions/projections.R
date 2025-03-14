# Copyright (c) [2023] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we calculate the compositional similarities
#' across studies

library(tidyverse)
library(compositions)
library(ROCR)

t_file <- "./results/compositions/t_res.csv"
comps_file <- "./results/compositions/sample_comps.csv"
meta_file <- "./results/compositions/comps_meta.csv"

meta_data_comps <- read_csv(meta_file) %>%
  dplyr::mutate(study = strsplit(study,"_") %>%
                  map_chr(., ~.x[[1]]))

t_vects <- read_csv(t_file) %>%
  column_to_rownames("cell_type") %>%
  as.matrix()

t_vects <- t_vects * -1 # To make healthy the top

comps <- read_csv(comps_file) %>%
  left_join(dplyr::select(meta_data_comps, sample_id, study)) %>%
  group_by(study) %>%
  nest() %>%
  dplyr::mutate(comp_proj = map(data, function(dat) {

    comp_mat <- dat %>%
      column_to_rownames("sample_id") %>%
      as.matrix() %>%
      clr() %>%
      scale()

    comp_mat <- comp_mat[, rownames(t_vects)]

    proj <- comp_mat %*% t_vects %>%
      as.data.frame() %>%
      rownames_to_column("sample_id") %>%
      pivot_longer(-sample_id, names_to = "train_study", values_to = "score") %>%
      left_join(meta_data_comps[,c("sample_id", "heart_failure")], by = "sample_id") %>%
      na.omit()
  })) %>%
  dplyr::select(study, comp_proj) %>%
  dplyr::mutate(performance = map(comp_proj, function(dat) {

    eval_per_study <- dat %>%
      group_by(train_study) %>%
      nest() %>%
      mutate(AUROC = map(data, function(s_dat){

        eval_dat <- s_dat %>%
          arrange(desc(score)) %>%
          dplyr::mutate(healthy = ifelse(heart_failure == "NF", 1, 0))

        PR_object <- ROCR::prediction(eval_dat$score,
                                      eval_dat$healthy) #Evaluate classification

        AUC_pr <- ROCR::performance(PR_object, measure = "auc")

        return(AUC_pr@y.values[[1]])

      }))

    return(eval_per_study)

  }))

# Make plots of the projections

plt_df <- comps %>%
  dplyr::select(study, performance) %>%
  unnest() %>%
  dplyr::select(-data) %>%
  unnest()


#How good are studies to predict
training_plt <- plt_df %>%
  dplyr::filter(study != train_study) %>%
  ggplot(aes(x = train_study, y = AUROC)) +
  geom_boxplot() +
  geom_point(aes(color = study)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggtitle("Training performance") +
  xlab("")

#How good are studies at being predicted
testing_plt <- plt_df %>%
  dplyr::filter(study != train_study) %>%
  ggplot(aes(x = study, y = AUROC)) +
  geom_boxplot() +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggtitle("Testing performance") +
  xlab("")

tile_proj <- plt_df %>%
  unnest() %>%
  ggplot(aes(x = study, y = train_study, fill = AUROC)) +
  geom_tile(color = "black") +
  theme_minimal() +
  theme(axis.title = element_text(size = 13),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11, angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_gradient2(low = "white", midpoint = 0.5, limits = c(0.5,1)) +
  coord_equal() +
  ylab("Train") +
  xlab("Test")

pdf("./results/compositions/projections.pdf", height = 3.1, width = 3.1)

plot(tile_proj)

dev.off()

plt_df %>%
  ungroup() %>%
  dplyr::filter(study != train_study) %>%
  summarize(mean(AUROC))

plt_df %>%
  dplyr::filter(study != train_study) %>%
  group_by(study) %>%
  summarize(mean(AUROC))

plt_df %>%
  dplyr::filter(study != train_study) %>%
  group_by(train_study) %>%
  summarize(mean(AUROC))































