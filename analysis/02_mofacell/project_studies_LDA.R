# Copyright (c) [2023] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we project each study to the multicellular space of
#' any other study.
#'
#' Additionally, we test if the multicellular space classifies NF samples
#' in the projected data
library(MOFAcellulaR)
library(tidyverse)
library(MASS)
library(ROCR)
source("./code/reheat2_pilot/aesthetics.R")

sample_class_plt_file <- "./results/models/study_class_LDA.pdf"
summary_class_plt_file <- "./results/models/study_class_LDA_summ.pdf"
sample_class_file <- "./results/models/study_class_LDA.csv"
sample_projection_plt <- "./results/models/projections_LDA_UMAPs.pdf"
sample_class_plt_file_main <- "./results/models/study_class_LDA_summ_main.pdf"

# These are the training objects
mofa_models_obj <-  tibble(file = list.files("./data/pbulk/")) %>%
  dplyr::mutate(study = gsub("_pbulk[.]csv","",file),
                mofa_out = paste0("./data/models/", gsub("_pbulk[.]csv","_mofa.hdf5",file))) %>%
  dplyr::select(-file) %>%
  dplyr::mutate(mofa_model = map(mofa_out, MOFA2::load_model)) %>%
  dplyr::select(study, mofa_model)

train_obj <- readRDS("./results/LDA/model_stats.rds") %>%
  dplyr::rename("study" = "name") %>%
  left_join(mofa_models_obj, by = "study")

# These are the test ones

test_obj <- tibble(file = list.files("./data/pbulk/")) %>%
  dplyr::mutate(study = gsub("_pbulk[.]csv","",file),
                coldata_file =  paste0("./data/coldata/", gsub("_pbulk[.]csv","_coldata.csv",file)),
                metadata_file =  paste0("./data/metadata_ext/", gsub("_pbulk[.]csv","_metadata.csv",file))) %>%
  dplyr::mutate(file = paste0("./data/pbulk/",file))

model_outs <- pmap(test_obj, function(study, file, coldata_file,  metadata_file) {

  print(study)

  test_study <- study

  # Importing pb data
  pb_data <- read_csv(file,
                      show_col_types = FALSE)

  colnames(pb_data)[1] <- "sample_id"

  pb_data <- pb_data %>%
    column_to_rownames("sample_id") %>%
    as.matrix() %>%
    t()

  # Importing coldata of the matrices - ignoring not annotated cell types
  coldat <- read_csv(coldata_file,
                     show_col_types = FALSE)[,-1]  %>%
    column_to_rownames("colname") %>%
    dplyr::rename(ncells = "counts") %>%
    dplyr::filter(cell_type != "none")

  pb_data <- pb_data[,rownames(coldat)]

  # Importing meta_data

  meta_data <- read_csv(metadata_file,
                        show_col_types = FALSE)

  # Defining cts
  cts <- coldat$cell_type %>%
    unique() %>%
    set_names()

  # Defining parameters
  n_samples <- nrow(meta_data)
  min_samples <- (n_samples * 0.4) %>% floor()

  # Creating MOFAcell object to make projection
  MOFAcell_obj <- MOFAcellulaR::create_init_exp(pb_data, coldat) %>%
    MOFAcellulaR::filt_profiles(pb_dat = .,
                                cts = cts,
                                ncells = 20,
                                counts_col = "ncells",
                                ct_col = "cell_type") %>%
    MOFAcellulaR::tmm_trns(pb_dat_list = .,
                           scale_factor = 1000000) %>%
    MOFAcellulaR::center_views(pb_dat_list = .) %>%
    MOFAcellulaR::pb_dat2MOFA(pb_dat_list = .,
                              sample_column = "sample_id")

  # Here we test across all studies for all cell_types

  pred_obj <- map(set_names(c("All",cts)), function(ct){

    study_eval_df <- train_obj %>%
      #dplyr::filter(study != test_study) %>%
      dplyr::mutate(projected_data = map(mofa_model, function(mf) {

        if(ct == "All") {

          test_projection <- project_data(model = mf,
                                          test_data = MOFAcell_obj)

        } else {

          test_projection <- project_data(model = mf,
                                          test_data = MOFAcell_obj %>%
                                            dplyr::filter(view == ct))

        }

        return(test_projection)
      })) %>%
      dplyr::mutate(AUPRC = map2(value, projected_data, function(LDAm, pd){

        RF_pd <- pd %>%
          as.data.frame() %>%
          tibble::rownames_to_column("sample_id") %>%
          left_join(meta_data[, c("sample_id", "heart_failure")], by = "sample_id") %>%
          dplyr::select(-sample_id) %>%
          dplyr::mutate(HF = ifelse(heart_failure == "NF", "no", "yes")) %>%
          dplyr::mutate(HF = factor(HF)) %>%
          dplyr::select(-heart_failure)

        real_labels <- RF_pd$HF

        pred_labels <- stats::predict(LDAm, RF_pd %>% dplyr::select(-HF))

        pred_labels <- pred_labels$posterior

        eval <- tibble(real_labels,
                       prob_healthy = pred_labels[,"no"]) %>%
          arrange(desc(prob_healthy)) %>%
          dplyr::mutate(real_labels = ifelse(real_labels == "no", 1, 0))

        PR_object <- ROCR::prediction(eval$prob_healthy,
                                      eval$real_labels) #Evaluate classification

        AUC_pr <- ROCR::performance(PR_object, measure = "auc")

        return(AUC_pr@y.values[[1]])

      })) %>%
      dplyr::select(study, projected_data, AUPRC) %>%
      dplyr::mutate("test_study" = test_study) %>%
      dplyr::rename("train_study" = study) %>%
      mutate(AUPRC = unlist(AUPRC))

    return(study_eval_df)

  })

  return(enframe(pred_obj, name = "cell_type") %>%
           unnest(c(value)))
})

study_pred <- model_outs %>%
  enframe() %>%
  dplyr::select(-name) %>%
  unnest(c(value)) %>%
  dplyr::select(-projected_data) %>%
  dplyr::rename("AUROC" = "AUPRC") %>%
  dplyr::filter(train_study != "Reichart2022_DCM",
                test_study != "Reichart2022_DCM")

write_csv(study_pred, sample_class_file)

#How good are studies to predict
training_plt <- study_pred %>%
  group_by(train_study, cell_type) %>%
  dplyr::filter(test_study != train_study) %>%
  summarize(mean_AUROC = mean(AUROC)) %>%
  ggplot(aes(x = cell_type, y = mean_AUROC)) +
  geom_boxplot() +
  geom_point(aes(color = train_study)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggtitle("Training performance") +
  xlab("")

#How good are studies at being predicted
testing_plt <- study_pred %>%
  group_by(test_study, cell_type) %>%
  dplyr::filter(test_study != train_study) %>%
  summarize(mean_AUROC = mean(AUROC)) %>%
  ggplot(aes(x = cell_type, y = mean_AUROC)) +
  geom_boxplot() +
  geom_point(aes(color = test_study)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggtitle("Testing performance") +
  xlab("")

sample_class_plt <- study_pred %>%
  dplyr::filter(cell_type != "All") %>%
  ggplot(aes(x = test_study, y = train_study, fill = AUROC)) +
  geom_tile(color = "black") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom") +
  scale_fill_gradient2(low = "white", midpoint = 0.5) +
  coord_equal() +
  facet_wrap(.~ cell_type, nrow = 1)

pdf(sample_class_plt_file, height = 5, width = 10)

plot(sample_class_plt)

dev.off()

pdf(summary_class_plt_file, height = 5, width = 4)

cowplot::plot_grid(testing_plt, training_plt, nrow = 2, align = "hv")

dev.off()


# For main figures

sample_class_plt_main <- study_pred %>%
  dplyr::filter(cell_type == "All",
                train_study != "Reichart2022_DCM",
                test_study != "Reichart2022_DCM") %>%
  dplyr::mutate(train_study = strsplit(train_study,"_") %>%
                  map_chr(., ~.x[[1]])) %>%
  dplyr::mutate(test_study = strsplit(test_study,"_") %>%
                  map_chr(., ~.x[[1]])) %>%
  ggplot(aes(x = test_study, y = train_study, fill = AUROC)) +
  geom_tile(color = "black") +
  theme_minimal() +
  theme(axis.title = element_text(size = 13),
    axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11, angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_gradient2(low = "white", midpoint = 0.5, limits = c(0.5,1)) +
  coord_equal() +
  xlab("Test") +
  ylab("Train")

pdf(sample_class_plt_file_main, height = 3.1, width = 3.1)

plot(sample_class_plt_main)

dev.off()



# Generate umaps of joint analysis
# First get joint meta-data
set.seed(1908)

all_meta <- tibble(file = list.files("./data/metadata_ext/")) %>%
  dplyr::mutate(study = gsub("_metadata[.]csv","",file)) %>%
  dplyr::filter(study != "Reichart2022_DCM") %>%
  dplyr::mutate(file = paste0("./data/metadata_ext/",file)) %>%
  dplyr::mutate(metadata = map(file, ~ read_csv(.x, show_col_types = FALSE))) %>%
  unnest(metadata) %>%
  dplyr::select(-file)

sample_projection_dat <- model_outs %>%
  enframe() %>%
  dplyr::select(-name) %>%
  unnest(c(value)) %>%
  dplyr::filter(train_study != "Reichart2022_DCM",
                test_study != "Reichart2022_DCM",
                cell_type == "All") %>%
  dplyr::select(-AUPRC) %>%
  group_by(train_study) %>%
  nest() %>%
  dplyr::mutate(UMAP_dat = map2(train_study, data, function(train_s, dat) {

    projected_dat <- bind_rows(map(dat$projected_data, as.data.frame))

    umap_dat <- as.matrix(projected_dat) %>%
      uwot::umap(.,n_neighbors = 30) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("sample_id") %>%
      left_join(all_meta, by = "sample_id")

    return(umap_dat)

  })) %>%
  dplyr::select(train_study, UMAP_dat)


all_scaters <- sample_projection_dat %>%
  unnest(c(UMAP_dat)) %>%
  mutate(type_sample = ifelse(train_study == study, "original", "projected")) %>%
  dplyr::mutate(train_study = strsplit(train_study,"_") %>%
                  map_chr(., ~.x[[1]])) %>%
  ggplot(aes(x = V1, y = V2, color = disease_code, shape = type_sample)) +
  geom_point(size = 2.7, alpha = 0.7) +
  theme_minimal(base_size = 15) +
  # Add a border around each facet
  theme(panel.border=element_rect(fill=NA, colour="grey40"),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size =12),
        legend.title = element_blank(), legend.position = "bottom") +
  facet_wrap(. ~ train_study, ncol = 4) +
  scale_color_manual(values = etiology_colors) +
  xlab("UMAP1") +
  ylab("UMAP2")

pdf(sample_projection_plt, height = 3.3, width = 8)
plot(all_scaters)
dev.off()

# Paper info

study_pred %>%
  dplyr::filter(train_study != "Reichart2022_DCM",
                test_study != "Reichart2022_DCM",
                train_study != test_study,
                cell_type == "All") %>%
  pull(AUROC) %>%
  mean()

study_pred %>%
  dplyr::filter(train_study != "Reichart2022_DCM",
                test_study != "Reichart2022_DCM",
                train_study != test_study,
                cell_type != "All") %>%
  group_by(cell_type) %>%
  summarize(mean_AUROC = mean(AUROC)) %>%
  dplyr::arrange(desc(mean_AUROC))



sample_projection_dat %>%
  unnest(c(UMAP_dat)) %>%
  mutate(type_sample = ifelse(train_study == study, "original", "projected")) %>%
  dplyr::mutate(train_study = strsplit(train_study,"_"))  %>%
                  write_csv("./Revision/figures/Figure2/Figure2G.csv")




study_pred %>%
  dplyr::filter(cell_type == "All",
                train_study != "Reichart2022_DCM",
                test_study != "Reichart2022_DCM") %>%
  dplyr::mutate(train_study = strsplit(train_study,"_") %>%
                  map_chr(., ~.x[[1]])) %>%
  dplyr::mutate(test_study = strsplit(test_study,"_") %>%
                  map_chr(., ~.x[[1]])) %>%
  write_csv("./Revision/figures/Figure2/Figure2H.csv")


