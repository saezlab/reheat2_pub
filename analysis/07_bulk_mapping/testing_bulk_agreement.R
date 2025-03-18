## MOFA cell for deconvolution
library(tidyverse)
library(MOFAcellulaR)
library(decoupleR)
library(cluster)
library(ROCR)

model_outfile <- "./results/meta_mofacell/metamodel_mofa.hdf5"
model <- MOFA2::load_model(model_outfile)
F_analyze <- c("Factor1", "Factor2")

# Calculate disease scores for all samples ----------------------------------
reheat <- readRDS("./data/reheat/METAheart2023.rds")

bulk_mats <- map(set_names(F_analyze), function(f_name) {

  print(f_name)

  # Get loadings
  factor_loadings <- MOFAcellulaR::get_geneweights(model = model,
                                                   factor = f_name)

  #factor_loadings <- factor_loadings %>%
  #  group_by(feature) %>%
  #  dplyr::mutate(n_times = n()) %>%
  #  dplyr::filter(n_times == 1) %>%
  #  dplyr::select(-n_times)

  #factor_loadings <- factor_loadings


  # Calculate disease scores for all samples ----------------------------------
  reheat <- readRDS("./data/reheat/METAheart2023.rds")

  dr_res <- map(reheat, function(x) {

    mat <- x$GEX %>%
      t() %>%
      scale() %>%
      t()

    dR_run <- decoupleR::run_wmean(mat = mat,
                                   network = factor_loadings,
                                   .source = ctype,
                                   .target = feature,
                                   .mor = "value") %>%
      dplyr::filter(statistic == "norm_wmean")

    dR_run <- dR_run  %>%
      left_join(x$TARGETS %>%
                  dplyr::select(Sample, HeartFailure),
                by = c("condition" = "Sample"))

  }) %>%
    enframe() %>%
    unnest()

  ## Make matrix with all scores

  sample_info <- dr_res %>%
    dplyr::select(condition, HeartFailure, name) %>%
    unique() %>%
    dplyr::mutate(sample_id = paste0(name, "_", condition)) %>%
    column_to_rownames("sample_id")

  ds_mat <- dr_res %>%
    dplyr::mutate(sample_id = paste0(name, "_", condition)) %>%
    dplyr::select(sample_id, source, score) %>%
    pivot_wider(names_from = source, values_from = score) %>%
    column_to_rownames("sample_id")

  ds_mat <- as.matrix(ds_mat) %>%
    scale()

  bulk_map <- list("sample_info" = sample_info,
                   "ct_score" = ds_mat)

  return(bulk_map)
})

## For each factor and cell-type test classification

bulk_mats <- map(bulk_mats, function(f_info) {

  sample_anns <- f_info$sample_info %>%
    as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    dplyr::select(-condition) %>%
    dplyr::mutate(response = ifelse(HeartFailure == "no", 1, 0))

  es_scores <- f_info$ct_score %>%
    as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    pivot_longer(-sample_id, names_to = "cell_type", values_to = "ES") %>%
    left_join(sample_anns) %>%
    arrange(name, desc(ES)) %>%
    group_by(name, cell_type) %>%
    nest() %>%
    dplyr::mutate(AUPRC_vals = map(data, function(dat) {

      PR_object <- ROCR::prediction(dat$ES,
                                    dat$response) #Evaluate classification

      AUC_pr <- ROCR::performance(PR_object, measure = "auc") #Change to ROC

      return(AUC_pr@y.values[[1]])

    })) %>%
    dplyr::select(name, cell_type, AUPRC_vals) %>%
    unnest()

  f_info[["AUPRCs"]] <- es_scores %>%
    pivot_wider(values_from = AUPRC_vals, names_from = cell_type) %>%
    column_to_rownames("name") %>%
    as.matrix()

  return(f_info)

})

## Plot heatmaps

cts <- colnames(bulk_mats$Factor1$ct_score)

lims <- max(c(max(bulk_mats$Factor1$ct_score),
            max(bulk_mats$Factor2$ct_score))) + 0.5

color_fun <- colorRamp2(c(lims * -1, 0, lims), c("darkblue", "white", "#d62728"))

# Factor 1 ----------------------------------------------------------------------

ra <- rowAnnotation(HF = bulk_mats$Factor1$sample_info[rownames(bulk_mats$Factor1$ct_score),
                                                           "HeartFailure"],
                        col = list(HF = c("yes" = "black", "no" = "darkgrey")),
                    gap = unit(2.5, "mm"),
                    border = TRUE)

ca <-  HeatmapAnnotation(AUROC =  anno_boxplot(bulk_mats$Factor1$AUPRCs[,cts],
                                                    height = unit(1.7, "cm"),
                                               axis_param=list(gp=gpar(fontsize = 11))))

F1_bulk_mapping <- ComplexHeatmap::Heatmap(bulk_mats$Factor1$ct_score[, cts], name = "ES",
                        right_annotation = ra,
                        top_annotation = ca,
                        cluster_columns = FALSE,
                        show_row_dend = FALSE,
                        show_row_names = FALSE,
                        border = TRUE,
                        col = color_fun,
                        column_names_gp = gpar(fontsize = 10))

# Factor 2 ----------------------------------------------------------------------

ra <- rowAnnotation(HF = bulk_mats$Factor2$sample_info[rownames(bulk_mats$Factor2$ct_score),
                                                       "HeartFailure"],
                    col = list(HF = c("yes" = "black", "no" = "darkgrey")),
                    gap = unit(2.5, "mm"),
                    border = TRUE)

ca <-  HeatmapAnnotation(AUROC =  anno_boxplot(bulk_mats$Factor2$AUPRCs[,cts],
                                               height = unit(1.7, "cm"),
                                               axis_param=list(gp=gpar(fontsize = 11))))


F2_bulk_mapping <- ComplexHeatmap::Heatmap(bulk_mats$Factor2$ct_score[, cts], name = "ES",
                                           right_annotation = ra,
                                           top_annotation = ca,
                                           cluster_columns = FALSE,
                                           show_row_dend = FALSE,
                                           show_row_names = FALSE,
                                           border = TRUE,
                                           col = color_fun,
                                           column_names_gp = gpar(fontsize = 10))




## S. scores

sillohuette_res <- map(set_names(F_analyze), function(f_name) {

 scores <- bulk_mats[[f_name]]$ct_score

 meta_info <- bulk_mats[[f_name]]$sample_info %>%
   dplyr::mutate(label_ix = .data[["HeartFailure"]] %>%
                   as.factor() %>%
                   as.integer())

 scores <- scores[rownames(meta_info),]

 patient_dists <- dist(scores,
                       method = "euclidean")


 si <- (silhouette(x = meta_info$label_ix, patient_dists)) %>%
   as.data.frame() %>%
   dplyr::select(cluster, sil_width) %>%
   left_join(unique(meta_info[,c("HeartFailure", "label_ix")]),
             by = c("cluster" = "label_ix")) %>%
   dplyr::mutate(factor_name = f_name)

})

# Testing differences between factors
batch_sw_plt <- bind_rows(sillohuette_res) %>%
  dplyr::mutate(HeartFailure = ifelse(HeartFailure == "no", "NF", "HF")) %>%
  ggplot(.,
         aes(y = HeartFailure, x = sil_width, color = factor_name)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  theme(plot.margin = unit(c(0.5,0,0,0.5), "cm")) + #
  scale_color_manual(values = c("black", "grey")) +
  xlab("silhouette \n width") +
  ylab("patient group")

# Export results

pdf("results/bulk_integration/reheat_mapF1.pdf", height = 3.1, width = 2.5)

plot(F1_bulk_mapping)

dev.off()


pdf("results/bulk_integration/reheat_mapF2.pdf", height = 3.1, width = 2.5)

plot(F2_bulk_mapping)

dev.off()


pdf("results/bulk_integration/reheat_sw.pdf",  height = 2, width = 3)

plot(batch_sw_plt)

dev.off()

# Paper plans

median(bulk_mats$Factor1$AUPRCs[,cts])

median(bulk_mats$Factor2$AUPRCs[,cts])

bind_rows(sillohuette_res) %>%
  group_by(factor_name) %>%
  summarize(median(sil_width))

# Make source data
# AUROC

F1_performance <- bulk_mats$Factor1$AUPRCs %>%
  as.data.frame() %>%
  rownames_to_column("bulk_study") %>%
  pivot_longer(-bulk_study, names_to = "cell_type", values_to = "Factor1_AUROC")

F2_performance <- bulk_mats$Factor2$AUPRCs %>%
  as.data.frame() %>%
  rownames_to_column("bulk_study") %>%
  pivot_longer(-bulk_study, names_to = "cell_type", values_to = "Factor2_AUROC")

left_join(F1_performance, F2_performance, by = c("bulk_study", "cell_type")) %>%
  write_csv("./Revision/figures/Figure5/Figure5A_up.csv")


#
F1_ES <- bulk_mats$Factor1$ct_score %>%
  as.data.frame() %>%
  rownames_to_column("bulk_sample") %>%
  pivot_longer(-bulk_sample,
               names_to = "cell_type",
               values_to = "Factor1_ES")

F2_ES <- bulk_mats$Factor2$ct_score %>%
  as.data.frame() %>%
  rownames_to_column("bulk_sample") %>%
  pivot_longer(-bulk_sample,
               names_to = "cell_type",
               values_to = "Factor2_ES")

left_join(F1_ES,F2_ES, by =c("bulk_sample", "cell_type")) %>%
  left_join(bulk_mats$Factor2$sample_info %>%
              as.data.frame() %>%
              rownames_to_column("bulk_sample")) %>%
write_csv("./Revision/figures/Figure5/Figure5A_down.csv")

