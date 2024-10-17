# Copyright (c) [2024] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Code for multicellular information network from MCPs

library(lme4)
library(lmerTest)
library(MuMIn)
library(decoupleR)
library(MOFAcellulaR)

# For association with HF we don't need these disease codes from Kuppe et. al.
input_df <- tibble(file = list.files("./data/pbulk/")) %>%
  dplyr::mutate(study = gsub("_pbulk[.]csv","",file),
                coldata_file =  paste0("./data/coldata/",
                                       gsub("_pbulk[.]csv","_coldata.csv",file)),
                metadata_file =  paste0("./data/metadata_ext/",
                                        gsub("_pbulk[.]csv","_metadata.csv",file))) %>%
  dplyr::mutate(file = paste0("./data/pbulk/",file)) %>%
  dplyr::filter(study %in% c("Chaffin2022_DCM", "Koenig2022_DCM",
                             "LVReichart2022_DCM", "Simonson2023_ICM"))


# First we need to get all markers together [Consensus]
all_mrks <- read_csv("./results/consensus_markers.csv", show_col_types = F) %>%
  dplyr::filter(adj_Fisher_p <= 0.0001, mean_LFC > 2) %>%
  dplyr::select(cell_type, gene) %>%
  unique() %>%
  group_by(cell_type) %>%
  nest() %>%
  dplyr::mutate(data = map(data, ~.x[[1]])) %>%
  deframe()

# We need all pbulks together
# Processing will be handled for each study individually
file = "./data/pbulk/Chaffin2022_DCM_pbulk.csv"
study = "Chaffin2022_DCM"
coldata_file = "./data/coldata/Chaffin2022_DCM_coldata.csv"
metadata_file = "./data/metadata_ext/Chaffin2022_DCM_metadata.csv"

all_pbs <- pmap(input_df, function(file, study, coldata_file, metadata_file) {

  print(study)

  # Importing pb data
  pb_data <- read_csv(file,
                      show_col_types = FALSE)

  colnames(pb_data)[1] <- "sample_id"

  pb_data <- pb_data %>%
    column_to_rownames("sample_id") %>%
    as.matrix() %>%
    t()

  # Importing meta_data
  meta_data <- read_csv(metadata_file,
                        show_col_types = FALSE)

  # Importing coldata of the matrices - ignoring not annotated cell types
  coldat <- read_csv(coldata_file,
                     show_col_types = FALSE)[,-1]  %>%
    column_to_rownames("colname") %>%
    dplyr::rename(ncells = "counts") %>%
    dplyr::filter(cell_type != "none") %>%
    dplyr::filter(sample_id %in% meta_data$sample_id) # This ensures that filtering happens

  pb_data <- pb_data[,rownames(coldat)]

  # Defining cts
  cts <- coldat$cell_type %>%
    unique() %>%
    set_names()

  # Defining parameters
  n_samples <- nrow(meta_data)
  min_samples <- (n_samples * 0.4) %>% floor()

  # Creating MOFAcell object to make manipulations
  MOFAcell_obj <- MOFAcellulaR::create_init_exp(pb_data, coldat) %>%
    MOFAcellulaR::filt_profiles(pb_dat = .,
                                cts = cts,
                                ncells = 20,
                                counts_col = "ncells",
                                ct_col = "cell_type") %>%
    MOFAcellulaR::filt_views_bysamples(pb_dat_list = .,
                                       nsamples = min_samples) %>%
    MOFAcellulaR::filt_gex_byexpr(pb_dat_list = .,
                                  min.count = 20,
                                  min.prop = 0.4) %>%
    MOFAcellulaR::filt_views_bygenes(pb_dat_list = .,
                                     ngenes = 50) %>%
    MOFAcellulaR::filt_samples_bycov(pb_dat_list = ., # Filtering of low quality samples
                                     prop_coverage = 0.97) %>%
    MOFAcellulaR::tmm_trns(pb_dat_list = .,
                           scale_factor = 1000000) %>%
    MOFAcellulaR::filt_gex_bybckgrnd(pb_dat_list = .,
                                     prior_mrks = all_mrks) %>%
    MOFAcellulaR::filt_views_bygenes(pb_dat_list = .,
                                     ngenes = 50) %>%
    MOFAcellulaR::center_views()

  return(MOFAcell_obj)

})

names(all_pbs) <- input_df$study

# Positive and negative enrichment requires the definition of two signatures

metamofacell <- MOFA2::load_model("./results/meta_mofacell/metamodel_mofa.hdf5")

cts <- MOFA2::views_names(metamofacell) %>% set_names()

factor_signs <- map(set_names(c("Factor1","Factor2")), function(Fact) {

  human_factor <- MOFAcellulaR::get_geneweights(metamofacell, Fact) %>%
    dplyr::filter(abs(value) > 0.1) %>%
    dplyr::mutate(direction = ifelse(value>0, "NF", "HF"))

  # This enriches all signatures

  enrichment_scores <- map(all_pbs, function(MOFAcell_obj) {

    signature_df <- purrr::map(cts, function(ct) {

      dat <- SummarizedExperiment::assay(MOFAcell_obj[[ct]], "logcounts")
      colnames(dat) <- SummarizedExperiment::colData(MOFAcell_obj[[ct]])[, "sample_id"]

      ct_net <- human_factor %>%
        dplyr::filter(ctype == ct) %>%
        dplyr::mutate(value = abs(value))

      enriched_sign <- decoupleR::run_ulm(mat = dat,
                                          network = ct_net,
                                          .source = "direction",
                                          .target = "feature",
                                          .mor = "value") %>%
        dplyr::select(source, condition, score) %>%
        pivot_wider(names_from = source, values_from = score) %>%
        column_to_rownames("condition") %>%
        as.matrix() %>%
        base::scale() %>%
        as.data.frame() %>%
        rownames_to_column("sample_id") %>%
        pivot_longer(-sample_id, names_to = "dir_set", values_to = "scaled_scores")
    }) %>%
      enframe("ct") %>%
      unnest()

  }) %>%
    enframe(., "study") %>%
    unnest()

})


mcell_infonets <- map(factor_signs, function(sign_mat) {

  # Per condition build an information network

  all_nets <- map(set_names(c("NF", "HF")), function(net_test) {

    net_data <- sign_mat %>%
      dplyr::filter(dir_set == net_test) %>%
      dplyr::select(sample_id, study, ct, scaled_scores) %>%
      pivot_wider(values_from = scaled_scores, names_from = ct) %>%
      dplyr::select(-sample_id)

    # Now make a linear mixed model per cell_type

    net_estimates <- map(cts, function(target_ct) {

      pred_cts <- cts[cts != target_ct]
      pred_cts <- paste(pred_cts, collapse = " + ")

      model <- lmerTest::lmer(as.formula(paste0(target_ct," ~ ",
                                                pred_cts, " + (1 | study)")),
                              data = net_data)

      # R2 of fixed effects - marginal (R2m)
      FE_R2 <- r.squaredGLMM(model)[,"R2m"]

      # Coefficients
      model_res <- summary(model)
      # Information edges
      inf_summary <- model_res$coefficients %>%
        as.data.frame() %>%
        rownames_to_column("predictor") %>%
        dplyr::filter(predictor != "(Intercept)") %>%
        dplyr::select(predictor, Estimate) %>%
        dplyr::mutate(target = target_ct,
                      R2 = FE_R2)

      return(inf_summary)

    }) %>%
      bind_rows()

    return(net_estimates)

  })

  all_nets <- all_nets %>%
    enframe("condition") %>%
    unnest() %>%
    dplyr::mutate(cor_estimate = Estimate * R2)

})

factor_outs <- map(set_names(c("Factor1","Factor2")), function(Fact) {

  MIN_file <- paste0("./results/MCP_characs/MIN_", Fact, ".csv")
  MIN_plot <- paste0("./results/MCP_characs/MIN_", Fact, ".pdf")
  DIF_plot <- paste0("./results/MCP_characs/MIN_", Fact,"_dif", ".pdf")

  all_nets <- mcell_infonets[[Fact]]

  write_csv(all_nets, MIN_file)
  # What's the actual relationships?

  MCIN <- all_nets %>%
    ggplot(aes(x = predictor, y = target, fill = cor_estimate)) +
    geom_tile() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    coord_equal() +
    scale_fill_gradient2(midpoint = 0) +
    facet_wrap(.~condition)

  pdf(MIN_plot, height = 4, width = 8)
  plot(MCIN)
  dev.off()

  # Do cells change their coordination roles on HF?
  HFvNF_diff <- all_nets %>%
    group_by(predictor) %>%
    nest() %>%
    dplyr::mutate(coord_comp = map(data, function(dat) {

      pairwise.t.test(dat$cor_estimate,
                      dat$condition,
                      p.adjust.method = "BH") %>%
        broom::tidy()

    })) %>%
    dplyr::select(predictor, coord_comp) %>%
    unnest() %>%
    ungroup() %>%
    dplyr::arrange((p.value))

  bic_coordinator <- all_nets %>%
    group_by(condition) %>%
    nest() %>%
    dplyr::mutate(coord_comp = map(data, function(dat) {

      pairwise.t.test(dat$cor_estimate,
                      dat$predictor,
                      p.adjust.method = "BH") %>%
        broom::tidy()

    })) %>%
    dplyr::select(condition, coord_comp) %>%
    unnest() %>%
    ungroup() %>%
    dplyr::arrange((p.value))

  difrntl_net <- all_nets %>%
    dplyr::select(condition, target, predictor, cor_estimate) %>%
    pivot_wider(names_from = condition, values_from = cor_estimate) %>%
    dplyr::mutate(diff = HF - NF)

  MIN_Dif <- ggplot(difrntl_net, aes(x = predictor, y = target, fill = diff)) +
    geom_tile() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    coord_equal() +
    scale_fill_gradient2(midpoint = 0)

  pdf(DIF_plot, height = 4, width = 8)
  plot(MIN_Dif)
  dev.off()

  list("HFvNF_diff" = HFvNF_diff,
       "bic_coordinator" = bic_coordinator)

})

# Compare coordinator roles across networks

difrntl_net_F <- mcell_infonets %>%
  enframe("Factor") %>%
  unnest() %>%
  dplyr::select(-c(Estimate, R2)) %>%
  group_by(condition) %>%
  nest() %>%
  dplyr::mutate(data = map(data, function(dat) {

    pivot_wider(dat, names_from = Factor, values_from = cor_estimate)

  })) %>%
  unnest() %>%
  dplyr::mutate(diff = Factor1 - Factor2)


MCP_comp <- difrntl_net_F %>%
  dplyr::filter(predictor != target) %>%
  ggplot(aes(x = predictor, y = target, fill = diff)) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  coord_equal() +
  scale_fill_gradient2(midpoint = 0,na.value = "grey") +
  facet_wrap(.~condition) +
  ggtitle("Factor1 - Factor2")


