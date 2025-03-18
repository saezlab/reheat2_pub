# Copyright (c) [2025] [Ricardo O. Ramirez Flores]
# flores@ebi.ac.uk

#' Supervised analysis of interaction genes with sex and age

library(SingleCellExperiment)
library(scater)
library(edgeR)
library(tidyverse)
library(lmerTest)
library(MOFA2)

# For association with HF we don't need these disease codes from Kuppe et. al.
input_df <- tibble(file = list.files("./data/pbulk/")) %>%
  dplyr::mutate(study = gsub("_pbulk[.]csv","",file),
                coldata_file =  paste0("./data/coldata/", gsub("_pbulk[.]csv","_coldata.csv",file)),
                metadata_file =  paste0("./data/metadata_ext/", gsub("_pbulk[.]csv","_metadata.csv",file))) %>%
  dplyr::filter(study != "Reichart2022_DCM") %>%
  dplyr::mutate(file = paste0("./data/pbulk/",file)) %>%
  dplyr::select(file, study, coldata_file, metadata_file)

# First we need to get all markers together [Consensus]
all_mrks <- read_csv("./results/consensus_markers.csv", show_col_types = F) %>%
  dplyr::filter(adj_Fisher_p <= 0.0001, mean_LFC > 2) %>%
  dplyr::select(cell_type, gene) %>%
  unique() %>%
  group_by(cell_type) %>%
  nest() %>%
  dplyr::mutate(data = map(data, ~.x[[1]])) %>%
  deframe()

# All meta-data
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
    MOFAcellulaR::pb_dat2MOFA(pb_dat_list = .,
                              sample_column = "sample_id") %>%
    dplyr::mutate(group = study)

})

all_pbs <- bind_rows(all_pbs)

all_pbs <- all_pbs %>%
  dplyr::mutate(feature = map_chr(feature, function(x) {
    strsplit(x, "_")[[1]][2]
  }))

# Generating the meta-data for all
all_meta <-  bind_rows(map(input_df$metadata_file, read_csv))

# Extracting the info from the model
metamofacell <- MOFA2::load_model("./results/meta_mofacell/metamodel_mofa.hdf5")
cts <- MOFA2::views_names(metamofacell) %>% set_names()

human_factor <- MOFAcellulaR::get_geneweights(metamofacell, "Factor1") %>%
  #dplyr::filter(abs(value) > 0.1) %>%
  dplyr::mutate(direction = ifelse(value>0, "NF", "HF"))

ct_stats <- map(cts, function(ct){

  print(ct)

  # Filtering data to consider only cell-type and genes modelled

  genes2test <- human_factor %>%
    dplyr::filter(ctype == ct) %>%
    pull(feature) %>%
    unique()

  data2test <- all_pbs %>%
    dplyr::filter(view == ct,
                  feature %in% genes2test)

  # Running linear mixed models on age and sex, looking for interactions

  lmer_results <- data2test %>%
    group_by(feature) %>%
    nest() %>%
    dplyr::mutate(sup_tests = map(data, function(dat){

      comp_dat <- left_join(dat, all_meta, by = c("sample"="sample_id"))

      n_studies <- length(unique(comp_dat$group))

      if(n_studies > 1) {

        sex_model <- suppressMessages(lmerTest::lmer("value ~ heart_failure + sex + heart_failure*sex + (1 | group)",
                                    data = comp_dat,verbose = 0))

        sex_model_res <- summary(sex_model)

        sex_model_res <- sex_model_res$coefficients %>%
          as.data.frame() %>%
          rownames_to_column("term") %>%
          dplyr::mutate(test= "sex",
                        studies = n_studies) %>%
          dplyr::rename("pval" = "Pr(>|t|)",
                        "tval" = "t value") %>%
          dplyr::filter(grepl(":",term)) %>%
          dplyr::select(-df)

        age_model <- suppressMessages(lmerTest::lmer("value ~ heart_failure + age + heart_failure*age + (1 | group)",
                                    data = comp_dat,
                                    verbose = 0))

        age_model_res <- summary(age_model)

        age_model_res <- age_model_res$coefficients %>%
          as.data.frame() %>%
          rownames_to_column("term") %>%
          dplyr::mutate(test= "age",
                        studies = n_studies) %>%
          dplyr::rename("pval" = "Pr(>|t|)",
                        "tval" = "t value") %>%
          dplyr::filter(grepl(":",term)) %>%
          dplyr::select(-df)

        return(rbind(sex_model_res, age_model_res))

        # Regular linear model
      } else {

        sex_model <- lm("value ~ heart_failure + sex + heart_failure*sex",
                        data = comp_dat)

        sex_model_res <- summary(sex_model)

        sex_model_res <- sex_model_res$coefficients %>%
          as.data.frame() %>%
          rownames_to_column("term") %>%
          dplyr::mutate(test= "sex",
                        studies = n_studies) %>%
          dplyr::rename("pval" = "Pr(>|t|)",
                        "tval" = "t value") %>%
          dplyr::filter(grepl(":",term))

        age_model <- lm("value ~ heart_failure + age + heart_failure*age",
                        data = comp_dat)

        age_model_res <- summary(age_model)

        age_model_res <- age_model_res$coefficients %>%
          as.data.frame() %>%
          rownames_to_column("term") %>%
          dplyr::mutate(test= "age",
                        studies = n_studies) %>%
          dplyr::rename("pval" = "Pr(>|t|)",
                        "tval" = "t value") %>%
          dplyr::filter(grepl(":",term))

        return(rbind(sex_model_res, age_model_res))

      }

    }))

  stats <-  lmer_results %>%
    dplyr::select(feature, sup_tests) %>%
    unnest() %>%
    dplyr::group_by(test) %>%
    dplyr::mutate(padj = p.adjust(pval)) %>%
    arrange(padj)

})

all_res <- ct_stats %>%
  enframe() %>%
  unnest()

write_csv(all_res, "./results/meta_mofacell/interacting_sexage.csv")

all_res %>%
  dplyr::filter(padj <= 0.15)

# This is to plot interacting genes in case this is needed

plot_data <- all_pbs %>%
  dplyr::filter(view == "CM",
                feature == "PRKAG2") %>%
  left_join(all_meta, by = c("sample"="sample_id"))


ggplot(plot_data, aes(x = sex, fill = heart_failure, y = value)) +
  geom_boxplot() +
  facet_wrap(.~group, scales = "free_y")

ggplot(plot_data, aes(x = age, color = heart_failure, y = value)) +
  #geom_point() +
  geom_smooth(method = "lm", se = T) +
  facet_wrap(.~group, scales = "free_y")


