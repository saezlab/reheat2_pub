library(tidyverse)
library(MOFAcellulaR)

# First get the clinical data processed

Cdata <- read_csv("./data/val_studies/clinical/LVADRec_clinical.csv") %>%
  dplyr::rename("patient_id" = "LVAD ID") %>%
  dplyr::mutate(patient_id = map_chr(patient_id, ~paste0("p_", .x)))

# We select numeric ones and do Pearson correlations
Cdata_numeric <- Cdata %>%
  dplyr::select(patient_id, where(is.numeric)) %>%
  pivot_longer(-patient_id) %>%
  group_by(name) %>%
  nest()

# We select categorical variables
Cdata_cat <- Cdata %>%
  dplyr::select(patient_id, where(is.character))

# We want categories with more than two levels
test_cats <- Cdata_cat %>%
  pivot_longer(-patient_id) %>%
  dplyr::select(-patient_id) %>%
  unique() %>%
  group_by(name) %>%
  summarize(n_cats = n()) %>%
  dplyr::filter(n_cats > 1) %>%
  pull(name)

# We want levels within used categories with at least 3 samples
test_cats <- Cdata_cat %>%
  pivot_longer(-patient_id) %>%
  dplyr::select(-patient_id) %>%
  dplyr::filter(name %in% test_cats) %>%
  group_by(name, value) %>%
  summarize(n_cats = n()) %>%
  group_by(name) %>%
  summarize(n_levels = sum(n_cats > 2)) %>%
  dplyr::filter(n_levels > 1) %>%
  pull(name)

# From those categories, we model
Cdata_cat <- Cdata_cat %>%
  dplyr::select_at(c("patient_id",test_cats)) %>%
  pivot_longer(-patient_id) %>%
  group_by(name) %>%
  nest()

# Then we generate Amrutes pbulk
# Data to be processed
meta_data <- read_csv("./data/val_studies/metadata/Armute2023_LVAD_metadata.csv",
                      show_col_types = FALSE)[,-1] %>%
  dplyr::mutate(patient_id = strsplit(sample_id,"-[A-Z]") %>%
                  map_chr(., ~.x[[1]]))

# Project data to multicellular space

pb_data <- read_csv("./data/val_studies/pbulk/Armute2023_LVAD_pbulk.csv",
                    show_col_types = FALSE)

colnames(pb_data)[1] <- "sample_id"

pb_data <- pb_data %>%
  column_to_rownames("sample_id") %>%
  as.matrix() %>%
  t()

# Importing coldata of the matrices - ignoring not annotated cell types
coldat <- read_csv("./data/val_studies/coldata/Armute2023_LVAD_coldata.csv",
                   show_col_types = FALSE)[,-1]  %>%
  column_to_rownames("colname") %>%
  dplyr::rename(ncells = "counts") %>%
  dplyr::filter(cell_type != "none")

pb_data <- pb_data[,rownames(coldat)]

# Processing the test data
# Defining cts
cts <- coldat$cell_type %>%
  unique() %>%
  set_names()

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

# Now congruent genes

useful_genes <- read_csv("./data/val_studies/clinical/top_enriched_features_mcp_supporting.csv") %>%
  dplyr::filter(relevant_feature == "yes",
                coef == "pre_r_v_pre_nr",
                Factor == "Factor1",
                adj.P.Val < 0.05) %>%
  dplyr::select(feature, ctype) %>%
  dplyr::mutate(feature = paste0(ctype, "_",feature)) %>%
  pull(feature)

# Processing pb so as to do associations

patients <- coldat %>%
  dplyr::filter(biopsy == "pre") %>%
  dplyr::select(sample_id, biopsy) %>%
  unique() %>%
  pull("sample_id")

MOFAcell_obj <- MOFAcell_obj %>%
  dplyr::filter(feature %in% useful_genes,
                sample %in% patients) %>%
  dplyr::mutate(patient_id = map_chr(sample, function(sid){
    vect <- strsplit(sid,"-") %>%
      unlist()
    return(paste0("p_",vect[2]))
  })) %>%
  dplyr::filter(patient_id %in% Cdata$patient_id) %>%
  dplyr::select(feature, value, patient_id) %>%
  dplyr::rename("expression" = "value")

# Now associations, per gene we do numerical and categorical

assocs_res <- MOFAcell_obj %>%
  group_by(feature) %>%
  nest() %>%
  dplyr::mutate(assoc_stats = map(data, function(gdat){

    Cdata_numeric_lm <- map2_df(Cdata_numeric$name, Cdata_numeric$data,
                                function(name, dat) {

                                  testdat <- dat %>%
                                    left_join(gdat, by = "patient_id")

                                  return(cor.test(testdat$expression, testdat$value) %>%
                                           broom::tidy() %>%
                                           dplyr::mutate(covar = name))

                                }) %>%
      dplyr::mutate(adj_p = p.adjust(p.value)) %>%
      dplyr::select(adj_p, p.value,covar) %>%
      dplyr::mutate(test = "correlation")

    Cdata_categoric_lm <- map2_df(Cdata_cat$name, Cdata_cat$data,
                                  function(name, dat) {

                                    testdat <- dat %>%
                                      left_join(gdat, by = "patient_id")

                                    return(broom::tidy(aov(expression ~ value, testdat)) %>%
                                             dplyr::filter(term == "value") %>%
                                             dplyr::mutate(covar = name))

                                  }) %>%
      dplyr::mutate(adj_p = p.adjust(p.value)) %>%
      dplyr::select(adj_p, p.value, covar) %>%
      dplyr::mutate(test = "anova")

    return(bind_rows(Cdata_numeric_lm, Cdata_categoric_lm))

  }))

assocs_res %>%
  dplyr::select(-data) %>%
  unnest() %>%
  dplyr::filter(adj_p < 0.05)

