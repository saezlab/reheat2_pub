library(tidyverse)

Fscores <- read_csv("./results/predict_project/LVAD/LVAD_factor_coords.csv")

pre_samples <- Fscores %>%
  dplyr::filter(heart_failure == "HF") %>%
  dplyr::mutate(patient_id = map_chr(sample_id, function(sid){
    vect <- strsplit(sid,"-") %>%
      unlist()
    return(paste0("p_",vect[2]))
  })) %>%
  dplyr::filter(biopsy == "pre")

Cdata <- read_csv("./data/val_studies/clinical/LVADRec_clinical.csv") %>%
  dplyr::rename("patient_id" = "LVAD ID") %>%
  dplyr::mutate(patient_id = map_chr(patient_id, ~paste0("p_", .x)))

# We select numeric ones and do Pearson correlations

Cdata_numeric <- Cdata %>%
  dplyr::select(patient_id, where(is.numeric)) %>%
  pivot_longer(-patient_id) %>%
  group_by(name) %>%
  nest()

Cdata_numeric_lm <- map2_df(Cdata_numeric$name, Cdata_numeric$data,
                            function(name, dat) {

  testdat <- dat %>%
    left_join(pre_samples, by = "patient_id")

  return(cor.test(testdat$Factor1, testdat$value) %>%
    broom::tidy() %>%
      dplyr::mutate(covar = name))

}) %>%
  dplyr::mutate(adj_p = p.adjust(p.value)) %>%
  dplyr::select(adj_p, covar) %>%
  dplyr::mutate(test = "correlation")

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

Cdata_categoric_lm <- map2_df(Cdata_cat$name, Cdata_cat$data,
                            function(name, dat) {

                              testdat <- dat %>%
                                left_join(pre_samples, by = "patient_id")

                              return(broom::tidy(aov(Factor1 ~ value, testdat)) %>%
                                dplyr::filter(term == "value") %>%
                                dplyr::mutate(covar = name))

                            }) %>%
  dplyr::mutate(adj_p = p.adjust(p.value)) %>%
  dplyr::select(adj_p, covar) %>%
  dplyr::mutate(test = "anova")

# Plot showing the results

bind_rows(Cdata_categoric_lm, Cdata_numeric_lm) %>%
  dplyr::filter(covar != "...6") %>%
  ggplot(aes(y = adj_p, x = covar)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  geom_bar(stat = "identity", color = "black", fill = "white") +
  geom_hline(yintercept = 0.05, color = "red") +
  facet_grid(. ~ test,scales = "free_x", space = "free_x") +
  ylab("adjusted_pvalue") +
  xlab("")

bind_rows(Cdata_categoric_lm, Cdata_numeric_lm) %>%
  dplyr::filter(covar != "...6") %>%
  ggplot(aes(x = covar, y = adj_p)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  geom_segment(aes(x = covar, xend = covar, y = 0, yend = adj_p),
               color = "grey50", linewidth = 1) +  # Line from 0 to value
  geom_point(size = 4, color = "black") +  # Points at the tip
  geom_hline(yintercept = 0.05, color = "red", linetype = "dashed") +  # Threshold line
  facet_grid(. ~ test, scales = "free_x", space = "free_x") +
  ylab("Adjusted p-value") +
  xlab("")
