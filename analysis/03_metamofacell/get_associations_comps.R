# Copyright (c) [2024] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we associate each covariate with
#' MCPs in the consensus analysis

library(MOFAcellulaR)
library(compositions)
library(tidyverse)
library(cowplot)
library(lmerTest)

assoc_csv <- "./results/MCP_characs/assocs_comps.csv"
assoc_plt <- "./results/MCP_characs/assocs_comps.pdf"

sample_compositions <- read_csv("./results/compositions/sample_comps.csv") %>%
  column_to_rownames("sample_id") %>%
  as.matrix() %>%
  compositions::clr() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id")

all_meta <- read_csv("./results/meta_mofacell/metamodel_meta.csv")[,c("study", "sample_id", "heart_failure")] %>%
  left_join(sample_compositions)

reheat_model <- MOFA2::load_model("./results/meta_mofacell/metamodel_mofa.hdf5")

# Run associations with MOFAcell model loadings
cts <- MOFA2::views_names(reheat_model) %>%
  set_names()

# Associations with linear mixed models
assocs_df <- map(set_names(c("Factor1", "Factor2")), function(Fact) {

  all_dat <- MOFAcellulaR::get_tidy_factors(model = reheat_model,
                                            metadata = all_meta,
                                            factor = Fact,
                                            group = T,
                                            sample_id_column = "sample_id") %>%
    dplyr::filter(heart_failure == "HF")

  map(cts, function(ct) {

    model <-  all_dat %>%
      lmerTest::lmer(as.formula(paste0("value ~ ", ct, " + (1|study)")), data = .)

    # R2 of fixed effects - marginal (R2m)

    FE_R2 <- r.squaredGLMM(model)[,"R2m"]

    # R2 of random effects

    model_res <- summary(model)

    RE_R2 <- model_res$varcor %>%
      as.data.frame() %>%
      dplyr::select(grp, vcov) %>%
      pivot_wider(names_from = grp, values_from = vcov) %>%
      dplyr::mutate(perc_studyvar = study/(study + Residual)) %>% # change study for
      pull(perc_studyvar)

    res_df = model_res$coefficients %>%
      as.data.frame() %>%
      rownames_to_column("ct") %>%
      dplyr::filter(ct != "(Intercept)") %>%
      dplyr::mutate("FE_R2" = FE_R2,
                    "RE_R2" = RE_R2)

    return(res_df)

  }) %>%
    enframe() %>%
    dplyr::select(-name) %>%
    unnest()

}) %>%
  enframe("Factor") %>%
  unnest() %>%
  group_by(Factor) %>%
  dplyr::mutate(adj_p = p.adjust(`Pr(>|t|)`))

write_csv(assocs_df, assoc_csv)

## Plots

MOFAcellulaR::get_tidy_factors(model = reheat_model,
                               metadata = all_meta,
                               factor = "all",
                               group = T,
                               sample_id_column = "sample_id") %>%
  dplyr::filter(heart_failure == "HF",
                Factor %in% c("Factor1","Factor2")) %>%
  ggplot(aes(x = value, y = CM, color = study)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  facet_wrap(.~Factor)

MOFAcellulaR::get_tidy_factors(model = reheat_model,
                               metadata = all_meta,
                               factor = "all",
                               group = T,
                               sample_id_column = "sample_id") %>%
  dplyr::filter(heart_failure == "HF",
                Factor %in% c("Factor1","Factor2")) %>%
  ggplot(aes(x = value, y = Myeloid, color = study)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  facet_wrap(.~Factor)


comp_assocs_plt <- assocs_df %>%
  dplyr::mutate(signif = ifelse(adj_p<= 0.1, "*", "")) %>%
  ggplot(aes(y = Factor, x = ct, fill = Estimate, label = signif)) +
  geom_tile() +
  geom_text(color = "white") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  coord_equal() +
  xlab("") +
  ylab("") +
  scale_fill_gradient2(midpoint = 0)

pdf(assoc_plt, height = 2.5, width = 2.5)
plot(comp_assocs_plt)
dev.off()





