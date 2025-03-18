# Copyright (c) [2024] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we associate extended covariates with
#' MCPs in the consensus analysis
library(MOFAcellulaR)
library(tidyverse)
library(lmerTest)

reheat_model <- MOFA2::load_model("./results/meta_mofacell/metamodel_mofa.hdf5")
comb_meta <- read_csv("./data/misc/combined_meta.csv")[,-1]

comb_meta_HF <- comb_meta %>%
  dplyr::filter(heart_failure == "HF") %>%
  dplyr::filter(!disease_code %in% c("ARVC","NCC")) %>%
  dplyr::select_at(c("study", "sample_id","sample_site",
                     "sample_accquisition", "genetic_HF", "disease_code",
                     "age", "sex"))

assocs_df <- map(set_names(c("Factor1", "Factor2")), function(Fact) {
  print(Fact)

  terms <- c("sample_site", "sample_accquisition",
             "genetic_HF", "age", "sex", "disease_code")

  map(set_names(terms), function(term) {

    print(term)


    if(term == "sample_accquisition") {

      data <- MOFAcellulaR::get_tidy_factors(model = reheat_model,
                                             metadata = comb_meta_HF,
                                             factor = Fact,
                                             group = T,
                                             sample_id_column = "sample_id") %>%
        dplyr::select_at(c("sample","study", "sample_accquisition", "value")) %>%
        dplyr::filter(sample_accquisition %in% c("Explanted", "LVAD"))

      model <- lmerTest::lmer(as.formula(paste0("value ~ ", term, " + (1|study)")),
                              data = data)

      return(broom::tidy(anova(model)))

      # model <- summary(model)
      #
      # res <- model$coefficients %>%
      #   as.data.frame() %>%
      #   rownames_to_column("coef") %>%
      #   dplyr::filter(coef != "(Intercept)") %>%
      #   dplyr:::mutate(test = term) %>%
      #   dplyr::rename("pval" = "Pr(>|t|)") %>%
      #   dplyr::select(test,pval)
      #return(res)

    } else {

      data <- MOFAcellulaR::get_tidy_factors(model = reheat_model,
                                             metadata = comb_meta_HF,
                                             factor = Fact,
                                             group = T,
                                             sample_id_column = "sample_id") %>%
        dplyr::select_at(c("sample","study", term, "value")) %>%
        na.omit()

      model <- lmerTest::lmer(as.formula(paste0("value ~ ", term, " + (1|study)")),
                              data = data)

      return(broom::tidy(anova(model)))

      #model <- summary(model)

      # res <- model$coefficients %>%
      #   as.data.frame() %>%
      #   rownames_to_column("coef") %>%
      #   dplyr::filter(coef != "(Intercept)") %>%
      #   dplyr:::mutate(test = term) %>%
      #   dplyr::rename("pval" = "Pr(>|t|)") %>%
      #   dplyr::select(test,pval)

      #return(res)

    }

  }) %>%
    enframe("term") %>%
    unnest()

}) %>%
  enframe("Factor") %>%
  unnest()

assocs_df <- assocs_df %>%
  group_by(Factor) %>%
dplyr::mutate(adj_p = p.adjust(p.value))

write_csv(assocs_df, "./results/meta_mofacell/HF_lmer_covars.csv")




