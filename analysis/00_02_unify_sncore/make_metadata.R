# Copyright (c) [2024] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Generate files of metadata for extended models

centralized_meta <- read_csv("./data/misc/combined_meta.csv")[,-1] %>%
  dplyr::group_by(study) %>%
  nest()

study_dictionary <- tibble(study = c("Chaffin2022","Koenig2022", "Simonson2023", "Reichart2022"),
                           study_id = c("Chaffin2022_DCM","Koenig2022_DCM", "Simonson2023_ICM", "Reichart2022_DCM"))

centralized_meta <- left_join(centralized_meta, study_dictionary)

walk2(centralized_meta$study_id, centralized_meta$data, function(study_label, dat) {

  fname <- paste0("./data/metadata_ext/", study_label, "_metadata.csv")

  write_csv(dat, fname)

})




