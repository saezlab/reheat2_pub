# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Transform ensenmbl IDs to gene names for the MOFA model
#' Kuppe and Reichart needed it

library(tidyverse)
library(biomaRt)

setwd("/Users/ricardoramirez/Dropbox/PostDoc/Research/ReHeaT2/")

input_df <- tibble(file = list.files("./data/pbulk/")) %>%
  dplyr::mutate(study = gsub("_pbulk[.]csv","",file)) %>%
  dplyr::mutate(file = paste0("./data/pbulk/",file)) %>%
  dplyr::filter(study %in% c("Reichart2022_DCM", "Kuppe2022_MI"))

# Importing pb data
pwalk(input_df, function(file, study,marker_csv) {
  
  pb_data <- read_csv(file,
                      show_col_types = FALSE)
  
  colnames(pb_data)[1] <- "sample_id"
  
  pb_data <- pb_data %>%
    column_to_rownames("sample_id") %>%
    as.matrix() %>%
    t()
  
  # Importing mart
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  genes <- rownames(pb_data)
  gene_dict <- getBM(filters= "ensembl_gene_id", 
                     attributes= c("ensembl_gene_id","hgnc_symbol"),
                     values=genes,mart= mart)
  
  # Exclude genes that have duplicated terms
  # Transforming
  pb_dat_red <- pb_data %>%
    as.data.frame() %>%
    rownames_to_column("ensembl_gene_id") %>%
    pivot_longer(-ensembl_gene_id) %>%
    left_join(gene_dict,  by = "ensembl_gene_id") %>%
    na.omit() %>%
    dplyr::select(-ensembl_gene_id) %>%
    group_by(name, hgnc_symbol) %>%
    summarise(sum_value = sum(value)) %>%
    pivot_wider(names_from = name, values_from = sum_value)
  
  # Saving is bugging
  pb_data <- pb_dat_red %>%
    dplyr::filter(hgnc_symbol != "") %>%
    column_to_rownames("hgnc_symbol") %>%
    as.matrix() %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("sample_id")
  
  write_csv(pb_data, file)
  
})




