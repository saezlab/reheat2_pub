# Copyright (c) [2023] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we calculate markers of cells
#' using edgeR and pseudobulk profiles of all samples

library(SingleCellExperiment)
library(scater)
library(edgeR)
library(tidyverse)

setwd("/Users/ricardoramirez/Dropbox/PostDoc/Research/ReHeaT2/")

# Make input parameter dataframe
input_df <- tibble(file = list.files("./data/pbulk/")) %>%
  dplyr::mutate(study = gsub("_pbulk[.]csv","",file),
                coldata_file =  paste0("./data/coldata/", gsub("_pbulk[.]csv","_coldata.csv",file)),
                metadata_file =  paste0("./data/metadata_ext/", gsub("_pbulk[.]csv","_metadata.csv",file)),
                marker_csv = paste0("./data/mrkrs/",
                                    gsub("_pbulk[.]csv","_mrkrs.csv",file))) %>%
  dplyr::mutate(file = paste0("./data/pbulk/",file)) %>%
  dplyr::select(study, file, coldata_file, marker_csv)


pwalk(input_df, function(file, study, coldata_file, marker_csv) {

  print(study)

  # Importing pb data
  pb_data <- read_csv(file,
                      show_col_types = FALSE)

  colnames(pb_data)[1] <- "sample_id"

  pb_data <- pb_data %>%
    column_to_rownames("sample_id") %>%
    as.matrix() %>%
    t()

  is_whole_var <- ((colSums(pb_data) %>% sum()) %% 1) == 0

  if(is_whole_var) {

    # Importing coldata of the matrices - ignoring not annotated cell types
    coldat <- read_csv(coldata_file,
                       show_col_types = FALSE)[,-1]  %>%
      column_to_rownames("colname") %>%
      dplyr::rename(ncells = "counts") %>%
      dplyr::filter(cell_type != "none")

    pb_data <- pb_data[,rownames(coldat)]

    # Defining cts
    cts <- coldat$cell_type %>%
      unique() %>%
      set_names()

    de_res <- map(cts, function(ct) {
      print(ct)

      ct_meta_data <- coldat %>%
        mutate(test_column = ifelse(cell_type == ct, ct, "rest"))

      dat <- DGEList(pb_data, samples = DataFrame(ct_meta_data))

      keep <- filterByExpr(dat, group = ct_meta_data$test_column)

      dat <- dat[keep,]

      dat <- calcNormFactors(dat)

      design <- model.matrix(~factor(test_column,
                                     levels = c("rest",ct)), dat$samples)

      colnames(design) <- c("int", ct)

      dat <- estimateDisp(dat, design)

      fit <- glmQLFit(dat, design, robust=TRUE)

      res <- glmQLFTest(fit, coef=ncol(design))

      de_res <- topTags(res, n = Inf) %>%
        as.data.frame() %>%
        rownames_to_column("gene")

      return(de_res)

    })

    de_res <- de_res %>%
      enframe() %>%
      unnest()

    de_res %>%
      dplyr::filter(logFC > 0) %>%
      arrange(name, FDR, - logFC) %>%
      write_csv(file = marker_csv)

  } else {
    print("Your study doesn't have counts, check your pseudobulk")
  }

})


