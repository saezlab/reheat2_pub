# HEADER --------------------------------------------
#
# Author:     Jan D. Lanzer
# Copyright   Copyright 2023 - Jan D. Lanzer
# Email:      lanzerjan@gmail.com
#
# Date:     2023-10-26
#
# Script Name:    
#
# Script Description:
# check enrihment of TYPE marker in STATE PB

library(tidyverse)
library(cowplot)
library(edgeR)
library(SummarizedExperiment)

pb= read.csv("output/fib_sub_analysis/all_state_pb_sub_reg.csv")%>% as_tibble()
col.dat= read.csv("output/fib_sub_analysis/all_state_coldata_sub_reg.csv")%>% as_tibble()## prepare coldata; add cell count:
col.dat= col.dat %>% 
  group_by(sample_state)%>% 
  summarise(cell_counts = n())%>%
  left_join(col.dat, by= "sample_state")%>% 
  select(sample_state, sample_id, cell_counts, leiden_0.5, cell_type, disease_code,heart_failure, batch)%>% 
  distinct()%>%
  dplyr::rename(cell_state= "leiden_0.5")

## prepare count data; 
# check for integer values 
# sum by gene
sum.g <- apply(pb[,colnames(pb)!= "X"], 2, sum)
# sum by sample
sum.s <- apply(pb[,colnames(pb)!= "X"], 1, sum)

# split the joint pb table into a list, separated by each study
df.split= pb %>% 
  dplyr::rename(sample_state= "X")%>% 
  left_join(col.dat%>% select(sample_state, batch))

dfs<-split.data.frame(df.split, df.split$batch)
dfs <- lapply(dfs, function(df) df[, -which(names(df) == "batch")])


#' Prepare pseudobulk data
#'
#' This function prepares pseudobulk RNA-Seq data by processing
#' a list of data frames containing expression counts and associated metadata.
#'
#' @param df_list A list of data frames, where each data frame contains expression counts.
#' @param col.dat metadata.
#' 
#' @return A list of processed pseudobulk profiles in long format together with coldata, one for each input data frame.
#'


prep_pb_data = function(df_list, col.dat){
  pb_dfs= map(names(df_list), function(x){
    print(x)
    
    counts= dfs[[x]] %>% 
      as.data.frame() %>% 
      column_to_rownames("sample_state")%>%
      as.matrix()%>% 
      t()
    
    coldata= col.dat %>% 
      filter(batch ==x) %>% 
      ungroup()%>%
      filter(sample_state %in% colnames(counts))%>%
      select(sample_state, everything()) %>%
      as.data.frame() %>% 
      column_to_rownames("sample_state")
    
    #check order
    print("dimensions of coldata and counts ")
    print(c(dim(coldata), dim(counts)))
    print("sample IDs are identical if all TRUE:")
    print(c(table(colnames(counts) %in% rownames(coldata)), table(rownames(coldata)%in% colnames(counts))
    )
    )
    
    #reorder counts: 
    counts= counts[,rownames(coldata)]
    
    # create summ exp object
    pb_dat <- 
      SummarizedExperiment(assays = list("counts" = counts), colData = DataFrame(coldata))
    colData(pb_dat)
    #filter for profiles and 
    ix <- which(colData(pb_dat)[,"cell_counts"] >= 20) # used to be 25
    
    useful_genes <- edgeR::filterByExpr(pb_dat,
                                        min.count = 10,
                                        min.prop = 0.1)
    
    pb_dat <- pb_dat[useful_genes, ix]
    colData(pb_dat)
    all_nf <- edgeR::calcNormFactors(pb_dat, method = "TMM")
    sfs <- all_nf$samples$lib.size * all_nf$samples$norm.factors
    pb <- sweep(assay(pb_dat, "counts"), MARGIN = 2, sfs, FUN = "/")
    assay(pb_dat, "logcounts") <- log1p(pb * 1000000)
    
    pb_dat
    
    dat <- assay(pb_dat, "logcounts") 
    
    rest_info <- colData(pb_dat) %>%
      as.data.frame() 
    
    long.dat= dat %>%
      as.data.frame() %>%
      tibble::rownames_to_column("feature") %>%
      pivot_longer(-feature, names_to = "sample_state", values_to = "value") %>%
      left_join(rest_info%>% rownames_to_column("sample_state"), by = "sample_state")
    
    long.dat
  })
  
  names(pb_dfs)= names(df_list)
  return(pb_dfs)
}

# run the processing and normalization
dfs2= prep_pb_data(dfs, col.dat = col.dat)

saveRDS(dfs2, "output/fib_sub_analysis/all_state_pb_processed_sub_reg.rds")

# process state marker to df ----------------------------------------------

## load state marker df
states= read.csv("output/fib_sub_analysis/cluster_markers_leiden_0.7_new_reg.csv")%>%
  as_tibble()
states
col_names = c("gene", "pval", "logfc", "cluster")
empty_df <- data.frame(matrix(ncol = length(empty_df <- data.frame(matrix(ncol = length(col_names), nrow = 0)))))
colnames(empty_df) <- col_names


for (i in 0:5){
  print(i)
  sub.df= states[,grepl(i, colnames(states))]%>% 
    mutate(cluster= paste0("Fib_", i))
  #df=rbind(df,states[,grepl(i, colnames(states))])
  print(sub.df)
  colnames(sub.df)= c("gene", "pval", "logfc", "cluster")
  empty_df= rbind(empty_df, sub.df)
}

states= empty_df %>% as_tibble()%>% drop_na()

write.csv(states, "output/fib_sub_analysis/cluster_markers_processed.csv")


