# HEADER --------------------------------------------
#
# Author:     Jan D. Lanzer
# Copyright   Copyright 2024 - Jan D. Lanzer
# Email:      lanzerjan@gmail.com
#
# Date:     2024-06-13
#
# Script Name:    
#
# Script Description: 
# perform differential gene expression analysis for LVAD and fetal expression
# to compare which part of the signature relates to Factor1 

library(limma)
library(edgeR)
library(tidyverse)



# AMRUTE ------------------------------------------------------------------
pb2 <- read_csv("data/val_studies/pbulk/Armute2023_LVAD_pbulk.csv")
amrute_target <- read_csv("data/val_studies/coldata/Armute2023_LVAD_coldata.csv")
# check if there are any 
amrute_target %>% filter(psbulk_n_cells<20)%>% print()

col<- amrute_target
dim(col)
df<-pb2%>% column_to_rownames('...1')%>%t()

col <- col%>% mutate(sample_cell = paste0(sample_id, "_", cell_type))

colnames(df)

df= df[,col$sample_cell]

celltypes <- unique(col$cell_type)

listed_M <- lapply(celltypes, function(x){
  col2 <- col %>% 
    filter(cell_type ==x)%>%
    filter(psbulk_n_cells >20)
  
  df2 <- df[, col2$sample_cell]
  
})

names(listed_M)<- celltypes

listed_M2 <- lapply(listed_M, function(x){
  
  meta <- col %>%filter(sample_cell %in% colnames(x))
  
  dge<- DGEList(counts=x, group= meta$disease)#, group=group)
  
  #detect and remove low expressed gene
  keep <- filterByExpr(dge, min.count	= 15, min.total.count= 20, min.prop = 0.75)
  table(keep)
  
  dge <- dge[keep,,keep.lib.sizes=FALSE]
  
  dge <- calcNormFactors(dge)
  v <- voom(dge, plot=TRUE)
  
})

names(listed_M2)<- celltypes


# DGE analysis w limma ----------------------------------------------------

col <- as_tibble(col)
M<-listed_M2$CM

fits<-lapply(listed_M2, function(M){
  
  target <- col %>% filter(sample_cell %in% colnames(M))
  
  group = (target$disease_code)
  pb= M$E[, target$sample_cell]
  
  design = target %>% select(condition, response, sample_cell, disease_code)%>%
    distinct()%>%
    # mutate(sex= as.numeric(as.factor(sex)), 
    #       disease_code= as.numeric(as.factor(disease_code)))%>%
    as.data.frame()%>% 
    column_to_rownames("sample_cell")
  
  mm= model.matrix(~ 0 + condition+ response,
                   design)
  
  fit = lmFit(pb, design = mm)
  
  cont.matrix <- makeContrasts(
    post_r_v_post_nr= conditionRpost- conditionNRpost, 
    pre_r_v_pre_nr= conditionRpre- conditionNRpre,
    pre_r_v_post_r = conditionRpre-conditionRpost, 
    
    levels = mm
  )
  fit2 = contrasts.fit(fit, cont.matrix)
  fit2 = eBayes(fit2)
  return(fit2)
})

De_res=  map(colnames(fits[[1]]$coefficients), function(coef){
  print(coef)
  map(names(fits), function(clust){
    print(clust)
    DE_results = as.data.frame(topTable(fits[[clust]],
                                        coef= coef,
                                        adjust.method = "BH",
                                        number = Inf)) %>%
      tibble::rownames_to_column("gene") %>%
      arrange(desc(abs(t))) %>%
      as_tibble()%>%
      mutate(sig= ifelse(adj.P.Val      <0.05, "*", ""),
             sig= ifelse(adj.P.Val      <0.01, "**", sig),
             sig= ifelse(adj.P.Val      <0.001, "***", sig),
             sig2 = ifelse(sig == "", "n", "y"),
             clust= clust, 
             coef= coef)
    
    p.pval.dist= DE_results %>% 
      ggplot(., aes(x= P.Value))+geom_histogram(bins = 100)
    
    p.volcano= DE_results%>%
      ggplot(., aes(x= logFC, y= -log10(P.Value), color= factor(sig)))+
      geom_point()
    
    tab= DE_results%>% group_by(sig2)%>% summarise(count= n())%>%
      mutate(clust= clust, coef= coef)
    
    p= cowplot::plot_grid(p.volcano+ggtitle(paste("cluster:", clust, "; contrast: ",coef)),
                          p.pval.dist)
    return(list("plot"=p, 
                "numb.DE"= tab, 
                "toptab"= DE_results)
    )
  })
})

De_res

toptables= map(De_res, function(x){
  #print(x)
  map(x, function(y){
    y$toptab
  })%>% do.call(rbind, .)
})%>% do.call(rbind, .)

toptables%>%
  saveRDS("data/additional_studies/DGE_celltypes_AMRUTE.rds")


# FETAL -------------------------------------------------------------------

pb <- readRDS("data/additional_studies/pb_data/pediatricDCM_pb_sample_celltype.rds")
pb$target
col <-pb$target %>% drop_na()%>% 
  rename(psbulk_n_cells= cell_count, 
         sample_cell= pb_id)

df <- pb$pb
df= df[,col$sample_cell]

celltypes <- unique(col$cell_type)

listed_M <- lapply(celltypes, function(x){
  col2 <- col %>% 
    filter(cell_type ==x)%>%
    filter(psbulk_n_cells >20)
  
  df2 <- df[, col2$sample_cell]
  
})

names(listed_M)<- celltypes

listed_M2 <- lapply(listed_M, function(x){
  
  meta <- col %>%filter(sample_cell %in% colnames(x))
  
  dge<- DGEList(counts=x, group= meta$disease)#, group=group)
  
  #detect and remove low expressed gene
  keep <- filterByExpr(dge, min.count	= 5, min.total.count= 20, min.prop = 0.75)
  table(keep)
  
  dge <- dge[keep,,keep.lib.sizes=FALSE]
  
  dge <- calcNormFactors(dge)
  v <- voom(dge, plot=TRUE)
  
})

names(listed_M2)<- celltypes


# DGE analysis w limma ----------------------------------------------------
# col <- as_tibble(col)
M <-listed_M2$vSMCs

fits<-lapply(listed_M2, function(M){
  
  target <- col %>% filter(sample_cell %in% colnames(M))
  
  group = (target$group)
  pb= M$E[, target$sample_cell]
  
  design = target %>% select(group, sample_cell)%>%
    distinct()%>%
    # mutate(sex= as.numeric(as.factor(sex)), 
    #       disease_code= as.numeric(as.factor(disease_code)))%>%
    as.data.frame()%>% 
    column_to_rownames("sample_cell")
  
  print(design)
  mm= model.matrix(~ 0 + group,
                   design)
  
  fit = lmFit(pb, design = mm)
  if(length(colnames(mm))==3){
    cont.matrix <- makeContrasts(
      DCM_v_CT= groupDCM-groupCT,
      FET_v_CT = groupFETAL-groupCT,
      FET_v_DCM =  groupFETAL-groupDCM,
      levels = mm
    )
    
  }else{
    colnames_mm <- colnames(mm)
    contrast_formula <- paste0(colnames_mm[2], " - ",  colnames_mm[1])
    print(contrast_formula)
    cont.matrix <- makeContrasts(
      FET_v_DCM = contrast_formula,
      levels = mm
    )
  }
  
  fit2 = contrasts.fit(fit, cont.matrix)
  fit2 = eBayes(fit2)
  return(fit2)
})

De_res=  map(colnames(fits[[1]]$coefficients), function(coef){
  print(coef)
  map(names(fits), function(clust){
    print(clust)
    if(!coef %in% colnames(fits[[clust]]$coefficients)){
      return(NULL)
    }
    DE_results = as.data.frame(topTable(fits[[clust]],
                                        coef= coef,
                                        adjust.method = "BH",
                                        number = Inf)) %>%
      tibble::rownames_to_column("gene") %>%
      arrange(desc(abs(t))) %>%
      as_tibble()%>%
      mutate(sig= ifelse(adj.P.Val      <0.05, "*", ""),
             sig= ifelse(adj.P.Val      <0.01, "**", sig),
             sig= ifelse(adj.P.Val      <0.001, "***", sig),
             sig2 = ifelse(sig == "", "n", "y"),
             clust= clust, 
             coef= coef)
    
    p.pval.dist= DE_results %>% 
      ggplot(., aes(x= P.Value))+geom_histogram(bins = 100)
    
    p.volcano= DE_results%>%
      ggplot(., aes(x= logFC, y= -log10(P.Value), color= factor(sig)))+
      geom_point()
    
    tab= DE_results%>% group_by(sig2)%>% summarise(count= n())%>%
      mutate(clust= clust, coef= coef)
    
    p= cowplot::plot_grid(p.volcano+ggtitle(paste("cluster:", clust, "; contrast: ",coef)),
                          p.pval.dist)
    return(list("plot"=p, 
                "numb.DE"= tab, 
                "toptab"= DE_results)
    )
  })
})

De_res

toptables= map(De_res, function(x){
  #print(x)
  map(x, function(y){
    y$toptab
  })%>% do.call(rbind, .)
})%>% do.call(rbind, .)


toptables%>%
  saveRDS("data/additional_studies/DGE_celltypes_FETAL.rds")


