# HEADER --------------------------------------------
#
# Author:     Jan D. Lanzer
# Copyright   Copyright 2024 - Jan D. Lanzer
# Email:      lanzerjan@gmail.com
#
# Date:     2024-05-31
#
# Script Name:    ~/R-projects/reheat2_pilot/process_studies/process_sarkoidosis.R
#
# Script Description:
# process sarkoidosis study

library(Seurat)
library(tidyverse)
library(SingleCellExperiment)
library(scDblFinder)
library(cluster)

#source("../neonatal_cardiac_hyperplasia/utils.R")

x= Seurat::Read10X(data.dir= "data/sarkoidosis/GSM622406/")

## D7
GSM6222406="data/sarkoidosis/GSM622406/"
GSM6222407="data/sarkoidosis/GSM6222407_TWCM-11-173/"
GSM6222408="data/sarkoidosis/GSM6222408_TWCM-11-193/"
GSM6222409="data/sarkoidosis/GSM6222409_TWCM-574/"
GSM6222410="data/sarkoidosis/GSM6222410_TWCM-676/"
GSM6222411="data/sarkoidosis/GSM6222411_TWCM-707/"
GSM6222412="data/sarkoidosis/GSM6222412_TWCM-761/"

samples_l= list(GSM6222406, GSM6222407,GSM6222408,GSM6222409,GSM6222410,GSM6222411, GSM6222412)

names(samples_l)<-c("GSM6222406", "GSM6222407", "GSM6222408", "GSM6222409", "GSM6222410", "GSM6222411", "GSM6222412")
i= 1

for (i in 1:length(names(samples_l))){
  print(samples_l[i])
  assign(names(samples_l)[i], Read10X(unlist(samples_l[i])))
}



samples= list.files("data/raw/Raw_Data_SCS/")

sample_names= map(str_split(samples, "_"), function(x){x[[1]]})%>% unlist()

paths= paste0("data/raw/Raw_Data_SCS/", samples, "/outs/filtered_feature_bc_matrix")


# 
# # test vars:
# input_data= ct1
# sample_name= "ct1"

dis_score_m= read_csv("../scell_hfpef/data/prior_knowledge/coregene_df-FALSE-v3_MURINE.csv")

dis_score_feats <- list("disease_score" = dis_score_m[1:250,] %>% pull(x))

#translate gene cycle scoring genes:
gene_translate= readRDS("../scell_hfpef/data/prior_knowledge/gene_translate.rds")
cc.genes.m = list("s.genes"= gene_translate %>%
                    filter(Gene.name %in% cc.genes.updated.2019$s.genes)%>%
                    pull(MGI.symbol),
                  "g2m.genes"= gene_translate %>%
                    filter(Gene.name %in% cc.genes.updated.2019$g2m.genes)%>%
                    pull(MGI.symbol)
)

# main function for sample wise processing: ---------------------------------------------------

process_data= function(sample_name, input_data,
                       dis_filter= F){
  
  # Initialize the Seurat object
  sample_seurat <- Seurat::CreateSeuratObject(counts = input_data,
                                              project = sample_name,
                                              min.cells = 10,
                                              min.features = 250)
  rm(input_data)
  DefaultAssay(sample_seurat)
  sample_seurat@assays$RNA
  # Get mitochondrial genes
  sample_seurat[["percent.mt"]] <- PercentageFeatureSet(sample_seurat, pattern = "^MT-")
  # Get ribo genes
  sample_seurat[["percent.rb"]] <- PercentageFeatureSet(sample_seurat, pattern = "^RB[sl]")
  hist(sample_seurat$percent.mt, breaks= 100)
  hist(sample_seurat$percent.rb, breaks= 100)
  #sample_seurat$percent.rb[sample_seurat$percent.rb>0.1]
  
  print("Filtering cells by gene expression and mitochondrial gene expression")
  
  # Get a less stringent quantile check (since we are taking doublets anyway)
  # I will take the 1%
  nfeature_filter <- quantile(sample_seurat$nFeature_RNA,
                              1-0.01)
  
  # Here we make a collection of plots needed ------------------------------
  # N features and ncount
  
  filt_p1 <- sample_seurat@meta.data %>%
    ggplot(aes(x = nCount_RNA, y = nFeature_RNA)) +
    geom_point() +
    theme_classic() +
    geom_hline(yintercept = 300) +
    geom_hline(yintercept = nfeature_filter) +
    geom_vline(xintercept = 500) +
    ggtitle(paste0("ncells ", ncol(sample_seurat)))
  
  filt_p2 <- sample_seurat@meta.data %>%
    ggplot(aes(x = nCount_RNA, y = percent.mt)) +
    geom_point() +
    theme_classic() +
    geom_hline(yintercept = 25) +
    geom_vline(xintercept = 500) +
    ggtitle(paste0("ncells ", ncol(sample_seurat)))
  
  filt_p3 <- sample_seurat@meta.data %>%
    ggplot(aes(x = nCount_RNA, y = percent.rb)) +
    geom_point() +
    # theme_classic() +
    # geom_hline(yintercept = 300) +
    # geom_hline(yintercept = nfeature_filter) +
    # geom_vline(xintercept = 500) +
    ggtitle(paste0("ncells ", ncol(sample_seurat)))
  
  
  # Get mitochondrial genes
  sample_seurat <- subset(sample_seurat,
                          subset = nFeature_RNA > 300 &
                            nFeature_RNA < nfeature_filter &
                            percent.mt < 0.1 & ## single nuc stricter filter. 
                            percent.rb< 0.1 &
                            nCount_RNA > 500)
  
  print("Identifying doublets")
  
  # Identify doublets
  
  sce <- SingleCellExperiment(assays = list(counts = as.matrix(sample_seurat@assays$RNA$counts)))
  doublets <- scDblFinder(sce =sce)
  sample_seurat$doublet_score <- doublets$scDblFinder.score
  sample_seurat$doublet <- doublets$scDblFinder.class
  
  print("Getting QC info")
  
  
  # Process the data --------------------------------------------------------------------
  sample_seurat <- NormalizeData(sample_seurat,
                                 normalization.method = 'LogNormalize',
                                 scale.factor = 10000,
                                 verbose = FALSE)
  
  sample_seurat <- FindVariableFeatures(sample_seurat,
                                        selection.method = 'vst',
                                        nfeatures = 3000,
                                        verbose = FALSE)
  
  sample_seurat <- ScaleData(sample_seurat,
                             verbose = FALSE,
                             features = rownames(sample_seurat))
  
  sample_seurat <- RunPCA(sample_seurat,
                          verbose = FALSE)
  
  sample_seurat <- RunUMAP(sample_seurat, reduction = 'pca', dims = 1:30, verbose = FALSE)
  
  quality_plt <- FeaturePlot(sample_seurat, features = c("doublet_score",
                                                         "percent.mt",
                                                         "nCount_RNA",
                                                         "nFeature_RNA"))
  
  quality_plt_bis <- DimPlot(sample_seurat, group.by = "doublet")
  
  # Cell Cycle Score:
  sample_seurat= CellCycleScoring(sample_seurat,
                                  s.features = cc.genes.updated.2019$s.genes,
                                  g2m.features = cc.genes.updated.2019$g2m.genes)
  
  quality_plt_cellcycle <- FeaturePlot(sample_seurat, features = "G2M.Score")
  quality_plt_cellcycle2 <- FeaturePlot(sample_seurat, features = "S.Score")
  # Dissociation score
  # sample_seurat <- AddModuleScore(sample_seurat,
  #                                 assay = DefaultAssay(sample_seurat),
  #                                 features = toupper(dis_score_feats),
  #                                 name = "dissociation_s")
  # 
  # quality_plt_bis_2 <- FeaturePlot(sample_seurat, features = 'dissociation_s1')
  # 
  # if(dis_filter){
  #   dis_filter <- quantile(sample_seurat$dissociation_s1,
  #                          1-0.01)
  #   
  # }else{
  #   dis_filter= max(sample_seurat$dissociation_s1)
  #   
  # }
  
  # remove Rb-proteins
  keep_genes= rownames(sample_seurat)[!grepl(pattern = "Rp[sl]", rownames(sample_seurat))]
  sample_seurat <- subset(sample_seurat, features = keep_genes)
  
  # Filter and do it again --------------------------------------------------------------------
  print("Getting variable features, scaling and low dimension representations")
  
  sample_seurat <- subset(sample_seurat,
                          subset = doublet == "singlet")
  
  sample_seurat <- FindVariableFeatures(sample_seurat,
                                        selection.method = 'vst',
                                        nfeatures = 3000,
                                        verbose = FALSE)
  
  sample_seurat <- ScaleData(sample_seurat,
                             verbose = FALSE,
                             features = rownames(sample_seurat))
  
  sample_seurat <- RunPCA(sample_seurat,
                          verbose = FALSE)
  
  sample_seurat <- RunUMAP(sample_seurat,
                           reduction = 'pca',
                           dims = 1:30,
                           verbose = FALSE)
  
  print("Optimizing clustering")
  
  sample_seurat <- FindNeighbors(sample_seurat, reduction = "pca", dims = 1:30)
  
  seq_res <- seq(0.3, 1.5, 0.1)
  
  sample_seurat <- FindClusters(sample_seurat,
                                resolution = seq_res,
                                verbose = F)
  
  
  cell_dists <- dist(sample_seurat@reductions$pca@cell.embeddings,
                     method = "euclidean")
  
  cluster_info <- sample_seurat@meta.data[,grepl("RNA_snn_res",
                                                 colnames(sample_seurat@meta.data))] %>%
    dplyr::mutate_all(as.character) %>%
    dplyr::mutate_all(as.numeric)
  
  silhouette_res <- apply(cluster_info, 2, function(x){
    si <- silhouette(x, cell_dists)
    if(!is.na(si)) {
      mean(si[, 'sil_width'])
    } else {
      NA
    }
  })
  
  sample_seurat[["opt_clust"]] <- sample_seurat[[names(which.max(silhouette_res))]]
  
  spam_cols <- grepl(paste0(DefaultAssay(sample_seurat), "_snn_res"),
                     colnames(sample_seurat@meta.data)) |
    grepl("seurat_clusters",colnames(sample_seurat@meta.data))
  
  sample_seurat@meta.data <- sample_seurat@meta.data[,!spam_cols]
  
  #DimPlot(sample_seurat, group.by = "opt_clust")
  
  final_embedding <- DimPlot(sample_seurat, group.by = "opt_clust") +
    ggtitle(paste0("n cells ", ncol(sample_seurat)))
  
  print("Generating outputs")
  
  # Save data
  Idents(sample_seurat) = "opt_clust"
  
  saveRDS(sample_seurat, file = paste0("data/sarkoidosis/processed/", sample_name, ".rds"))
  
  # Plot QC files
  
  pdf(file = paste0("data/sarkoidosis/processed/", sample_name, ".pdf"), width = 8, height = 8)
  
  plot(cowplot::plot_grid(nrow = 1, ncol = 2, filt_p1, filt_p2))
  plot(quality_plt)
  plot(quality_plt_bis)
  #plot(quality_plt_bis_2)
  plot(final_embedding)
  plot(quality_plt_cellcycle2)
  dev.off()
  
}



# run for each sample -------------------------------------------------------------------------


for (name in names(samples_l)) {
  print(name)
  # Call the function using get to dynamically reference the variable
  process_data(sample_name = name, input_data = get(name))
}
