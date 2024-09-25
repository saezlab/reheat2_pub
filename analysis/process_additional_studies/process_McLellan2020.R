# HEADER --------------------------------------------
#
# Author:     Jan D. Lanzer
# Copyright   Copyright 2024 - Jan D. Lanzer
# Email:      lanzerjan@gmail.com
#
# Date:     2024-05-29
#
# Script Name:    
#
# Script Description:
# process AngII sc


library(Seurat)
library(tidyverse)

seu <- readRDS("~/R-projects/sc-exploration/data/circ_MI/integrated_seurat_circMI_processed.rds")
seu

DimPlot(seu)
seu$integrated_snn_res.0.3
Idents(seu)= "integrated_snn_res.0.2"
DimPlot(seu, group.by = "integrated_snn_res.0.2")

cmark=read.csv("data/consensus_markers.csv")%>% as_tibble()
#prepare the ctype marker, select n top genes
n_top_genes= 100

net= cmark%>%
  group_by(cell_type)%>% 
  filter(mean_LFC>0)%>%
  arrange(adj_Fisher_p, desc(mean_LFC))%>%
  slice_head(n= n_top_genes)%>%
  dplyr::rename(target= "gene", 
                source= "cell_type", 
  )%>%
  mutate(mor= 1)%>%
  mutate(target = str_to_title(target))
Idents(seu)

table(seu$orig.ident, seu$integrated_snn_res.0.2)

mat <- as.matrix(seu@assays$RNA@data)

library(decoupleR)
# Run ulm
acts <- run_ulm(mat=mat, net=net, .source='source', .target='target',
                .mor='mor', minsize = 5)
acts

#Extract ulm and store it in tfsulm in pbmc
seu[['ctype_ulm']] <- acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = seu) <- "ctype_ulm"

# Scale the data
seu <- ScaleData(seu)
seu@assays$ctype_ulm@data <- seu@assays$ctype_ulm@scale.data


Idents(seu)= "integrated_snn_res.0.3"

DimPlot(seu, label = T) | 
  DotPlot(seu,assay= "ctype_ulm",  features = unique(net$source))

#integrated_data@meta.data$RNA_snn_res.0.3
c_marks= c("Vcam1", "Vwf", "Cd3e","Cd4", "Cd8a", "Pecam1", "Tie1", "S100a8", "S100a9")
DimPlot(seu, label = T) | DotPlot(seu,assay= "RNA",  features = c_marks)


seu@meta.data<- 
  seu@meta.data%>%
  #rename(clust2= integrated_snn_res.0.3)%>%
  mutate(cell_type = ifelse(clust2 %in% c(2,7,4), "Endo", NA), 
         cell_type = ifelse(clust2 %in% c(0,3,6,17,10), "Fib", cell_type),
         cell_type = ifelse(clust2 %in% c(5), "CM", cell_type),
         cell_type = ifelse(clust2 == 12, "Lymphoid", cell_type),
         cell_type = ifelse(clust2 %in% c(1,9), "Myeloid", cell_type),
         cell_type = ifelse(clust2 == 8, "PC", cell_type),
         cell_type = ifelse(clust2 == 13, "vSMCs", cell_type))
library(cowplot)
plot_grid(DimPlot(seu, label = T, group.by = "cell_type"), 
          DimPlot(seu, label = F, group.by = "group"), 
          DimPlot(seu, label = F, group.by = "clust2"),
          DotPlot(seu,assay= "ctype_ulm",  features = unique(net$source))
)

cowplot::plot_grid(
  VlnPlot(seu, "nCount_RNA", group.by = "clust2"), 
  # VlnPlot(integrated_data, "nFeature_RNA", group.by = "cell_type")
)



# pseudobulking -----------------------------------------------------------


get_pseudobulk<- function(seu, 
                          cluster_col= "cell_type", 
                          sample_col= "orig.ident", 
                          exclude_samples= c(),
                          exclude_clusters= c()
                          
){
  
  ctypes<- unique(seu@meta.data[[cluster_col]])
  samps <- unique(seu@meta.data[[sample_col]])
  
  ctypes <- ctypes[!ctypes %in% exclude_clusters]
  samps <- samps[!samps %in% exclude_samples]
  
  x= map(samps, function(x){
    
    pb <-  sapply(ctypes, function(y){
      #x= samps[1]
      #y = ctypes[1]
      cells <- rownames(seu@meta.data %>% filter(.data[[sample_col]] == x, 
                                                 .data[[cluster_col]] == y)
      )#%>%
      length(cells)
      sample_1 <- seu@assays$RNA$counts[, cells ]
      dim(sample_1)
      print(c(x, y))
      print(dim(sample_1))
      if(is.null(dim(sample_1))){return(rep(0, dim(seu@assays$RNA$counts)[1]))}
      pb <- rowSums(sample_1)
      #rownames(pb)<- rownames(mtx.df)
      return(pb)
      
      
      
    })
    
    #rownames(pb)<- rownames(mtx.df)
    colnames(pb)<- paste0(x, sep ="_" ,ctypes)
    return(pb)
  })
  
  pb.matrix <- do.call(cbind, x)
}


pb.matrix<- get_pseudobulk(seu,cluster_col = "cell_type", 
                           sample_col = "orig.ident", 
                           exclude_clusters = NA)


target_file <- seu@meta.data %>%
  group_by(cell_type, orig.ident)%>% 
  summarise(cell_count= n())%>%
  left_join(seu@meta.data%>% 
              distinct(cell_type,
                       orig.ident,  
                       group))%>%
  mutate(pb_id = paste0(orig.ident, sep="_", cell_type))

saveRDS(list("pb"= pb.matrix, 
             "target"= target_file), 
        "data/AngII_pb_sample_celltype.rds")

