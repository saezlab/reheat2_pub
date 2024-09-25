## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2021-09-06
##
## Copyright (c) Jan D. Lanzer, 2021
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## Take single file ouput form run_sample_wise_preprocessing.R and integrate samples
## via harmony!

library(tidyverse)
library(Seurat)
library(harmony)
library(cluster)
library(clustree)
library(cowplot)
library(purrr)


path= "data/sarkoidosis/processed/"
slide_files= list.files(path)
slide_files<-slide_files[!grepl(pattern = ".pdf", slide_files)]
def_assay= "RNA"
# Because of incompatibility with objects objects should be appended manually
slide_files_path <- set_names(paste0(path,slide_files), gsub(pattern = "[.]rds",
                                                             replacement = "",
                                                             slide_files))

integrated_data <- map(slide_files_path, readRDS)

print("You managed to load everything")

# Calculate HVG per sample - Here we assume that batch and patient effects aren't as strong
# since cell-types and niches should be greater than the number of batches
hvg_list <- map(integrated_data, function(x) {

  DefaultAssay(x) <- def_assay

  x <- FindVariableFeatures(x, selection.method = "vst",
                            nfeatures = 3000)
  print(x)
  x@assays[[def_assay]]@meta.data$var.features

}) %>% unlist()

hvg_list <- table(hvg_list) %>%
  sort(decreasing = TRUE)

gene_selection_plt <- hvg_list %>% enframe() %>%
  group_by(value) %>%
  mutate(value = as.numeric(value)) %>%
  summarize(ngenes = length(name)) %>%
  ggplot(aes(x = value, y = ngenes)) +
  geom_bar(stat = "identity")

gene_selection <- hvg_list[1:3000] %>% names()

# Create merged object:
integrated_data <- purrr::reduce(integrated_data,
                          merge,
                          merge.data = TRUE)

print("You managed to merge everything")
print("Object size")
print(object.size(integrated_data))

DefaultAssay(integrated_data) <- def_assay

integrated_data= subset(integrated_data,
       subset =  percent.mt < 0.1 )

# Process it before integration -----------------------
integrated_data <- integrated_data %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(features = gene_selection,
         npcs = 30,
         verbose = FALSE)


original_pca_plt <- DimPlot(object = integrated_data,
                            reduction = "pca",
                            pt.size = .1,
                            group.by = "orig.ident")


# Add batch info
# must contain orig.ident
# batch_info <- read_csv(batch_file)
# batch_covars <- colnames(batch_info)
#
# tmp_meta <- integrated_data@meta.data %>%
#   left_join(batch_info, by = "orig.ident")
#
# integrated_data@meta.data <- bind_cols(integrated_data@meta.data, tmp_meta[, batch_covars[-1]])

# Integrate the data -----------------------
integrated_data <- RunHarmony(integrated_data,
                              group.by.vars = c("orig.ident"),
                              plot_convergence = TRUE,
                              assay.use = def_assay,
                              max.iter.harmony = 20)

corrected_pca_plt <- DimPlot(object = integrated_data,
                             reduction = "harmony",
                             pt.size = .1,
                             group.by = "orig.ident")

# Create the UMAP with new reduction -----------
integrated_data <- integrated_data %>%
  RunUMAP(reduction = "harmony", dims = 1:30,
          reduction.name = "umap_harmony") %>%
  RunUMAP(reduction = "pca", dims = 1:30,
          reduction.name = "umap_original")

integrated_data <- FindNeighbors(integrated_data,
                                 reduction = "harmony",
                                 dims = 1:30)

DefaultAssay(integrated_data) <- def_assay



optimize= T
if(optimize) {

  # Clustering and optimization -------------------------
  print("Optimizing clustering")

  seq_res <- seq(0.4, 1.5, 0.1)

  integrated_data <- FindClusters(object = integrated_data,
                                  resolution = seq_res,
                                  verbose = F)

  clustree_plt <- clustree(integrated_data,
                           prefix = paste0(DefaultAssay(integrated_data), "_snn_res."))

  # Optimize clustering ------------------------------------------------------
  cell_dists <- dist(integrated_data@reductions$harmony@cell.embeddings,
                     method = "euclidean")


  cluster_info <- integrated_data@meta.data[,grepl(paste0(DefaultAssay(integrated_data),"_snn_res"),
                                                   colnames(integrated_data@meta.data))] %>%
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

  integrated_data[["opt_clust_integrated"]] <- integrated_data[[names(which.max(silhouette_res))]]

  Idents(integrated_data) = "opt_clust_integrated"

  # Reduce meta-data -------------------------------------------------------------------------
  spam_cols <- grepl(paste0(DefaultAssay(integrated_data), "_snn_res"),
                     colnames(integrated_data@meta.data)) |
    grepl("seurat_clusters",colnames(integrated_data@meta.data))

  integrated_data@meta.data <- integrated_data@meta.data[,!spam_cols]

} else {

  print("Not Optimizing clustering")

  seq_res <- seq(0.4, 1.6, 0.2)

  integrated_data <- FindClusters(integrated_data,
                                  resolution = seq_res,
                                  verbose = F)

  clustree_plt <- clustree(integrated_data,
                           prefix = paste0(DefaultAssay(integrated_data),
                                           "_snn_res."))
  default_resolution= 0.3
  integrated_data <- FindClusters(integrated_data,
                                  resolution = default_resolution,
                                  verbose = F)

  integrated_data[["opt_clust_integrated"]] <- integrated_data[["seurat_clusters"]]

  spam_cols <- grepl(paste0(DefaultAssay(integrated_data), "_snn_res"),
                     colnames(integrated_data@meta.data)) |
    grepl("seurat_clusters",colnames(integrated_data@meta.data))

  integrated_data@meta.data <- integrated_data@meta.data[,!spam_cols]

}

DimPlot(integrated_data, group.by = "opt_clust_integrated")

DimPlot(integrated_data)
# add sample level meta data ----------------------------------------------

meta<- read_csv("~/Downloads/SraRunTable.txt")
meta<- meta %>% select('Library Name', 
                disease_state, tissue)%>%
  dplyr::rename(orig.ident = 'Library Name')
x= integrated_data@meta.data

x<-x %>% left_join(meta, by= "orig.ident")

integrated_data= AddMetaData(integrated_data, metadata = x[,])


# Save object ------------------------------------------------------
#save intermediate results (optional)
saveRDS(integrated_data, "data/additional_studies/sarkoidosis/integrated_v1.rds")
integrated_data= readRDS("data/additional_studies/sarkoidosis/integrated_v1.rds")

# get main cell annotation ------------------------------------------------------

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
  mutate(mor= 1)
Idents(integrated_data)

table(integrated_data$orig.ident, integrated_data$opt_clust_integrated)

mat <- as.matrix(integrated_data@assays$RNA$data)

library(decoupleR)
# Run ulm
acts <- run_ulm(mat=mat, net=net, .source='source', .target='target',
                .mor='mor', minsize = 5)
acts

#Extract ulm and store it in tfsulm in pbmc
integrated_data[['ctype_ulm']] <- acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = integrated_data) <- "ctype_ulm"

# Scale the data
integrated_data <- ScaleData(integrated_data)
integrated_data@assays$ctype_ulm@data <- integrated_data@assays$ctype_ulm@scale.data


Idents(integrated_data)= "opt_clust_integrated"

pdf("data/sarkoidosis/marker_enriched.pdf", 
    width = 10, 
    height= 6)
  DimPlot(integrated_data, label = T) | 
  DotPlot(integrated_data,assay= "ctype_ulm",  features = unique(net$source))
dev.off()

#integrated_data@meta.data$RNA_snn_res.0.3
c_marks= c("Vcam1", "Vwf", "Cd3e","Cd4", "Cd8a", "Pecam1", "Tie1", "S100a8", "S100a9")
c_marks= c("CCL21","CD3E", "CD4", "MZB1",  "KIT", "PLIN1", "LEPR", "NRXN1", "MZB1", "C1QC") 
DimPlot(integrated_data, label = T) | 
  DotPlot(integrated_data,assay= "RNA",  features = unique(toupper(c_marks)))+
                                                theme(axis.text.x = element_text(angle = 45, hjust= 1))

########

integrated_data@meta.data<- 
  integrated_data@meta.data%>%
  dplyr::rename(clust2= opt_clust_integrated)%>%
  mutate(cell_type = ifelse(clust2 %in% c(3, 5), "Endo", NA), 
         cell_type = ifelse(clust2 %in% c(0,11,4), "Fib", cell_type),
         cell_type = ifelse(clust2 %in% c(1), "CM", cell_type),
         cell_type = ifelse(clust2 %in% c(9), "Lymphoid", cell_type),
         cell_type = ifelse(clust2%in% c(2,16, 13, 15), "Myeloid", cell_type),
         cell_type = ifelse(clust2 %in% c(6), "PC", cell_type),
         cell_type = ifelse(clust2 %in% c(12), "vSMCs", cell_type)
         )
#unlabeled cells: 7, adipos; 8 LymphEC;
plot_grid(DimPlot(integrated_data, label = T, group.by = "cell_type"), 
          DimPlot(integrated_data, label = F, group.by = "disease_state"), 
          DimPlot(integrated_data, label = F, group.by = "clust2"),
          DotPlot(integrated_data,assay= "ctype_ulm",  features = unique(net$source))
)

saveRDS(integrated_data, "data/additional_studies/sarkoidosis/integrated_v2.rds")

# pseudobulking -----------------------------------------------------------

#we will do first pseudobulking per patient and cell type (for study projection)
#then we will save sample compositions and pseudobulk per patient (for deconvolution)

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


# per patient and cell type
pb.matrix<- get_pseudobulk(integrated_data,cluster_col = "cell_type", 
                           sample_col = "orig.ident", 
                           exclude_clusters = NA)


target_file <- integrated_data@meta.data %>%
  group_by(cell_type, orig.ident)%>% 
  summarise(cell_count= n())%>%
  left_join(integrated_data@meta.data%>% 
              distinct(cell_type,
                       orig.ident,  
                       disease_state))%>%
  mutate(pb_id = paste0(orig.ident, sep="_", cell_type))

saveRDS(list("pb"= pb.matrix, 
             "target"= target_file), 
        "data/Sarcoidosis_pb.rds")




