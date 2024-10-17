
library(tidyverse)
library(Seurat)
library(harmony)
library(cluster)
library(clustree)
library(cowplot)
library(purrr)


path= "data/pediatric_DCM/custom_processed/"
slide_files= list.files(path)
slide_files <- slide_files[grepl(".rds", slide_files)]
def_assay= "RNA"

# Because of incompatibility with objects objects should be appended manually
slide_files_path <- set_names(paste0(path,slide_files), gsub(pattern = "[.]rds",
                                                             replacement = "",
                                                             slide_files))

integrated_data <- map(slide_files_path, readRDS)

# Calculate HVG per sample - Here we assume that batch and patient effects aren't as strong
# since cell-types and niches should be greater than the number of batches
hvg_list <- map(integrated_data, function(x) {
  
  DefaultAssay(x) <- def_assay
  
  x <- FindVariableFeatures(x, selection.method = "vst",
                            nfeatures = 3000)
  VariableFeatures(x)
  #x@assays[[def_assay]]$var.features
  
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

# integrated_data= subset(integrated_data,
#                         subset =  percent.mt < 0.1 )

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

DimPlot(integrated_data, group.by = "orig.ident")

integrated_data <- FindNeighbors(integrated_data,
                                 reduction = "harmony",
                                 dims = 1:30)

DefaultAssay(integrated_data) <- def_assay

#save intermediate results (optional)
saveRDS(integrated_data, "data/pediatric_DCM/integrated_v1.rds")
