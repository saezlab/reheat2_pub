# HEADER --------------------------------------------
#
# Author:     Jan D. Lanzer
# Copyright   Copyright 2024 - Jan D. Lanzer
# Email:      lanzerjan@gmail.com
#
# Date:     2024-05-31
#
#
# Script Description:
# process hcm data from dimmler

library(Seurat)
library(tidyverse)

source("process_studies/utils_process.R")
seu<- readRDS("data/additional_studies/seurat_object_hca_as_harmonized_AS_SP_nuc_refined_cells.rds")

seu$celltype%>% unique()

# filter smaples to get only HCM samples, map cell types to our onotlogy then pb

#create the mapping manually: 
#from
seu@meta.data%>% distinct(cell_type)
#to 
cell_dic$real

#
seu@meta.data <- seu@meta.data%>% 
  mutate(cell_type_translated = ifelse(cell_type %in% c("EC"), "Endo", NA), 
         cell_type_translated = ifelse(cell_type %in% c("CM"), "CM", cell_type_translated),
         cell_type_translated = ifelse(cell_type %in% c("PC"), "PC", cell_type_translated),
         cell_type_translated = ifelse(cell_type %in% c("SMC"), "vSMCs", cell_type_translated),
         cell_type_translated = ifelse(cell_type %in% c("Lympoid"), "Lymphoid", cell_type_translated),
         cell_type_translated = ifelse(cell_type %in% c("Myeloid"), "Myeloid", cell_type_translated),
         cell_type_translated = ifelse(cell_type %in% c("FB"), "Fib", cell_type_translated)
  )
seu= subset(seu, subset = (condition== "AS" & 
                             cell_type_translated %in% cell_dic$real)
            )
seu$cell_type_translated
DimPlot(seu, group.by = "cell_type_translated")
# pseudobulk --------------------------------------------------------------

pb.matrix<- get_pseudobulk(seu,cluster_col = "cell_type_translated", 
                           sample_col = "orig.ident", 
                           exclude_clusters = NA)

target_file <- seu@meta.data %>%
  group_by(cell_type_translated, orig.ident)%>% 
  summarise(cell_count= n())%>%
  left_join(seu@meta.data%>% 
              distinct(cell_type_translated,
                       orig.ident,  
                       condition,
                       treatment))%>%
  mutate(pb_id = paste0(orig.ident, sep="_", cell_type_translated))

saveRDS(list("pb"= pb.matrix, 
             "target"= target_file), 
          "data/additional_studies/HCM_dimmeler_pb.rds")



