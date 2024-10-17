# HEADER --------------------------------------------
#
# Author:     Jan D. Lanzer
# Copyright   Copyright 2024 - Jan D. Lanzer
# Email:      lanzerjan@gmail.com
#
# Date:     2024-05-31
#
# Script Name:    ~/R-projects/reheat2_pilot/process_studies/provess_Hill_CHD.R
#
# Script Description:
# process HCa

directory = "./"
source(paste0(directory, "process_studies/utils_process.R"))
# sc ref: 
seu = readRDS(paste0(directory, "data/HCA_seu_filtered_nuc.rds"))

print(colnames(seu@meta.data))
print(unique(seu@meta.data$cell_type))
print(wanted_cells)

cell_dic<- read_csv("output/celltype_dictionary_HCA_ReHeaT2.csv")

seu@meta.data$cell_type_translated<- cell_dic[match(seu@meta.data$cell_type, cell_dic$pred),"real"]%>% pull(real)

# make sure to only use samples of interest in reference
#seu = Seurat::subset(seu, cell_type %in% wanted_cells)
seu= subset(seu, subset = (region %in% c("LV", "AX)") & 
                             cell_type_translated %in% cell_dic$real &
                             source == "Nuclei" & 
                             Used == "Yes" &
                             nCount_RNA >200 &  
                             nFeature_RNA > 200)
)

saveRDS(object = seu, "data/HCA_seu_filtered_nuc_translated.rds")

gex.count<- seu@assays$RNA$counts
print(unique(seu@meta.data$cell_type_translated))
print("input loaded")


# pbking: -----------------------------------------------------------------


pb.matrix<- get_pseudobulk(seu,cluster_col = "cell_type_translated", 
                           sample_col = "donor", 
                           exclude_clusters = NA)

target_file <- seu@meta.data %>%
  group_by(cell_type_translated, donor)%>% 
  summarise(cell_count= n())%>%
  left_join(seu@meta.data%>% 
              distinct(cell_type_translated,
                       donor,  
                       gender, 
                       age_group, 
                       region))%>%
  mutate(pb_id = paste0(donor, sep="_", cell_type_translated))

saveRDS(list("pb"= pb.matrix, 
             "target"= target_file), 
        "data/HCA_pb.rds")


