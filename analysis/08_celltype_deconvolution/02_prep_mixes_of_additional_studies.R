# HEADER --------------------------------------------
#
# Author:     Jan D. Lanzer
# Copyright   Copyright 2024 - Jan D. Lanzer
# Email:      lanzerjan@gmail.com
#
# Date:     2024-06-11
#
# Script Name:    
#
# Script Description:
# prepare mixes of the additional studies (outside HF-core)
# we have to calculate cell type compositions for each study from the single cell 
# object. Next, we summarize to pseudobulks on sample level. 

library(tidyverse)

# load the p
studies<-list.files(path = "data/additional_studies/pb_data/", full.names = T)

# Function to read RDS file and assign it to a variable named after the file (without extension)
read_and_assign <- function(file_path) {
  # Read the Seurat object from the file
  print(file_path)
  study <- readRDS(file_path)
  
  # Extract the file name without extension
  file_name <- str_remove(basename(file_path), "\\.rds$")
  
  # Dynamically assign the Seurat object to a variable named after the file
  assign(file_name, study, envir = .GlobalEnv)
}
# Apply the function to all files in the list
x= lapply(studies, read_and_assign)

names(x)<- str_remove(basename(studies), "\\.rds$")

## add LVAD and MI study 

x$Amrute$pb <- read_csv("data/val_studies/pbulk/Armute2023_LVAD_pbulk.csv")%>%
  as.data.frame()%>% column_to_rownames('...1')%>%as.matrix()%>% t()
x$Kuppe$pb  <- read_csv("data/val_studies/pbulk/Kuppe2022_MI_pbulk.csv")%>%
  as.data.frame()%>% column_to_rownames("sample_id")%>%as.matrix()%>% t()

x$Amrute$target  <- read_csv("data/val_studies/coldata/Armute2023_LVAD_coldata.csv")%>% 
  rename(orig.ident = sample_id,
         cell_count = psbulk_n_cells)
x$Kuppe$target <- read_csv("data/val_studies/coldata/Kuppe2022_MI_coldata.csv")%>% 
  rename(orig.ident = sample_id, 
         cell_count = psbulk_n_cells)

# manually renaming some of the columns ..
x$CHD_HILL_pb$target<-x$CHD_HILL_pb$target%>% 
  select(-orig.ident)%>%
  dplyr::rename(cell_type= cell_type_translated ,
                orig.ident= orig_ident)

x$HCA_pb$target<-x$HCA_pb$target%>% 
  dplyr::rename(cell_type= cell_type_translated, 
                orig.ident = donor)

x$HCM_dimmeler_pb$target<-x$HCM_dimmeler_pb$target%>% 
  dplyr::rename(cell_type= cell_type_translated)

#calculate compositions for each sample:
comp_matrices<-lapply(x, function(study){
  print(study$target)
  comp.matrix<- study[["target"]]%>% 
    group_by(orig.ident)%>%
    filter(!is.na(cell_type))%>%
    mutate(total_cells = sum(cell_count))%>% 
    rowwise()%>%
    mutate(prop= cell_count/total_cells)%>%
    select(orig.ident, prop, cell_type)%>%
    pivot_wider(names_from = cell_type, values_from = prop)%>%
    dplyr::rename(sample_id= orig.ident)
    
})
comp_matrices

lapply(comp_matrices, colnames) %>% unlist()%>% table()
comp_matrice_combined<-do.call(rbind, comp_matrices)


#clean target file: 
lapply(x, function(y){
  
})

x$CHD_HILL_pb$target<-
x$CHD_HILL_pb$target%>% 
  mutate(heart_failure= ifelse(DEid == "CHD", "HF", "NF"), 
         disease_code= DEid)

x$HCM_dimmeler_pb$target<-
x$HCM_dimmeler_pb$target%>%
  mutate(heart_failure= ifelse(condition == "AS", "HF", "NF"),
         disease_code= condition)

x$HCA_pb$target<-x$HCA_pb$target%>% mutate(heart_failure= "NF", 
                                           disease_code= "NF")

x$pediatricDCM_pb_sample_celltype$target<-
x$pediatricDCM_pb_sample_celltype$target%>%
  mutate(disease_code = group, 
         heart_failure=  ifelse(disease_code == "CT", "NF", "HF"))

x$Sarcoidosis_pb$target<-x$Sarcoidosis_pb$target%>% 
  mutate(disease_code = disease_state ,
         heart_failure= "HF")

x$Kuppe$target

x$Amrute$target

comb_meta<-lapply(names(x), function(y){
  print(y)
  x[[y]]$target%>% select(orig.ident, cell_type, cell_count, heart_failure, disease_code)%>%
    mutate(study=y)
})%>% do.call(rbind,.)



#summarize pseudobulks
pb.summed <- lapply(names(x), function(y){
  print(y)
  gex<- x[[y]]$pb %>%as.data.frame()%>% as.matrix()%>% t()
  gex.summed <- gex%>% as.data.frame()%>% rownames_to_column("sample_id")%>%
    mutate(sample_id = gsub("_[^_]*$", "", sample_id))%>%
    group_by(sample_id)%>% summarise_all(sum)
  gex.summed%>%as.data.frame()%>%
    column_to_rownames("sample_id")%>% 
    as.matrix()%>% 
    t()
  # #remove the cell_type label to get patient ID
  # pb$sample_ids <- gsub("_[^_]*$", "", pb[,1])
  # #group by patient Id and summarise counts across cell types for each gene
  # counts.df <- pb[,-1] %>% group_by(sample_ids)%>% summarise_all(sum)
})
names(pb.summed)<- names(x)

pb.summed %>% saveRDS("data/val_studies_samplebp.rds")
comp_matrice_combined%>% write_csv("data/val_studies_composition.csv")
comb_meta%>% write_csv("data/val_studies_meta.csv")


pb= readRDS(paste0(directory, "data/HF_studiespb.rds"))

#Normalize with TPM

pb2= lapply(names(x), function(y){
  gex<- x[[y]]$pb
  bulk= Normalization(gex, method = "TPM", local = "local")
  bulk <- bulk %>% as.data.frame()%>% rownames_to_column("GeneSymbol")
  
  write.table(bulk,
              paste0(directory, "output/deconvo/mixtures/mix_", x, ".txt"),
              append = FALSE, sep = "\t", dec = ".", row.names = F, col.names = T, quote =F)
  
  return(bulk)
})
