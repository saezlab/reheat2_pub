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
# process hill et al chd data

library(Seurat)
library(tidyverse)

directory= "data/additional_studies/hill_chd/SCP1852/expression/6282f075771a5b02c609b8b5/"


source("process_studies/utils_process.R")

expression_matrix <- ReadMtx(
  mtx = paste0(directory,"AllNuclei_snRNA_counts.mtx.gz"), 
  features = paste0(directory,"AllNuclei_snRNA_counts_rownames.txt.gz"),
  cells =paste0(directory,"AllNuclei_snRNA_counts_colnames.txt.gz"),
  feature.column = 1
)
  
seu <- CreateSeuratObject(counts = expression_matrix)

rm(expression_matrix)

# seu@meta.data %>% 
#   ggplot(aes(x= nFeature_RNA, y=nCount_RNA ))+
#   geom_point()
# already filtered 

meta.data <- read_csv(paste0(directory,"../../metadata/AllNuclei_snRNA_metadata.csv"))

table(meta.data$NAME %in% colnames(seu))
table( colnames(seu) %in% meta.data$NAME)

meta.data.sorted<-meta.data[match( colnames(seu),meta.data$NAME),]

x= cbind(seu@meta.data, meta.data.sorted)
table(rownames(x)== x$NAME)

seu <- AddMetaData(seu,metadata = meta.data.sorted )

## now summarize cell type labels: 
seu@meta.data<- seu@meta.data%>% 
  mutate(cell_type_translated = ifelse(MainCellType=="CM", "CM", NA), 
         cell_type_translated = ifelse(MainCellType %in% c("Mac", "Mast"), "Myeloid", cell_type_translated), 
         cell_type_translated = ifelse(MainCellType %in% c("Tcells", "Mast"), "Lymphoid", cell_type_translated),
         cell_type_translated = ifelse(MainCellType %in% c("CF"), "Fib", cell_type_translated),
         cell_type_translated = ifelse(MainCellType %in% c("SMC"), "vSMCs", cell_type_translated),
         cell_type_translated = ifelse(MainCellType %in% c("PeriC"), "PC", cell_type_translated),
         cell_type_translated = ifelse(MainCellType %in% c("Endo"), "Endo", cell_type_translated),
         )

seu@meta.data%>% distinct(labID, orig.ident, orig_ident, Diagnosis, biosample_id,organ__ontology_label,region, batch_indices,vers10X,organ, DEid, cell_type_translated)#%>% filter(region =="LV")
seu@meta.data%>% as_tibble()%>% select(-nCount_RNA,-nFeature_RNA, -NAME)%>% distinct()

pb.matrix<- get_pseudobulk(seu,
                           cluster_col = "cell_type_translated", 
                           sample_col = "orig_ident", 
                           exclude_clusters = NA)



target_file <- seu@meta.data %>%
  group_by(cell_type_translated, orig.ident)%>% 
  summarise(cell_count= n())%>%
  left_join(seu@meta.data%>% 
              distinct(labID, orig.ident, orig_ident,age, gender,  Diagnosis,procedure, region, batch_indices,vers10X,organ, DEid, cell_type_translated)
            )%>%
  mutate(pb_id = paste0(orig_ident, sep="_", cell_type_translated))

saveRDS(list("pb"= pb.matrix, 
             "target"= target_file), 
        "data/additional_studies/pb_data/CHD_HILL_pb.rds")

