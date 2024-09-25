# HEADER --------------------------------------------
#
# Author:     Jan D. Lanzer
# Copyright   Copyright 2024 - Jan D. Lanzer
# Email:      lanzerjan@gmail.com
#
# Date:     2024-04-10
#
# Script Name:    
#
# Script Description:
# create bulk mixtures from reheat2 
library(tidyverse)

pb.files= list.files("data/pb_data", full.names = T)
studies= map(str_split(list.files("data/pb_data", full.names = F),"_"), function(x) x[[1]])%>% unlist()

pb.list = map(pb.files, function(x){
  (read.csv(x,header = T, sep = ",", quote = ""))
})
names(pb.list) <- studies
#cell_dic <- read.csv("output/celltype_dictionary_HCA_ReHeaT2.csv")

## we will summarise each patient to count vector 
pb.summed <- lapply(pb.list, function(pb){
  #remove the cell_type label to get patient ID
  pb$sample_ids <- gsub("_[^_]*$", "", pb[,1])
  #group by patient Id and summarise counts across cell types for each gene
  counts.df <- pb[,-1] %>% group_by(sample_ids)%>% summarise_all(sum)
})

pb.matrix <- lapply(pb.summed, function(pb){
  pb %>% as.data.frame()%>%
    column_to_rownames("sample_ids")%>% 
    as.matrix()%>% 
    t()
})

pb.matrix %>% saveRDS("data/all_HF_studies_smaple_pb.rds")


pb= readRDS(paste0( "data/HF_studiespb.rds"))

pb2 = readRDS("data/all_HF_studies_smaple_pb.rds")
pb2
# get proportions  
sc<- read_csv("data/sample_comps.csv")
cm<- read_csv("data/comps_meta.csv")
    