# HEADER --------------------------------------------
#
# Author:     Jan D. Lanzer
# Copyright   Copyright 2024 - Jan D. Lanzer
# Email:      lanzerjan@gmail.com
#
# Date:     2024-02-05
#
# Script Name:    
#
# Script Description:
# process HF progressive data
library(tidyverse)

summary_mtx <- read.csv("data/mouse_data/GSE1200064/GSE120064_TAC_clean_cell_info_summary.txt", 
                        sep='\t')
summary_mtx

mtx <- read.csv("data/mouse_data/GSE1200064/TAC_raw_umi_matrix.csv", 
                        sep=',')
head(mtx)
mtx.df <- mtx %>% as.data.frame  %>% column_to_rownames("X")

summary_mtx$CellID %in% colnames(mtx.df)
summary_mtx <- summary_mtx%>% drop_na()

summary_mtx <-   summary_mtx%>% 
  mutate(celltype = ifelse(CellType =="MP", "Myeloid", NA), 
         celltype = ifelse(CellType =="FB", "Fib",celltype),
         celltype = ifelse(CellType =="EC", "Endo", celltype),
         celltype = ifelse(CellType =="T", "Lymphoid", celltype),
         celltype = ifelse(CellType =="GN", "Myeloid" ,celltype), 
         celltype = ifelse(CellType =="CM", "CM" ,celltype))
unique(g.load$ctype)

samps <- unique(summary_mtx$sample)
ctypes <- unique(summary_mtx$celltype)
x= samps[1]
y= ctypes[1]


##summarize to pseudo bulks<- 
x= map(samps, function(x){
  
 pb <-  sapply(ctypes, function(y){
    cells <- summary_mtx %>% filter(sample == x, 
                                    celltype == y)%>%
      pull(CellID)
    
    sample_1 <- mtx[, cells ]
    print(c(x, y))
    print(dim(sample_1))
    if(is.null(dim(sample_1))){return(rep(0, dim(mtx)[1]))}
    pb <- rowSums(sample_1)
    #rownames(pb)<- rownames(mtx.df)
    return(pb)
    
    
    
})
  rownames(pb)<- rownames(mtx.df)
  colnames(pb)<- paste0(x, sep ="_" ,ctypes)
  return(pb)
})

pb.matrix <- do.call(cbind, x)+
  geom_point()

## count number of cells
cell_count = sapply(samps, function(x){
  
  sapply(ctypes, function(y){
    cells <- summary_mtx %>% filter(sample == x, 
                                    celltype == y)%>%
      pull(CellID)
    
    return(length(cells))
    
    
  })
})

target_file <- summary_mtx %>%
  group_by(celltype, sample)%>% 
  summarise(n= n())%>%
  left_join(summary_mtx%>% distinct(celltype, sample,  condition))%>%
  mutate(pb_id = paste0(sample, sep="_", celltype))

saveRDS(list("pb"= pb.matrix, 
             "target"= target_file), 
        "data/mouse_data/GSE1200064/GSE1200064_pb.rds")



target_file%>% print(n=100)

target_file %>% 
  ggplot(aes(x= condition,y= n,  fill = sample))+
  geom_point(aes(color=celltype))+
  geom_boxplot()

# filter and normalize the pb ---------------------------------------------
pbs<- readRDS( "data/additional_studies/mouse_data/GSE1200064/GSE1200064_pb.rds")

pbs$target %>% pull(sample) %>% unique() %>% length()

