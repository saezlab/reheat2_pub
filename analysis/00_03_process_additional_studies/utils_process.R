
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
      cells <- rownames(seu@meta.data %>% as.data.frame() %>%dplyr::filter(.data[[sample_col]] == x, 
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


cell_dic<-read_csv("output/celltype_dictionary_HCA_ReHeaT2.csv")
