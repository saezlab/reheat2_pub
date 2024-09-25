# HEADER --------------------------------------------
#
# Author:     Jan D. Lanzer
# Copyright   Copyright 2023 - Jan D. Lanzer
# Email:      lanzerjan@gmail.com
#
# Date:     2023-10-24
#
# Script Name:    
#
# Script Description:
# utility functions


process_HMARKS= function(path= "data/h.all.v2023.2.Hs.symbols.gmt"){
  
  ## get msig genes
  lines <- readLines(path)
  # Initialize an empty list to store gene sets
  gene_sets <- list()
  
  # Process the lines and create a list of gene sets
  for (line in lines) {
    elements <- strsplit(line, "\t")[[1]]
    gene_set_name <- elements[1]
    genes <- elements[-1]
    gene_sets[[gene_set_name]] <- genes
  }
  
  gene_sets
  
  
  net.msig= enframe(gene_sets, name = "source", value = "target")%>% 
    unnest(target)%>% 
    filter(!grepl("https", target))%>% 
    mutate(mor= 1)
}


process_MSIG2= function(path= "data/c3.all.v2023.2.Hs.symbols.gmt"){
  
  ## get msig genes
  lines <- readLines(path)
  # Initialize an empty list to store gene sets
  gene_sets <- list()
  
  # Process the lines and create a list of gene sets
  for (line in lines) {
    elements <- strsplit(line, "\t")[[1]]
    gene_set_name <- elements[1]
    genes <- elements[-1]
    gene_sets[[gene_set_name]] <- genes
  }
  
  gene_sets
  
  
  net.msig= enframe(gene_sets, name = "source", value = "target")%>% 
    unnest(target)%>% 
    filter(!grepl("https", target))%>% 
    mutate(mor= 1)
}




plot_exp= function(genes,pb_long){
  lapply(genes, function(gene){
    mod.df= pb_long %>% filter(feature ==gene)
    p.exp= ggplot(mod.df, aes(x= batch, y= (value),color = leiden_0.7, fill = heart_failure))+
      geom_boxplot()+
      scale_fill_manual(values= c("white", "darkgrey"))+
      ggtitle(paste0("Expression of ", gene))+
      theme_cowplot()+
      theme(axis.text.x = element_text(angle =45, hjust= 1))
    
    
  }) %>% plot_grid(plotlist = .)
  
}


get_pk_nets<- function(){
  require(msigdbr)
  require(decoupleR)
  
  net_ct <- get_collectri(organism='human', split_complexes=FALSE)
  net <- get_dorothea()
  net_prog <- get_progeny()
  
  immune = read_csv("data/prior_knowledge/immune.csv")
  ILs <- immune$source %>% unique()
  ILs <- enframe(ILs, value = "target", name= "source")%>%
    mutate(source= "Cytokine_espression", weight=0.5 )%>%
    select(target, source, weight)
  
  net_immune<- rbind(immune, 
                     ILs)
  
  net_c2 <- msigdbr(category = "C2" )
  net_c2 <- net_c2%>% 
    filter(gs_subcat != "CGP")%>%
    dplyr::rename(target = gene_symbol,source=gs_name)%>% distinct(source, target)
  
  net_c5 <- msigdbr(category = "C5")%>% 
    filter(grepl("GO", gs_subcat))%>%
    dplyr::rename(target = gene_symbol,source=gs_name)%>% distinct(source, target)
  
  net_naba<-net_c2 %>% filter(grepl("NABA", source))%>% distinct(source, target)
  net_h <- msigdbr(category = "H")%>% 
    dplyr::rename(target = gene_symbol,source=gs_name)%>% distinct(source, target)
  
  net_list = list("tfs" = net,
                  "tfs_ct"= net_ct,
                  "cytokines"= net_immune,
                  "progeny"= net_prog, 
                  "naba" = net_naba, 
                  "hmarks"= net_h, 
                  "GO_MF"= net_c5%>% filter(grepl("GOMF", source)), 
                  "GO_BP"= net_c5%>% filter(grepl("GOBP", source)),
                  "cannonicalPways"= net_c2)
}

