# HEADER --------------------------------------------
#
# Author:     Jan D. Lanzer
# Copyright   Copyright 2024 - Jan D. Lanzer
# Email:      lanzerjan@gmail.com
#
# Date:     2024-04-04
#
# Script Name:    
#
# Script Description:
# We will test whether the preselection of genes towards markers that are not
# not subjected to strong molecular changes, will improve deconvolution results

set.seed(20)

library(Matrix.utils)
library(tidyverse)
library(Matrix)
library(Seurat)
library(ComplexHeatmap)
library(cowplot)

directory = "./"

source(file = "celltype_deconvolution/CIBERSORT.R")
source(paste0(directory, "celltype_deconvolution/utils_decon.R"))

## load data: 

# sc ref: 
seu = readRDS(paste0(directory, "data/additional_studies/HCA_seu_filtered_nuc.rds"))

print(colnames(seu@meta.data))
print(unique(seu@meta.data$cell_type))
print(wanted_cells)

# make sure to only use samples of interest in reference
#seu = Seurat::subset(seu, cell_type %in% wanted_cells)
seu= subset(seu, subset = (region %in% c("LV", "AX)") & 
                             cell_type %in% wanted_cells &
                             source == "Nuclei" & 
                             Used == "Yes" &
                             nCount_RNA >200 &  
                             nFeature_RNA > 500)
)
gex.count<- seu@assays$RNA$counts
print(unique(seu@meta.data$cell_type))
print("input loaded")

#Idents(seu)= "cell_type"

# consensus marker:
cons<- read.csv("data/consensus_markers.csv")%>%
  as_tibble()

marker.genes <- cons %>% 
  group_by(cell_type)%>%
  filter(adj_Fisher_p< 10e-50, mean_LFC>=2)%>%
  #top_n(300, mean_LFC)%>% 
  pull(gene)

# compositional info on genes:
comp_genes = read.csv("data/reheat2_stats_bulk.csv")%>% 
  filter(fisher_pvalue<0.05)
length(unique(comp_genes$gene))
  

# desing different features list based on mofa model

feature_spaces = list("allGenes"= comp_genes %>% pull(gene), 
                      #"comp"= comp_genes %>% filter(scell_consistency == "compositional")%>% pull(gene), 
                      "mol"= comp_genes %>% filter(scell_consistency %in% c("molecular"))%>% pull(gene), 
                      "comp"= comp_genes %>% filter(scell_consistency %in% c("comp and mol", "compositional"))%>% pull(gene)
                      #"unknown" = comp_genes %>% filter(scell_consistency %in% c("unknown"))%>% pull(gene)
) 
#add marker genes that are not annotated:
table(marker.genes %in% feature_spaces$allGenes)
#feature_spaces$markers = marker.genes[!marker.genes %in% feature_spaces$allGenes]


lapply(feature_spaces, length)

feature_spaces_reheat = list("comp"= comp_genes %>% 
                               filter(scell_consistency %in% c("comp and mol", "compositional")) %>% 
                               pull(gene)
)

# prepare a sig matrix ----------------------------------------------------

#' @description function to create signature matrixes for cibersort input by sampling random cells
#' n times from the referenced single cell count matrix and summing the counts to cell lineage pseudobulks
#'# then TPM normalizing and finally subsetting features to the ones in marker.genes AND all different sets
#'from feature_spaces
#' @param n, number of signature matrices with random sampling

create_sigs = function(n= 1,
                       marker.genes, 
                       feature_spaces, 
                       gex.count, 
                       path_to_sig_matrices, 
                       create_unfiltered=T){
  
  # loop until n to create different signature matrices, to test robustness of deconv
  # results
  for (i in 1:n){
    set.seed(i+1)
    cell_ids= colnames(gex.count)
    cell_ids_sample= sample(cell_ids, size = 0.3*length(cell_ids), replace = F)
    groups = seu@meta.data %>% select(cell_type)
    groups <- groups[cell_ids_sample, ]
    
    gex.count2 = gex.count[ , cell_ids_sample]
    dim(gex.count2)
    
    pb <- aggregate.Matrix(t(gex.count2),
                           groupings = groups,
                           fun= "sum")
    # feature selection
    # this is the tricky part, we will select genes based on 
    # 1. min expression (done by filtering seurat obj above)
    # 2. consensus marker
    # 3. filter further for comp and mol etc. 
    
    pb <- t(pb)
    
    # select genes if at least 50 counts
    pb <- pb[apply(pb, 1, max)>10 , ]
    
    bulk= Normalization(pb, method = "TPM", local = "local")
    
    table(marker.genes %in% rownames(bulk))
    
    bulk <- bulk +1
    #x= "comp"
    map(names(feature_spaces), function(x){
      print(x)
      #filter for marker genes and captured in bulk
      marker.genes.filt <- intersect(marker.genes, rownames(bulk))
      #filter for genes from mofa
      marker.genes.filt <- intersect(marker.genes.filt, feature_spaces[[x]])
      
      print(length(marker.genes.filt))
      
      bulk_filt <- as.data.frame(bulk[marker.genes.filt,])%>%
        rownames_to_column("GeneSymbol")
      
      write.table(bulk_filt,
                  paste0(path_to_sig_matrices, "Matrix_",x,"_",i, ".txt"),
                  append = FALSE, sep = "\t", dec = ".", row.names = F, col.names = T, quote = F)
    })
    
    # here we add a signature matrix without any filtering
    marker.genes.filt <- intersect(marker.genes, rownames(bulk))
    #remove any deregulated marker
    marker.genes.filt <- marker.genes.filt[!marker.genes.filt %in% unlist(feature_spaces)]
    print(length(marker.genes.filt))
    
    bulk_filt <- as.data.frame(bulk[marker.genes.filt,])%>%
      rownames_to_column("GeneSymbol")
    
    write.table(bulk_filt,
                paste0(path_to_sig_matrices, "Matrix_not_reg_",i, ".txt"),
                append = FALSE, sep = "\t", dec = ".", row.names = F, col.names = T, quote = F)
  }
}

# signature matrices for deconvo testing: 
create_sigs(n= 3, marker.genes = marker.genes, 
            feature_spaces, 
            gex.count, 
            path_to_sig_matrices = "output/deconvo/signature_Ms/")

create_sigs(n= 3, marker.genes = marker.genes, feature_spaces_reheat, gex.count, 
            path_to_sig_matrices = "output/deconvo/signature_Ms_reheat/")
#groups= seu@meta.data[, c("cell_type")]
#groups$sample <- str_replace_all(groups$sample, "_", ".") 
#groups$cell_type <- str_replace_all(groups$cell_type, "_", ".") 


## here we will test ctype marker genes for their regulation in HF

# merge annotations
marker.genes.df <- cons %>% 
  group_by(cell_type)%>%
  filter(adj_Fisher_p< 10e-50, mean_LFC>=2)%>%
  left_join(comp_genes %>% 
              filter(fisher_pvalue < 0.05)%>%
              select(gene,scell_consistency, fisher_pvalue))
  

# replace 0 pvals with next lowest pval
min.p = marker.genes.df%>%filter(adj_Fisher_p>0)%>% pull(adj_Fisher_p) %>% min()

# test wether marker genes that are also regulated in mofa have different
# pval distribtution means
marker.genes.df<- marker.genes.df%>% 
  mutate(mofa= ifelse(is.na(scell_consistency), "no", "yes"),
         p_corr= ifelse(adj_Fisher_p==0, min.p, adj_Fisher_p), 
         lop = -log10(p_corr)
         )

marker.genes.df%>%
  
  group_by(cell_type)%>%
  nest()%>%
  mutate(pval= map(data, function(x){
    res= t.test(lop ~ mofa, data = x)
    res$p.value
  }))%>% unnest(pval)%>%
  mutate(stat= map(data, function(x){
    res= t.test(lop ~ mofa, data = x)
    res$statistic
  }))%>% unnest(stat)#%>% mutate(pval= unlist(pval))

# plot 
p.bars_dist= marker.genes.df%>%
  ggplot(aes(x= cell_type , fill =scell_consistency))+
  geom_bar(color="black")+
  #ggtitle("cell type marker by regulation in HF")+
  theme(axis.text.x = element_text(angle= 45, hjust= 1))+
  labs(x= "")

p.marker.sig = marker.genes.df%>%
  ggplot(aes(y= lop, x =scell_consistency))+
  geom_violin()+
  facet_wrap(~cell_type)+
  theme(axis.text = element_text(angle= 90, hjust= 1))+
  labs(x= "", y= "log(p-val)")
p.bars_dist  
p.marker.sig

plot_grid(p.bars_dist, p.marker.sig)

## calculate top 50% of cell type whether they are regulated in HF
x= cons %>% 
  group_by(cell_type)%>%
  filter(adj_Fisher_p< 0.01, mean_LFC>=2)%>%
  left_join(comp_genes %>% 
              filter(fisher_pvalue <0.1) %>%
              select(gene,scell_consistency, fisher_pvalue))%>% 
  mutate(scell_consistency = ifelse(scell_consistency=="comp and mol", "compositional", scell_consistency))%>%
  group_by(cell_type)

min.p = x%>%filter(adj_Fisher_p>0)%>% pull(adj_Fisher_p) %>% min()

x<-x%>% 
  mutate(mofa= ifelse(is.na(scell_consistency), "no", "yes"),
         p_corr= ifelse(adj_Fisher_p==0, min.p, adj_Fisher_p), 
         lop = -log10(p_corr)
  )%>%
  arrange(desc(lop))%>%
  #mutate(scell_consistency = ifelse(is.na(scell_consistency), "not_regulated", "regulated"))
  mutate(scell_consistency = ifelse(is.na(scell_consistency), "not_regulated", scell_consistency))


marker.reg <- prop.table(table(x$cell_type, x$scell_consistency),margin= 1 )
colnames(marker.reg)<- c( "dereg_comp", "dereg_mol", "unreg", "dereg_unkown")
p_state_marker_composition<- marker.reg %>% 
  as.data.frame() %>%
  mutate(Var2= factor(Var2, levels= sort(colnames(marker.reg))))%>%
  ggplot(., aes(x= Var1, y= Freq, fill = Var2))+
  geom_col(col= "black")+
  theme_cowplot()+
  scale_fill_brewer(type = "qual", direction = -1)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  theme(axis.text.x = element_text(angle= 90,vjust= 0.5,  hjust= 1))+
  labs(fill = "Cell type markers\nin HF are", y= "", 
       x= "")
p_state_marker_composition
1-min(marker.reg[,"unreg"])
1-max(marker.reg[,"unreg"])

pdf("output/figures/deconv_compare_ctype_markerscomp.pdf", 
    width= 4, 
    height= 3.5)
p_state_marker_composition
dev.off()

pdf("output/figures/deconv_compare_ctype_markers.pdf", 
    width= 9, 
    height= 6)
print(plot_grid(p.bars_dist, p.marker.sig))
dev.off()

# Prepare Mixtures -----------------------------------------------------------------

#load psuedobulked data matrices
pb_core= readRDS(paste0(directory, "data/all_HF_studies_smaple_pb.rds"))
pb_val= readRDS(paste0(directory, "data/val_studies_samplebp.rds"))
pb= append(pb_core, pb_val)

#Normalize with TPM
pb2= lapply(names(pb), function(x){
  bulk= Normalization(pb[[x]], method = "TPM", local = "local")
  bulk <- bulk %>% as.data.frame()%>% rownames_to_column("GeneSymbol")
  
  write.table(bulk,
              paste0(directory, "output/deconvo/mixtures/mix_", x, ".txt"),
              append = FALSE, sep = "\t", dec = ".", row.names = F, col.names = T, quote =F)
  
  return(bulk)
})


## prepare the Mixtures of real bulk

#to do : update with expanded resource!
##meta_heart<-readRDS("data/MetaHeart_RNA_counts.rds")
meta_heart <-readRDS("output/reheat1/Metaheart_counts2023.rds")

lapply(names(meta_heart), function(x){
  bulk= Normalization(meta_heart[[x]]$gex, method = "TPM", local = "local")
  bulk <- bulk %>% as.data.frame()%>% rownames_to_column("GeneSymbol")
  
  write.table(bulk,
              paste0(directory, "output/deconvo/mixtures_reheat/", x, ".txt"),
              append = FALSE, sep = "\t", dec = ".", row.names = F, col.names = T, quote =F)
  
  return(bulk)
})


# deconvolute Mixtures -------------------------------------------------------------

# load signature matrices
sig.matrices_paths <- list.files("output/deconvo/signature_Ms", full.names = T)
sig.matrices_names <-str_replace_all(sig.matrices <- list.files("output/deconvo/signature_Ms", full.names = F), ".txt", "")

#optional 
#sig.matrices_names<-sig.matrices_names[grepl("markers", sig.matrices_names)]
#sig.matrices_paths<-sig.matrices_paths[grepl("markers", sig.matrices_paths)]
#load mixture data
mixtures_paths <- list.files("output/deconvo/mixtures", full.names = T)
mixtures_names <- str_replace_all(list.files("output/deconvo/mixtures", full.names = F),".txt", "")

mixtures_paths<-mixtures_paths[!grepl("mix_Reichart2022", mixtures_paths)]
mixtures_paths<-mixtures_paths[!grepl("HCA", mixtures_paths)]
mixtures_names<-mixtures_names[!grepl("HCA", mixtures_names)]
mixtures_names<-mixtures_names[!grepl("mix_Reichart2022", mixtures_names)]

# will run cibersort with all sig.matrices on all mixtures:
total_res = map2(mixtures_paths, mixtures_names, function(m.p, m.n){
  res= map2(sig.matrices_paths, sig.matrices_names, function(x.p, x.n){
    
    print(paste0("deconvolute with ",x.n , " mixture of ", m.n))
    # run cibesort here
    res <- CIBERSORT(sig_matrix = x.p,
                     mixture_file = m.p , QN = F, perm = 0)
  })
  
  names(res)<- sig.matrices_names
  return(res)
})

names(total_res)= mixtures_names

saveRDS(total_res, "output/deconvo/cibersort_results_bench.rds")


# deconvolute ReHeaT1 -----------------------------------------------------

# load signature matrix
#sig.matrices_paths <- list.files("output/deconvo/signature_Ms/", full.names = T)
#sig.matrices_names <- str_replace_all(sig.matrices <- list.files("output/deconvo/signature_Ms_reheat", full.names = F), ".txt", "")

sig.matrices_paths <- list.files("output/deconvo/signature_Ms", full.names = T)
sig.matrices_names <-str_replace_all(sig.matrices <- list.files("output/deconvo/signature_Ms", full.names = F), ".txt", "")

sig.matrices_paths_reheat<- sig.matrices_paths[grepl("comp", sig.matrices_paths)]
sig.matrices_names_reheat<- sig.matrices_names[grepl("comp", sig.matrices_paths)]

#load mixture data
mixtures_paths <- list.files("output/deconvo/mixtures_reheat", full.names = T)
mixtures_names <- str_replace_all(list.files("output/deconvo/mixtures_reheat", full.names = F),".txt", "")

# 
total_res = map2(mixtures_paths, mixtures_names, function(m.p, m.n){
  res= map2(sig.matrices_paths_reheat, sig.matrices_names_reheat, function(x.p, x.n){
    
    print(paste0("deconvolute with ",x.n , " mixture of ", m.n))
    # run cibesort here
    res <- CIBERSORT(sig_matrix = x.p,
                     mixture_file = m.p , QN = F, perm = 0)
  })
  
  names(res)<- sig.matrices_names_reheat
  return(res)
})

names(total_res)= mixtures_names

saveRDS(total_res, "output/deconvo/cibersort_results_reheat.rds")


