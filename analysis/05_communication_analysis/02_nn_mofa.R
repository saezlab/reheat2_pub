# HEADER --------------------------------------------
#
# Author:     Jan D. Lanzer
# Copyright   Copyright 2023 - Jan D. Lanzer
# Email:      lanzerjan@gmail.com
#
# Date:     2023-11-16
#
# Script Name:    
#
# Script Description:
# explore run Nichenet 
#devtools::install_github("saeyslab/nichenetr")

library(nichenetr)
library(cowplot)
library(tidyverse)
library(igraph)

g.loads= read_csv("output/mofa/gene_loadings.csv")

liana_df <- readRDS( "output/communication/lr_meta_mofa.rds")

#updated PK

lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
#ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))

#
ligand_target_matrix = readRDS("~/Downloads/ligand_target_matrix_nsga2r_final.rds")

#use liana pre run for nichenet ------------------------------------------
  
  run_cross_nichnet_w_liana=function(g.loads, 
                                     neg_scores= T,
                                     ctypes,
                                     top_n_geneset_oi= 100, 
                                     fact= "Factor1",
                                     explain_top_genes= 0.1, 
                                     liana_df,
                                     prop_background= 0.3){
    g.loads <- g.loads%>% 
      filter(Factor == fact)
    
    # pariwise mapping 
    #receiver_c = "CM"
    #sender_c ="Endo"
    df.niche= map(ctypes, function(sender_c){
      map(ctypes, function(receiver_c){
        
        ##1
        # select the target genes (receiver), that will be explained by ligands(sender)
        # we run this for binary classes associated with a neg factor scores or positive factor scores
        # (neg_scores T or F)
        if(neg_scores){
          
          expressed_genes_sender <- g.loads %>% 
            filter(ctype == sender_c,
                   (value)<0
            )%>% 
            pull(feature)
          
          expressed_genes_receptor <- g.loads %>% 
            filter(ctype == receiver_c, 
                   (value)<=0
            )%>% 
            pull(feature)
          
          ligands_sig <- liana_df %>% 
            filter(group =="HF" & sig=="s" & sender == sender_c & receiver == receiver_c)%>% 
            pull(source_genesymbol)%>% unique()%>% sort()
          
          receptors_sig <-liana_df %>% 
            filter(group =="HF", 
                   sig=="s",
                   sender == sender_c, 
                   receiver == receiver_c)%>% 
            pull(target_genesymbol)%>% unique()
          
        }else{
          
          expressed_genes_sender <- g.loads %>% 
            filter(ctype == sender_c,
                   (value)>0
            )%>% 
            pull(feature)
          
          expressed_genes_receptor <- g.loads %>% 
            filter(ctype == receiver_c, 
                   (value)>=0)%>% 
            pull(feature)
          
          ligands_sig <- liana_df %>% 
            filter(group =="NF", 
                   sig=="s",
                   sender == sender_c, 
                   receiver == receiver_c)%>% 
            pull(source_genesymbol)%>% unique()%>% sort()
          ligands_sig
          receptors_sig <-liana_df %>% 
            filter(group =="NF", 
                   sig=="s",
                   sender == sender_c, 
                   receiver == receiver_c)%>% 
            pull(target_genesymbol)%>% unique()
        }
        
        
        geneset_oi <- g.loads %>% 
          filter(ctype == receiver_c) %>%
          mutate(valuem= abs(mean(value) - value))%>% 
          arrange(desc(valuem))%>% 
          slice_head(prop= explain_top_genes)%>%
          pull(feature)
        
        # for the expression genes of the receiver we will consider 
        # features around 0 
        
        expressed_genes_receiver <- g.loads %>% 
          filter(ctype == receiver_c) %>%
          mutate(valuem= abs(mean(value) - value))%>% #here we identify the geneset oi as 
          arrange(valuem)%>% 
          slice_head(prop= prop_background)%>%
          pull(feature)
        
        print(paste("sender cell:", sender_c, "expressed genes:", length(expressed_genes_sender)))
        print(paste("receiver cell:", receiver_c, "expressed genes:", length(expressed_genes_receiver)))
        
        ##2  
        
        background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
        background_expressed_genes <- background_expressed_genes[!background_expressed_genes %in% geneset_oi]
        head(background_expressed_genes)
        
        ##3 
        # ligands = lr_network %>% pull(from) %>% unique()
        # expressed_ligands = intersect(ligands,expressed_genes_sender)
        # 
        # receptors = lr_network %>% pull(to) %>% unique()
        # expressed_receptors = intersect(receptors,expressed_genes_receptor)
        # 
        # lr_network_expressed = lr_network %>% filter((from %in% ligands_sig) & (to %in% receptors_sig))
        # head(lr_network_expressed)
        # 
        # potential_ligands = lr_network_expressed %>% pull(from) %>% unique()
        # head(potential_ligands)
        potential_ligands = intersect(colnames(ligand_target_matrix), ligands_sig)
        print(paste("ligands n : ", length(potential_ligands)))
        if(length(potential_ligands)==0 | length(receptors_sig)== 0){
          return(NULL)
        }
        ##4 calc ligand activities:
        
        ligand_activities = predict_ligand_activities(geneset = geneset_oi, 
                                                      background_expressed_genes = background_expressed_genes, 
                                                      ligand_target_matrix = ligand_target_matrix, 
                                                      potential_ligands = potential_ligands)
        ligand_activities <- ligand_activities %>% 
          arrange(-aupr_corrected) %>%
          mutate(sender= sender_c, 
                 receiver = receiver_c, 
                 nbackground= length(background_expressed_genes), 
                 nsender= length(expressed_genes_sender),
                 nreceiver= length(expressed_genes_receiver), 
                 expl_genes= list(geneset_oi))
        
        active_ligand_target_links_df = ligand_activities$ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
        # 
        # k = 3 # 3-fold
        # n = 2 # 2 rounds
        # best_upstream_ligands = ligand_activities %>% top_n(50, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)
        # pemt_gene_predictions_top30_list = seq(n) %>% 
        #   lapply(assess_rf_class_probabilities, folds = k, geneset = geneset_oi, background_expressed_genes = background_expressed_genes, 
        #                                                      ligands_oi = best_upstream_ligands, 
        #                                                      ligand_target_matrix = ligand_target_matrix)
        # 
        #target_prediction_performances_cv = pemt_gene_predictions_top30_list %>% lapply(classification_evaluation_continuous_pred_wrapper) %>% bind_rows() %>% mutate(round=seq(1:nrow(.)))
        active_ligand_target_links_df = potential_ligands %>% 
          lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% 
          bind_rows()
        
        ligand_activities <- ligand_activities %>% left_join(active_ligand_target_links_df, by= join_by(test_ligand == ligand))
        
        
        # active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df,
        #                                                                  ligand_target_matrix = ligand_target_matrix, 
        #                                                                  cutoff = 0.1)
        # 
        # order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
        # order_targets = active_ligand_target_links_df$target %>% unique()
        # vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
        # 
        # p_ligand_target_network = vis_ligand_target %>% 
        #   make_heatmap_ggplot("Prioritized CAF-ligands","p-EMT genes in malignant cells",
        #                       color = "purple",
        #                       legend_position = "top",
        #                       x_axis_position = "top",
        #                       legend_title = "Regulatory potential") + 
        #   #scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01)) + 
        #   theme(axis.text.x = element_text(face = "italic"))
        # 
        # p_ligand_target_network
        # 
        
        
      })%>% 
        do.call(rbind, .)
    })%>% do.call(rbind,.)
    
  }

df.hf= run_cross_nichnet_w_liana(g.loads = g.loads,
                                 neg_scores = T,
                                 ctypes= unique(g.loads$ctype),
                                 top_n_geneset_oi =  0.2, 
                                 fact = "Factor1",
                                 explain_top_genes = 0.2, 
                                 liana_df = liana_df%>% mutate(sig ="s"),
                                 prop_background= 0.5                  
)


df.ct= run_cross_nichnet_w_liana(g.loads = g.loads,
                                 neg_scores = F,
                                 ctypes= unique(g.loads$ctype),
                                 top_n_geneset_oi =  0.2, 
                                 fact = "Factor1",
                                 explain_top_genes = 0.2, 
                                 liana_df = liana_df%>% mutate(sig ="s"),
                                 prop_background= 0.5
)


df.comb = df.hf %>%
  mutate(group = "hf")%>% rbind( df.ct %>% mutate(group = "ct"))%>% 
  mutate(Factor= "Factor1")

df.comb %>% 
  saveRDS("output/communication/nn_meta_mofa_fact_liana.rds")
