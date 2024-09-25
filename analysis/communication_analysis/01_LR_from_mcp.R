library(tidyverse)
library(OmnipathR)
library(liana)
library(purrr)
library(magrittr)
library(igraph)
library(cowplot)

g.load= read.csv("output/mofa/gene_loadings.csv") %>% as_tibble()

consensus_omni <- select_resource("Consensus")[[1]] %>%
  glimpse


# calc_LR_interactions ----------------------------------------------------

LR_interactions_from_mofacell= function(g.load,
                                       
                                        consensus_omni,
                                        Factor_oi = "Factor1"
){
  
  set.seed(seed)
  ctypes = unique(g.load$ctype)
  
  df_res<- lapply(ctypes, function(x){
    print("For Sender: ")
    lapply(ctypes, function(y){
      
      print(paste0((c(x, y)), collapse = "_>_")) 
      
      # subset mofa loadings to sender and receiver cell type and factor of interest
      # these will be used to calculate the real interaction score later and the
      # to sample a null
      
      df_sender= g.load%>% 
        filter(ctype ==x)%>%
        filter(Factor == Factor_oi)%>%
        rename(source_genesymbol = feature,
               value_ligand= value)%>%
        select(value_ligand,source_genesymbol)
      
      df_receiver=g.load%>% 
        filter(ctype ==y)%>%
        filter(Factor == Factor_oi)%>%
        rename(target_genesymbol = feature,
               value_receptor= value)%>%
        select(value_receptor,target_genesymbol)
      
      ## with PK on LR
      df <- consensus_omni%>% 
        left_join(df_sender, 
                  by= "source_genesymbol")%>%
        left_join(df_receiver, 
                  by= "target_genesymbol")%>%
        select(-source, 
               -target, 
               -sources, 
               -references, 
               -category_intercell_source, 
               -database_intercell_source, 
               -category_intercell_target,
               -database_intercell_target
        )
    
      ## Calculate REAL interaction scores by joining sender and receiver loadings
  
      df <- df%>%
        filter(sign(value_ligand) == sign(value_receptor))%>%
        filter(abs(value_ligand)>0.1 & abs(value_receptor)>0.1)%>%
        rowwise()%>% 
        mutate(
          pair= paste0((c(x, y)), collapse = "_>_"),
          interact = sum(value_ligand, value_receptor)/2,
          pair_2= pair)%>%
        separate(pair_2, into = c("sender", "receiver"), sep = "_>_")%>%
        mutate(Factor = Factor_oi)%>%
        ungroup()
      
      #hist(df$interact)
      
      return(df)
    })%>% do.call(rbind, .)
  })%>% do.call(rbind, .)
  
  return(df_res)
}


# calculate LR data frame -------------------------------------------------



df_mcp1 <- LR_interactions_from_mofacell(g.load = g.load,
                                   
                                    Factor_oi = "Factor1",
                                    consensus_omni = consensus_omni
)

df_mcp1%>% arrange(abs(interact))

df_mcp2 <- LR_interactions_from_mofacell(g.load = g.load, 
                                     Factor_oi = "Factor2",
                                     consensus_omni = consensus_omni
)


df<- rbind(df_mcp1, df_mcp2)
df<- df%>% 
  #mutate(sig= ifelse(p.adj<0.2, "sig", "ns"))%>%
  mutate(group= ifelse(interact >0, "NF", "HF"))

df %>% saveRDS(., "output/communication/lr_meta_mofa.rds")
df <- readRDS("output/communication/lr_meta_mofa.rds")

# PLOT --------------------------------------------------------------------

#summmarize to cell type level 
df.graph <- df%>% 
  select( receiver, sender,everything())%>%
  group_by(sender, receiver, group, Factor)%>%
  summarise(n= n(), 
            w= sum(interact, na.rm = T))%>%
  group_by(group,)%>% 
  mutate(scaled_weight= as.numeric(scale(w)))

# plot hmap of n LR 

pls_n= map(unique(df.graph$Factor), function(x){
  p.interactionM= 
    ggplot(df.graph%>% filter(Factor== x), aes(x= sender, y= receiver, 
                                               fill = abs(n)))+
    geom_tile(color ="white", size= 0.5)+
    scale_fill_gradient(low= "white" ,high = "red")+
    facet_grid(~group)+
    coord_equal()+
    theme_cowplot()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          axis.line = element_blank())+
    ggtitle(x)
  p.interactionM
})

pls_n

# subtract 
df.grap.net <- df.graph %>% 
  select(-w, -scaled_weight) %>% 
  pivot_wider(  names_from= group, values_from =n)%>%
  mutate(NF= ifelse(is.na(NF), 0, NF),
         HF= ifelse(is.na(HF), 0, HF))%>%
  mutate(diff= HF-NF)


p.interactionM_net = 
  df.grap.net %>% 
  ggplot(., aes(x= sender, y= receiver, 
                fill = (diff)))+
  geom_tile(color ="black")+
  scale_fill_gradient2(low= "blue" ,mid ="white", high = "red")+
  facet_grid(~Factor)+
  coord_equal()+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(fill = "Difference\ninteraction\nscores")


pls_n[[1]] | p.interactionM_net


df.grap.net
df_gossip= df.graph%>% 
  group_by(sender, Factor, group)%>%
  mutate(sum_sender= sum(n))%>%
  group_by(receiver, Factor, group)%>%
  mutate(sum_receiver= sum(n))%>%
  filter(sender==receiver)%>%
  select(-n)%>%
  rowwise()%>%
  mutate(sum_net = sum_sender- sum_receiver)
df_gossip%>%
  ggplot(.,aes(x=sender, y= sum_net, fill =group))+
  geom_hline(yintercept = 0, col="black")+
  geom_col(position= "dodge", col= "black",width=0.8)+
  facet_grid( ~Factor)+
  theme_cowplot()
  theme(axis.text.x = element_text(angle= 45,hjust= 1), 
        legend.position = "none")+
  labs(x= "", y="Gossip score", fill ="")

