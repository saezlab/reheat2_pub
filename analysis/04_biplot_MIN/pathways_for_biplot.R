# HEADER --------------------------------------------
#
# Author:     Jan D. Lanzer
# Copyright   Copyright 2024 - Jan D. Lanzer
# Email:      lanzerjan@gmail.com
#
# Date:     2024-07-18
#
# Script Name:    
#
# Script Description:
# TRY TO 

library(tidyverse)
X<-read_csv("data/multicell_sets_biplot.csv")
Y<-read_csv("data/singlecell_sets_biplot.csv")

multicell_pways= c(## fact1 & fact2
                    "BIOCARTA_VIP_PATHWAY", 
                   "BIOCARTA_HDAC_PATHWAY", 
                   ## fact1 
                   "BIOCARTA_CALCINEURIN_PATHWAY", 
                   "GOBP_POSITIVE_REGULATION_OF_NITRIC_OXIDE_METABOLIC_PROCESS"
                   "HP_VENTRICULAR_HYPERTROPHY"
                   "HALLMARK_FATTY_ACID_METABOLISM"
                   
                   ## fact2
                   "HALLMARK_HYPOXIA"
                   "BIOCARTA_IGF1MTOR_PATHWAY"
                   "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
                   
                   )


singlecell_pways= c(##CM
                   "HP_ABNORMAL_SARCOMERE_MORPHOLOGY", 
                   "HP_ABNORMAL_Z_DISC_MORPHOLOGY",
                   "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                   "KEGG_FATTY_ACID_METABOLISM", 
                   
                   "KEGG_CARDIAC_MUSCLE_CONTRACTION",
                   "KEGG_CALCIUM_SIGNALING_PATHWAY",
                   
                   ##FIB
                   "NABA_COLLAGENS",
                   "BIOCARTA_RECK_PATHWAY",
                   
                   ## ENDO 
                   "GOBP_BLOOD_VESSEL_ENDOTHELIAL_CELL_DIFFERENTIATION",
                   "GOBP_LYMPHANGIOGENESIS", 
                   "HALLMARK_INTERFERON_GAMMA_RESPONSE", 
                   "BIOCARTA_NO1_PATHWAY", 
                   "WP_ANGIOGENESIS",
                   
                   ##PC
                   "GOBP_RESPONSE_TO_NITROSATIVE_STRESS",
                   
                   
                   ##Myeloid
                   "HALLMARK_IL6_JAK_STAT3_SIGNALING"
                   
                   )

X %>% 
  filter(sign(F1_coord)!= sign(F2_coord))%>%
  arrange(desc(F1_coord))%>% print(n=100)

X %>% 
  filter(sign(F1_coord)== sign(F2_coord))%>%
  arrange(desc(F1_coord))%>% print(n=100)
Y$ct %>% unique()
Y %>% 
  filter(sign(F1_coord)== sign(F2_coord),
         ct =="Myeloid")%>%
  arrange(desc(F1_coord))%>% print(n=100)

Y %>% 
  filter(sign(F1_coord)!= sign(F2_coord),
         ct =="Myeloid")%>%
  arrange(desc(F1_coord))%>% print(n=100)


