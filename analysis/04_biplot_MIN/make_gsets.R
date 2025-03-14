# Copyright (c) [2024] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we calculate an overall annotation of ReHeaT
#'

library(tidyverse)
library(qusage)

# Gene set database

#Cytokines - Myeloid/Lymphoid
cytosig <- read_tsv("./data/gsets/cytosig_centroid.csv") %>%
  pivot_longer(-gene, names_to = "source", values_to = "weight") %>%
  dplyr::rename("target" = gene)

write_csv(cytosig, "./data/gsets/pathophysiol_processed/immune.csv")

#Cytokine dictionary
#cytoAtlas <- read_tsv('./data/gsets/Atlascyt__deresults.tsv') %>%
#  pivot_longer(-ID) %>%
#  separate(name, into = c('colid', 'filtering', 'de', 'contrast', 'dataset', 'study'), sep = '__') %>%
#  pivot_wider(names_from = colid, values_from = value)

#NABA - Fib/PC/Endo/vSMCs + metalloproteinases
NABA_gsets <- qusage::read.gmt("./data/gsets/MSDB_NABA.gmt") %>%
  enframe(name = "source", value = "target") %>%
  unnest() %>%
  mutate(weight = 1)

MMPs_gsets <- qusage::read.gmt("./data/gsets/MSDB_MMPs.gmt") %>%
  enframe(name = "source", value = "target") %>%
  unnest() %>%
  mutate(weight = 1)

mural_gsets <- bind_rows(NABA_gsets, MMPs_gsets) %>%
  unique()

write_csv(mural_gsets, "./data/gsets/pathophysiol_processed/mural.csv")

#Cardiomyocytes
# Calcium channels
# adherens, desmosomes and gap junctions
# (gap AND junction) OR desmosome OR adheren OR sarcomere OR (sarcotubular AND system) OR sarcolemma OR (sarcoplasmic AND reticulum) OR (hypertrophy)
# (CALCIUM AND SIGNALING) OR (ION AND CHANNEL) OR (calcium AND cycling)
# APOPTOSIS OR NECROSIS
# (Lipotoxicity) OR (Mitochondrial AND dysfunction) OR (Fatty AND acid AND METABOLISM) OR (Glucose AND METABOLISM)
# (Mitochondrial AND dysfunction)
# OR (adrenergic AND receptor) (Renin AND ANGIOTENSIN AND aldosterone)

CM_gsets <- qusage::read.gmt("./data/gsets/MSDB_CM.gmt") %>%
  enframe(name = "source", value = "target") %>%
  unnest() %>%
  mutate(weight = 1) %>%
  unique()

write_csv(CM_gsets, "./data/gsets/pathophysiol_processed/cm.csv")

# Vasculature
# Endothelin OR (Nitric AND oxide) OR Angiogenesis OR (vessel AND proliferation) OR vasculature OR angiogenesis OR (blood AND vessel)

vasc_gsets <- qusage::read.gmt("./data/gsets/MSDB_vasc.gmt") %>%
  enframe(name = "source", value = "target") %>%
  unnest() %>%
  mutate(weight = 1) %>%
  unique()

write_csv(vasc_gsets, "./data/gsets/pathophysiol_processed/vasc.csv")
