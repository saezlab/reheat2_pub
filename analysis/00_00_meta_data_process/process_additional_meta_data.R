# HEADER --------------------------------------------
#
# Author:     Jan D. Lanzer
# Copyright   Copyright 2024 - Jan D. Lanzer
# Email:      lanzerjan@gmail.com
#
# Date:     2024-01-05
#
# Script Name:    ~/R-projects/reheat2_pilot/process_additional_meta_data.R
#
# Script Description:
# process additional patient level meta data

library(tidyverse)

# HELPER ------------------------------------------------------------------

add_age_range <- function(meta1){
  meta1 <- meta1%>%
    mutate(age_range = cut(age, breaks = seq(0, 90, 10), include.lowest = TRUE, 
                           labels = c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"))
    )
  
}

# REICHART  ---------------------------------------------------------------
#1. we combine the two supplemental files with relevant columns:

meta1= read.csv("data/metadata/raw/Reichart_NIHMS1835045-supplement-Supp__table_1_HF_mod.csv",header = T)%>% 
  as_tibble()
meta2= read.csv("data/metadata/raw/Reichart_NIHMS1835045-supplement-Supp__table_1.csv")%>% 
  as_tibble()

meta1 <- meta1 %>% 
  dplyr::rename(age_range= Age, 
                Donor= Patient, 
                Ethnic.Origin= Ancestry)%>%
  mutate(Age = as.numeric(str_sub(age_range, start = 0, end= 1))*10 +5) # we take the mean of the age range

meta2 <- meta2 %>%
  mutate(Age = ifelse(is.na(Age), as.numeric(str_sub(Age.range , start = 0, end= 1))*10 +5, Age))%>%
  rename(age_range= Age.range)
  # rename(LVEF = Echocardiography..LVEF...)%>%
  # mutate(LVEF = str_replace_all(LVEF, ">", ""))

int= colnames(meta1)[colnames(meta1) %in% colnames(meta2)]

meta.comb<-rbind(meta1[,int], meta2[, int])

#2. we have to change some sample IDs to be mappable to our andata obj.
# from https://github.com/heiniglab/DCM_heart_cell_atlas/blob/main/Additional_note_samples.xlsx
mapping_table <- tibble(
  Donor = c("H2", "H3", "H4", "H5", "H6", "H7"),
  sample_id = c("ED_H26", "ED_H25", "ED_H37", "ED_H15", "ED_H20", "ED_H35")
)

# update names
meta.comb <-  meta.comb  %>%
  mutate(Donor = str_replace_all(Donor, "^(.*)_IC$", "IC_\\1"))%>% # there is an additional change here from suffix to prefix
  left_join(mapping_table, by = c("Donor"))%>%
  mutate(Donor= ifelse(!is.na(sample_id), sample_id, Donor))%>%
  select(-sample_id)


         
# load meta data from andata
coldata= read.csv("data/metadata/Reichart2022_DCM_metadata.csv")%>% 
  as_tibble()

#join both
coldata2= coldata%>% left_join(meta.comb %>% 
                       rename(sample_id = Donor)%>%
                       select(-Sex),
                     by= "sample_id")%>%
  rename(age= Age,
         age_range= age_range,
         race= Ethnic.Origin)%>%
  mutate(age_range = factor(age_range), 
         study= "Reichart2022")

# we will add the information that only affects HF samples as well :
coldata2 <- coldata2 %>% 
  left_join(meta1 %>% 
              rename(sample_id= Donor)%>%
              select(-Ethnic.Origin, -Sex,-age_range,-Age, -BMI, -LVEF,),
            by= "sample_id")

## genotype info
# create a factor for genetic or pv negative:
coldata2 <- coldata2%>% mutate(genetic_HF = ifelse(Primary.Genetic.Diagnosis %in% c("control", "PVneg"), "no", "yes"))
coldata2$genetic_HF

coldata2<-coldata2 %>% mutate(pv = ifelse(Primary.Genetic.Diagnosis %in% c("control", "PVneg"), "no_pv",Primary.Genetic.Diagnosis ))
table(coldata2$pv, coldata2$genetic_HF)
## explant info
table(coldata2$Tissue_Source, coldata2$Sample.site)
# sample accquisition info
coldata2 <- coldata2 %>%
  mutate(sample_accquisition = ifelse(Tissue_Source=="Apical core", "LVAD", 
                                      ifelse(Tissue_Source=="Explanted", "Explanted", NA)))
## sample location
#coldata2 <- coldata2%>% mutate(sample_site = ifelse())
coldata2$Sample.site
coldata2<- coldata2 %>% mutate(sample_site = ifelse(Sample.site=="Apex", "LV_Apex", 
                                         ifelse(Sample.site=="LV", "LV_FreeWall", Sample.site)))
table(coldata2$disease, coldata2$Sample.site)

samples_RV_S= coldata2 %>% filter(disease != "normal" & !Sample.site %in% c("Apex", "LV"))%>% pull(sample_id)

write.csv(samples_RV_S, "data/metadata/Reichart2022_DCM_RVsamples.csv")

coldata2 %>% 
  select(study, sample_id, disease_code, heart_failure, sex, age, age_range, BMI, LVEF,race, everything(),-X)%>%
  write.csv(., "data/metadata/Reichart2022_DCM_metadata_expanded.csv")

# KOENIG ------------------------------------------------------------------
meta1= read.csv("data/metadata/raw/Koenig_44161_2022_28_MOESM3_ESM.csv",header = T)%>% 
  as_tibble()

# meta data from sc file:
coldata= read.csv("data/metadata/Koenig2022_DCM_metadata.csv")%>% as_tibble()

meta1 <- meta1 %>%   rename(age= Age,
                            sex= Sex, 
                            sample_id = Samples,
                            LVEF = Pre.Op.EF,
                            race= Race)

meta1 <- add_age_range(meta1)

coldata2= coldata %>% 
  left_join(meta1%>%
              select(-Condition), by= "sample_id")%>%
  mutate(study= "Koenig2022")%>%
    select(study, sample_id, disease_code, heart_failure, sex, age, age_range,BMI, LVEF,race, everything())

coldata2<- coldata2 %>% 
  mutate(genetic_HF = ifelse(grepl("familial", Etiology.of.HF ), "yes", "no"),
         pv = "no_pv", 
         sample_accquisition= "LVAD_or_Explanted", 
         sample_site ="LV_Apex"
        )

coldata2 %>%
  write.csv(., "data/metadata/Koenig2022_DCM_metadata_expanded.csv")


# CHAFFIN -----------------------------------------------------------------

meta1 <- read_csv(file = "data/metadata/raw/Chaffin_2021-02-03277C-ST1.ClinicalCharacteristics_mod.csv")
Cha <- read.csv("data/metadata/Chaffin2022_DCM_metadata.csv")

meta1 <-meta1%>% 
  mutate(sample_id = paste0("P", Individual))%>%
  rename(BMI = `Body Mass\nIndex\n(kg/m2)`,
         LVEF = `Left Ventricular\nEjection Fraction\n(%)`)%>%
  mutate(genetic_HF = ifelse((`Cardiomyopathy gene mutation` == "--" & Disease!= "NF") | is.na(`Cardiomyopathy gene mutation`) , "no", "yes"))%>%
  mutate(sample_site ="LV_FreeWall", 
         sample_accquisition= "Explanted")%>%
  select(sample_id, LVEF, BMI, genetic_HF, everything())

table(meta1$`Cardiomyopathy gene mutation`, meta1$Disease)

## PV
meta1$pv <- str_extract(meta1$`Cardiomyopathy gene mutation`, "^[^_]+")
meta1 <- meta1 %>% mutate(pv = ifelse(pv =="--" | is.na(pv), "no_pv", pv))
table(meta1$pv, meta1$genetic_HF)
coldata <- Cha %>% left_join(meta1, by= "sample_id")

coldata %>%
  write.csv(., "data/metadata/Chaffin2022_DCM_metadata_expanded.csv")

# SIMONSON ----------------------------------------------------------------
# there is no additional meta data on sample level 
coldata= read.csv("data/metadata/Simonson2023_ICM_metadata.csv")%>% as_tibble()

coldata <- coldata %>% 
  mutate(pv = "no_pv", 
         sample_accquisition= "Explanted", 
         sample_site ="LV_FreeWall",
         genetic_HF= "no"
)
coldata%>% 
  write.csv(., "data/metadata/Simonson2023_ICM_metadata_expanded.csv")
