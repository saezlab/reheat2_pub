# HEADER --------------------------------------------
#
# Author:     Jan D. Lanzer
# Copyright   Copyright 2024 - Jan D. Lanzer
# Email:      lanzerjan@gmail.com
#
# Date:     2024-01-05
#
# Script Name:    
#
# Script Description:
# gather joint meta data

library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

col_list<- readRDS("color_list_figures.rds")

Koe <- read.csv("data/metadata/Koenig2022_DCM_metadata_expanded.csv")
Rei <- read.csv("data/metadata/Reichart2022_DCM_metadata_expanded.csv")
Sim <- read.csv("data/metadata/Simonson2023_ICM_metadata_expanded.csv")
Cha <- read.csv("data/metadata/Chaffin2022_DCM_metadata_expanded.csv")

meta_data= list()
# helper ------------------------------------------------------------------

add_age_range <- function(meta1){
  meta1 <- meta1%>%
    mutate(age.range = cut(age, breaks = seq(0, 90, 10), include.lowest = TRUE, 
                           labels = c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"))
    )
  
}

# Function to calculate the mean value of a range
calculate_mean_from_range <- function(x) {
  if (grepl("-", x)) {
    range_values <- as.numeric(strsplit(x, "-")[[1]])
    return(mean(range_values, na.rm = TRUE))
  } else {
    return(as.numeric(x))
  }
}


# -------------------------------------------------------------------------


Sim <- add_age_range(Sim) %>% mutate(study = "Simonson2023")
Cha <- add_age_range(Cha) %>% mutate(study = "Chaffin2022")
Koe<- add_age_range(Koe) %>% mutate(study = "Koenig2022")
Rei <- add_age_range(Rei) %>% mutate(study = "Reichart2022")

df.list <- list( Cha, Koe, Sim, Rei)

#fill in missing columns
col.names <- c("BMI", "LVEF", "race")

res= lapply(df.list, function(x){
  print(unique(x$study))
  for (i in col.names){
    print(i)
    if(!i %in% colnames(x)){
      x[i] <- NA
    }
  }
  return(x)
})

lapply(res, colnames)

df <- lapply(res, function(x){
  print(unique(x$study))
  x %>% select(study, sample_id, disease_code, heart_failure, sex, age, 
               age.range, LVEF, BMI, race, sample_site, sample_accquisition,
               genetic_HF, pv)
})%>% do.call(rbind,. )

df <-df %>%
  as_tibble()%>% 
  mutate(age.range = factor(age.range))


# clean up joint meta data ------------------------------------------------------

# race
unique(df$race)
df <- df %>% 
  mutate(race = ifelse(race == "W", "Caucasian", race),
         race = ifelse(race == "AA", "African American", race))

#lvef
unique(df$LVEF)
df <- df %>% 
  mutate(LVEF = ifelse(LVEF %in% c("UNK", "No echo", "N/A"), NA, LVEF),
         LVEF = as.numeric(LVEF))
df

#BMI
unique(df$BMI)

# some BMIs are ranged, we transform to numeric by also calculating the mean
df <- df %>%
  mutate(BMI = sapply(BMI, calculate_mean_from_range))%>%
  rename(age_range = age.range)


#save 
write.csv(df, "data/metadata/combined_meta.csv")

df<- read.csv("data/metadata/combined_meta.csv")

#remove the RV and septal samples form reichart<
samples_RV_S<-read.csv("data/metadata/Reichart2022_DCM_RVsamples.csv")
df <- df %>% filter(!sample_id %in% samples_RV_S$x)

meta_data$target_sc <- df

# calc_stats--------------------------------------------------------------------

genetic_hf <- df %>%
  filter(heart_failure=="HF")%>%
  group_by( genetic_HF, disease_code)%>% 
  count()%>% group_by(disease_code) %>%mutate(prop= n/sum(n))
genetic_hf
genetic_hf%>%
  ggplot(aes(fill= genetic_HF, y =n, x= disease_code))+
  geom_col(col= "black", position= position_dodge2(width = 0.3, preserve = "single"),
           width= 0.6)+
  scale_fill_manual(values= names(etio_vec))+
  theme_cowplot()+
  labs(fill =  "genetic variant\nor familial HF",x= "")

genetic_hf%>%
  ggplot(aes(fill= genetic_HF,y= prop,   x= disease_code))+
  geom_col(col= "black"  , width= 0.8)+
  scale_fill_manual(values= names(etio_vec))+
  theme_cowplot()+
  labs(fill =  "genetic variant\nor familial HF",x= "")
genetic_hf
  
df %>% filter(is.na(sample_site))


df %>%
  group_by(heart_failure)%>%
  summarise(m.age= median(age, na.rm = T),
            iqr.age = IQR(age, na.rm = T),
            prop_female = mean(sex == "female", na.rm = TRUE)) 


# plot distributions per study:  ----------------------------------------------------

p1.age = df%>% 
  ggplot(., aes(x= study,age, fill = heart_failure))+
  geom_jitter(col = "darkgrey")+
  #geom_violin(alpha= 0.6)+
  geom_boxplot(alpha= 0.5, width = 0.9)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 45, hjust= 1))+
  labs(y= "age (years)",
       x= "")
p1.age

p1.lvef = df%>% 
  ggplot(., aes(x= study,LVEF, fill = heart_failure))+
  geom_jitter(col = "darkgrey")+
  #geom_violin(alpha= 0.6)+
  geom_boxplot(alpha= 0.5, width = 0.9)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 45, hjust= 1))+
  labs(y= "LVEF (%)",
       x= "")

p1.bmi = df%>% 
  ggplot(., aes(x= study,BMI, fill = heart_failure))+
  geom_jitter(col = "darkgrey")+
  #geom_violin(alpha= 0.6)+
  geom_boxplot(alpha= 0.5, width = 0.9)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 45, hjust= 1))+
  labs(y= "BMI (kg/m²)",
       x= "")
legend <- get_legend(
  # create some space to the left of the legend
  p1.age + theme(legend.box.margin = margin(0, 0, 0, 12))
)
p.cont <- plot_grid(p1.age+theme(legend.position = "none"),
          p1.bmi+theme(legend.position = "none"),
          p1.lvef+theme(legend.position = "none"),
          legend)


p11 = df%>% 
  ggplot(., aes(x= study,age, fill = sex))+
  geom_jitter(col = "darkgrey")+
  #geom_violin(alpha= 0.6)+
  geom_boxplot(alpha= 0.5, width = 0.9)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 45, hjust= 1))+
  labs(y= "age (years)",
       x= "")
p11
p2= ggplot(df, aes(x = age.range, fill = study)) +
  geom_bar(position = "dodge", color = "black") +
  labs(title = "Distribution of age.range per Study",
       x = "Age Range",
       y = "Count") +
  theme_cowplot()

p3= ggplot(df, aes(x = sex, fill = age.range)) +
  geom_bar(position = "dodge", color = "black") +
  facet_wrap(~ study) +
  labs(title = "Distribution of age.range per Study",
       x = "Age Range",
       y = "Count") +
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 45, hjust= 1))
  


p4= ggplot(df, aes(x = sex, fill = age.range)) +
  geom_bar(position = "dodge", color = "black") +
  facet_wrap(~ study) +
  labs(title = "Distribution of age.range per Study",
       x = "Age Range",
       y = "Count") +
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 45, hjust= 1))


p5= ggplot(df, aes(x = study, fill = heart_failure)) +
  geom_bar(position = "stack", color = "black") +
  labs(title = "",
       x = "",
       y = "Count") +
  theme_cowplot()+
  coord_flip()
p.eti= ggplot(df, aes(x = study, fill = disease_code)) +
  geom_bar(position = "stack", color = "black") +
  labs(title = "",
       x = "",
       y = "Count") +
  theme_cowplot()+
  coord_flip()

p.race = df %>%
  mutate(race = factor(race, levels=c("Caucasian", "African American", "South Asian")))%>% 
  ggplot( aes(x = study, fill = race)) +
    geom_bar(position = "fill",color="black")+
    labs(title = "",
         x = "",
         y = "Proportion") +
    theme_cowplot()+
    coord_flip()

p.sex = df %>% 
  ggplot( aes(x = study, fill = sex)) +
  geom_bar(position = "fill",color="black")+
  labs(title = "",
       x = "",
       y = "Proportion") +
  theme_cowplot()+
  coord_flip()

p.cat <- plot_grid(p.race, p.sex, rel_widths = c(1 , 0.8))
p.cat

p.sample_site = df %>% 
  ggplot( aes(x = study, fill = sample_site)) +
  geom_bar(position = "fill",color="black")+
  labs(title = "",
       x = "",
       y = "Proportion") +
  theme_cowplot()+
  facet_grid(~heart_failure)+
  coord_flip()
p.sample_site
# plot distributions total.  ----------------------------------------------

p6= ggplot(df, aes(x = sex, fill = age.range)) +
  geom_bar(position = "dodge", color = "black") +
  labs(title = "Distribution of age.range per Study",
       x = "Age Range",
       y = "Count") +
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 45, hjust= 1))

p7= df%>% 
  ggplot(., aes(x= sex,age))+
  geom_jitter(col = "darkgrey")+
  geom_violin(alpha= 0.6)+
  geom_boxplot(alpha= 0.5, width = 0.2)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 45, hjust= 1))+
  labs(y= "age (years)",
       x= "")

p8= df%>% 
  ggplot(., aes(x= heart_failure, age))+
  geom_jitter(col = "darkgrey")+
  geom_violin(alpha= 0.6)+
  geom_boxplot(alpha= 0.5, width = 0.2)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 45, hjust= 1))+
  labs(y= "age (years)",
       x= "")


p9= df%>% 
  ggplot(., aes(x=  heart_failure, fill = sex))+
  geom_bar(position = "stack", col = "black") +
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 45, hjust= 1))+
  labs(y= "count",
       x= "")
p91= df%>% 
  ggplot(., aes(x=  heart_failure, fill = age.range))+
  geom_bar(position = "stack", col = "black") +
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 45, hjust= 1))+
  labs(y= "count",
       x= "")
p91

p92= df%>% 
  mutate(disease_code = ifelse(study=="Reichart2022" & disease_code== "DCM",
                               "DCM_genetic", 
                               disease_code))%>%
  ggplot(., aes(x=  sex, fill = disease_code))+
  geom_bar(position = "stack", col = "black") +
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 45, hjust= 1))+
  labs(y= "count",
       x= "")
p92

p10= df%>% 
  ggplot(., aes(x=  heart_failure, fill = study))+
  geom_bar(position = "stack", col = "black") +
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 45, hjust= 1))+
  labs(y= "count",
       x= "")



pdf("output/figures/meta_infos.pdf")
p1
p2
p3
p4
p5
p6
p7
p8
p9
p10
dev.off()




# plot NAs ----------------------------------------------------------------

columns_to_analyze <- c("sex",  "age" ,"LVEF" ,"BMI" ,  "race" , "sample_site", 
                        "sample_accquisition", "genetic_HF")

na_counts <- df %>%
  group_by(study) %>%
  summarise(across(all_of(columns_to_analyze), ~sum(is.na(.))))

sample_counts <- df %>%
  group_by(study) %>%
  summarise(x= length(unique(sample_id)))

# Reshape the data for plotting
na_counts_long <- na_counts %>%
  pivot_longer(-study, names_to = "Column", values_to = "NA_Count")

# Calculate the proportion of NAs
na_counts_long <- na_counts_long %>%
  left_join(sample_counts)%>%
  mutate(Proportion = NA_Count / x)

# Plot the stacked bar plot
ggplot(na_counts_long, aes(x = reorder(Column,Proportion), y = Proportion, fill = factor(Column))) +
  geom_bar(stat = "identity") +
  labs(title = "Proportion of NAs",
       x = "",
       y = "Proportion of NAs") +
  scale_fill_discrete(name = "Columns") +
  theme_minimal()+
  theme(legend.position ="none")
  
ggplot(na_counts_long, aes(y = study, fill = Proportion, 
                           x = factor(Column), 
                           label= round(Proportion, 1))) +
  geom_tile(color="black") +
  geom_text()+
  labs(title = "",
       x = "",
       y = "",
       fill = "Proportion of NAs") +
  scale_fill_gradient(low="red", high ="white")+
 # scale_fill_discrete(name = "Columns") +
  theme_cowplot()+
  theme(axis.text.x=element_text(angle=45,hjust= 1))

ggplot(df, aes(x=study))

# compare stats between etiologies ----------------------------------------

df %>% ggplot(.,aes(x= disease_code, y= LVEF, fill = study))+
    geom_boxplot()
  
  
  table(df$disease_code)
  
  

# add bulk meta data ------------------------------------------------------
read.csv("../HF_meta-analysis/data/clinical_description/Tables for MetaHeart Manuscript - Table 2 - Clinical Characteristics.csv")
read.csv("../HF_meta-analysis/data/clinical_description/Tables for MetaHeart Manuscript - Table 3 - Tissue Preparation.csv")
  
  
## ReHeaT 1 data->
load("../HF_meta-analysis/data/clinical_description/ClinicalCharacteristics_age.ro")
load("../HF_meta-analysis/data/clinical_description/ClinicalCharacteristics_gender.ro")
load("../HF_meta-analysis/data/clinical_description/ClinicalCharacteristics_ef.ro")

gender<- gender%>% rownames_to_column("study")%>% 
  pivot_longer(cols= c(female_HF, female_NF), names_to = "heart_failure", values_to = "proportion_female")%>%
  mutate(heart_failure =ifelse(grepl("HF", heart_failure), "HF", "NF"))

ages<- ages %>% rownames_to_column("study")%>%
  pivot_longer(cols= c(HFMeanAge, NFMeanAge ), names_to = "heart_failure", values_to = "age_mean")%>%
  mutate(heart_failure =ifelse(grepl("HF", heart_failure), "HF", "NF"))%>%
  pivot_longer(cols= c(HFSDAge, NFSDAge), names_to = "heart_failure2", values_to = "age_sd")%>%
    mutate(heart_failure2 =ifelse(grepl("HF", heart_failure2), "HF", "NF"))%>%
  filter(heart_failure==heart_failure2)%>%
  select(-heart_failure2)
ages 
ef<- ef%>% rownames_to_column("study")%>%
  pivot_longer(cols= c(HFMeanEF, NFMeanEF ), names_to = "heart_failure", values_to = "lvef_mean")%>%
  mutate(heart_failure =ifelse(grepl("HF", heart_failure), "HF", "NF"))%>%
  pivot_longer(cols= c(HFSdEF, NFSdEF), names_to = "heart_failure2", values_to = "lvef_sd")%>%
  mutate(heart_failure2 =ifelse(grepl("HF", heart_failure2), "HF", "NF"))%>%
  filter(heart_failure==heart_failure2)%>%
  select(-heart_failure2)

#merge to one data frame
bulk_df <- ages %>% 
  left_join(ef, by= c("study", "heart_failure"))%>%
  left_join(gender, by= c("study", "heart_failure"))%>%
  mutate(modality="bulk")

## add the new studies (bulk):

METAheart<- readRDS(file = "output/reheat1/METAheart2023.rds")

#only flam19 and hua19 has additional meta data->
flam19_df<- METAheart$Flam19$TARGETS%>% 
  mutate(heart_failure= ifelse(HeartFailure== "yes", "HF", "NF"))%>%
  group_by(heart_failure)%>%
  summarize(age_mean= mean(Age,na.rm= T), 
            age_sd = sd(Age,na.rm= T), 
            proportion_female = mean(Gender == "female", na.rm= T)*100)%>%
  mutate(study= "Flam19", lvef_mean= NA, lvef_sd= NA, modality = "bulk")%>%
  select(colnames(bulk_df))

hua19_df<- METAheart$Hua19$TARGETS%>% 
  mutate(heart_failure= ifelse(HeartFailure== "yes", "HF", "NF"))%>%
  group_by(heart_failure)%>%
  summarize(age_mean= mean(Age,na.rm= T), 
            age_sd = sd(Age,na.rm= T), 
            proportion_female = mean(Gender == "female", na.rm= T)*100)%>%
  mutate(study= "Hua19", lvef_mean= NA, lvef_sd= NA, modality = "bulk")%>%
  select(colnames(bulk_df))

bulk_df <- rbind(bulk_df, flam19_df, hua19_df)  

#add empty data for the other new studies
bulk_df<- rbind(bulk_df, 
      map(c("Forte22", "Wang22", "Rao21"), function(x){
        flam19_df %>%mutate(study= x)%>%
          mutate(across(-c(study,modality,  heart_failure), ~ NA)
                 )
        })%>% do.call(rbind, .)
      )


### get sc ->
sc_df<- df%>% 
  group_by(study, heart_failure)%>%
  summarize(age_mean= mean(age,na.rm= T), 
            age_sd = sd(age,na.rm= T), 
            proportion_female = mean(sex == "female", na.rm= T)*100, 
            lvef_mean = mean(LVEF, na.rm= T),
            lvef_sd= sd(LVEF,na.rm= T))%>%
  mutate(modality= "sc")

# merge both-> 
dfcomb<- rbind(bulk_df, sc_df)
dfcomb<- dfcomb%>% 
  rowwise()%>%
  mutate(year=as.numeric(str_extract_all(study, "\\d+")[[1]]), 
         year= ifelse(year<100, year+2000, year))
##

meta_data$demo <- dfcomb

dfcomb %>% 
  ggplot(aes(x=heart_failure, y= age_mean, fill = modality))+
  geom_boxplot()

x= dfcomb%>%filter(modality=="sc")%>% pull(study)%>% unique()
y= dfcomb%>%filter(modality=="bulk")%>% pull(study)%>% unique()
dfcomb %>% 
  arrange(modality)%>%
  mutate(study= factor(study, levels= c( y,x)))%>%
  ggplot(aes(x=study,  y= proportion_female, fill = heart_failure))+
  geom_col(position = "dodge", width = 0.7, color="black")+
  coord_flip()+
  geom_hline(yintercept = 50)+
  theme_cowplot()+
  scale_fill_manual(values= c("black", "grey"))+
  labs(x="", y= "female sex (%)")

dfcomb %>% 
  arrange(modality)%>%
  mutate(study= factor(study, levels= c( y,x)))%>%
  ggplot(aes(x=study,  y= proportion_female, fill = heart_failure))+
  geom_col(position = "dodge", width = 0.7, color="black")+
  coord_flip()+
  geom_hline(yintercept = 50,linetype=3)+
  theme_cowplot()+
  scale_fill_manual(values= c("black", "grey"))+
  labs(x="", y= "female sex (%)")


p.age
p.sex<- dfcomb%>%
  mutate(label= ifelse(heart_failure=="NF","",study ))%>%
  ggplot(aes(x=year, y=age_mean,  color=heart_failure))+
  geom_hline(yintercept = mean(dfcomb%>% filter(heart_failure=="HF")%>% pull(age_mean), na.rm = T),
             linetype=3,linewidth=0.8, color="black")+
  geom_hline(yintercept = mean(dfcomb%>% filter(heart_failure!="HF")%>% pull(age_mean), na.rm = T), 
             linetype=3, linewidth=0.8, color="darkgreen")+
  #geom_errorbar(aes(ymin=age_mean-age_sd, ymax=age_mean+age_sd),alpha= .4, 
   #             linetype=1,
    #            linewidth=.4)+
  ggrepel::geom_label_repel(aes(label=label, fill= modality),color="white",
                            max.overlaps = 2000,
                            show.legend = T, 
                            alpha= 0.8  )+
  geom_smooth(formula = y~x, na.rm = T, alpha= 0.2, se= F)+
  scale_color_manual(values= c("#E83151", "#757780"))+
  scale_fill_manual(values= c("darkblue", "black"))+
  #geom_point()+
  theme_cowplot()+
  ylim(c(25,80))+
  xlim(c(2005,2023))+
  labs(y="Age\n(study mean ± sd)", x="Year",  fill="Modality",color= "Cohort")+
  theme(panel.background = element_rect(fill="lavender"))


pdf("output/figures/meta_data_presentation/age_sex_plot.pdf", 
    height= 5, width= 7)

  plot_grid(p.age+theme(legend.position = "none"), p.sex, ncol = 1, align= "v")

dev.off()

dfcomb%>% group_by(heart_failure, modality)%>%
  summarise(m.age= mean(age_mean, na.rm = T),
            iqr.age = IQR(age_mean, na.rm = T),
            prop_female = mean(proportion_female , na.rm = TRUE)) 
dfcomb%>% group_by( heart_failure)%>%
  summarise(m.age= mean(age_mean, na.rm = T),
            iqr.age = IQR(age_mean, na.rm = T),
            prop_female = mean(proportion_female , na.rm = TRUE)) 



METAheart = readRDS(file = "output/reheat1/METAheart2023.rds") #main object

experiments = names(METAheart)
names(experiments) = experiments

new_study_ids= c("Forte22", "Flam19", "Wang22", "Rao21", "Hua19")
study_order <- c(old_study_ids, new_study_ids)
old_study_ids <- as.character(experiments[!experiments %in% new_study_ids])
a <- c(rep("black", length(old_study_ids)), rep("darkred", length(new_study_ids)))

# For labeling
experiment_size = sort(unlist(lapply(METAheart,
                                     function(x) ncol(x$GEX))),
                       decreasing = T)

plot(experiment_size)
sum(experiment_size)

# create a sample size overview:

df_blk= lapply(names(METAheart), function(x){
  METAheart[[x]]$TARGETS%>% dplyr::select(HeartFailure)%>% mutate(study= x)
})%>% do.call(rbind, .)

df_blk <- df_blk %>% group_by(study, HeartFailure)%>% summarise(x= n())%>% 
  mutate(new= ifelse(study %in% new_study_ids, "y", "n"),
         HeartFailure= factor(HeartFailure, levels= c("yes", "no")))


df_sc <- df%>% group_by(study, heart_failure)%>% count

df_blk <- df_blk%>% rename(heart_failure= HeartFailure, n=x)

df_comb_counts<- rbind(df_sc %>% mutate(new= "y"), 
                       df_blk %>% 
                         mutate(heart_failure =ifelse(grepl("yes", heart_failure), "HF", "NF")
                                )
                       )
df_comb_counts$study<-str_replace_all(df_comb_counts$study,pattern = "\\.",replacement =  "-")


dfcomb2<- dfcomb %>% left_join(df_comb_counts%>% ungroup(), by= c("study", "heart_failure"))%>% 
  ungroup()
meta_data$sample_size<- dfcomb2

pdf("output/figures/meta_data_presentation/sample_size_overview.pdf" , 
    width = 10, height= 4.7)
p.samplesize
dev.off()

#p.samplesize2 <-
  dfcomb2%>%
  ggplot(., aes(x= year ,y= n))+
    geom_point(aes(colour=heart_failure), size=12) + 
    geom_point(shape = 1,size = 12,colour = "black")
  geom_point(aes(fill= heart_failure, color= new), size= 4)
  scale_color_manual(values= rev(c("black", "white")))+
  scale_fill_manual(values = c("#eb4034", "#0b4f32"))+
  coord_flip()+
  labs(y= "", x="")+
  theme_cowplot()


dfcomb%>% 
  group_by(modality, heart_failure)%>%
  summarise(female= mean(proportion_female, na.rm=T))%>%
  mutate(modal= paste0(modality, "_", heart_failure),
         male= 100-female)%>%
  pivot_longer(cols= c(male, female), names_to = "sex", values_to = "percent")%>%
  ggplot(aes(x=modal, y= percent, fill = sex))+
  geom_col(position = "fill")

dfcomb%>%
  ggplot(., aes(x= study, y= ))
  

# compare sampling location and sample accquistion ------------------------

bulk_reheat1_tissueinfo<-read.csv("../HF_meta-analysis/data/clinical_description/Tables for MetaHeart Manuscript - Table 3 - Tissue Preparation.csv")
tissueinfo1<- bulk_reheat1_tissueinfo%>% select(Study.ID, Ventricle., Location, Location.2, reason.of.biopsy...diseased)

bulk_location<- tissueinfo1%>% 
  rename(sample_accquisition = reason.of.biopsy...diseased, 
         sample_site= Location)%>%
  mutate(sample_site = ifelse(grepl("free wall", sample_site), "LV_FreeWall", 
                                      ifelse(grepl("apex", sample_site), "LV_Apex", NA)))%>%
  mutate(study= ifelse(Study.ID == "Tarazón14", "Tarazon14", Study.ID ),
         modality= "bulk")

df_location<- rbind(bulk_location%>%full_join(dfcomb2)%>%
                      filter(modality=="bulk")%>%
  select(study, heart_failure, sample_site, n, modality),
  df%>% group_by(study, heart_failure, sample_site)%>%
  count()%>% mutate(                  modality= "sc")
  )

# for new studies we add location based on methods sections
#flam19 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9910091/)
# all explanted+ freewall
#wang22 (https://www.frontiersin.org/journals/cardiovascular-medicine/articles/10.3389/fcvm.2022.919355/full)
#explanted free wall
# forte21
# No info
# rao21(https://link.springer.com/article/10.1007/s00395-021-00897-1#Sec2)
# explanted and but no location info

df_location<- df_location%>% mutate(sample_site = ifelse(study=="Flam19",
                                           "LV_FreeWall", 
                                           sample_site), 
                      sample_site = ifelse(study=="Wang22",
                                           "LV_FreeWall", 
                                           sample_site)
                      
                      )


meta_data$location= df_location


# sampole_accquisition ----------------------------------------------------


bulk_sample_acc_df<- lapply(names(METAheart), function(x){
  if("HTx" %in% colnames(METAheart[[x]]$TARGETS)){
    METAheart[[x]]$TARGETS%>% dplyr::select(HTx, HeartFailure)%>% mutate(study= x)
  }
})%>% do.call(rbind, .)%>% group_by(study, HTx, HeartFailure)%>% summarise(n= n())%>% 
  rename(sample_accquisition= HTx,
         heart_failure =HeartFailure)%>%
  mutate(sample_accquisition = ifelse(sample_accquisition=="yes", "Explanted", "LVAD"),
         heart_failure = ifelse(heart_failure=="yes", "HF", "NF"),
         modality= "bulk")%>% ungroup()

add_row(bulk_sample_acc_df, study = "Flam19", sample_accquisition= "Explanted", modality="bulk")
add_row(bulk_sample_acc_df, study = "Wang22", sample_accquisition= "Explanted", modality="bulk")
add_row(bulk_sample_acc_df, study = "Rao21", sample_accquisition= "Explanted", modality="bulk")

sc_sample_acc_df <- df%>% group_by(study, sample_accquisition, heart_failure)%>% summarise(n= n())%>%
  mutate(sample_accquisition= ifelse(sample_accquisition=="LVAD_or_Explanted", NA, sample_accquisition),
         modality= "sc")
df_acc<- rbind(bulk_sample_acc_df, sc_sample_acc_df)

meta_data$accquisition <- df_acc
# compare etiologies ------------------------------------------------------

#reheat1
df_blk_etiologies <- lapply(names(METAheart), function(x){
  print(x)
  if(!"DCM" %in% colnames(METAheart[[x]]$TARGETS)){
    METAheart[[x]]$TARGETS<-METAheart[[x]]$TARGETS%>% mutate(DCM= "no")
  }
  METAheart[[x]]$TARGETS%>% dplyr::select(HeartFailure, DCM)%>% mutate(study= x)
})%>% do.call(rbind, .)

df_blk_etiologies<- df_blk_etiologies%>% 
  mutate(
    disease_code = ifelse(DCM =="yes" & HeartFailure == "yes", "DCM","ICM"),
    disease_code = ifelse(HeartFailure == "no", "NF",disease_code),
    disease_code = ifelse(study=="Flam19" & DCM == "no", "HCM", disease_code),
    disease_code = ifelse(study=="Rao21", NA, disease_code))%>%
  rename(heart_failure = HeartFailure)%>%
  dplyr::select(study, heart_failure, disease_code)%>%
  mutate(modality= "bulk")
df_sc_etiologies<- df%>% dplyr::select(study, heart_failure, disease_code)%>%
  mutate(modality= "sc")
df_etio <- rbind(df_blk_etiologies, df_sc_etiologies)

df_etio
meta_data$etiologiess <-df_etio

## waffle plots: 
p.etio <- df_etio %>% 
  group_by(disease_code, modality)%>%
  count()%>%
  ggplot(data = ., 
         aes(fill=disease_code, values=n/5)) +
  geom_waffle(color = "white")+
  scale_fill_manual(values = c("#387780", "#E83151", "#D2CCA1", "#757780"))+
  facet_grid(cols= vars(modality),  space="free", scales="free")+
  theme_void()

p.etio <- df_etio %>% 
  group_by(disease_code, modality, study)%>%
  count()%>%
  ggplot(data = ., 
         aes(fill=disease_code, values=n/5)) +
  geom_waffle(color = "whdite")+
  scale_fill_manual(values = c("#387780", "#E83151", "#D2CCA1", "#757780"))+
  facet_grid(cols= vars(modality),  space="free", scales="free")+
  theme_void()


# save meta_data ----------------------------------------------------------
saveRDS(meta_data, "output/meta_data_for_plotting.rds")
