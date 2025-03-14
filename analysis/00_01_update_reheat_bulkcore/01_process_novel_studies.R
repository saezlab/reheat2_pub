# HEADER --------------------------------------------
#
# Author:     Jan D. Lanzer
# Copyright   Copyright 2023 - Jan D. Lanzer
# Email:      lanzerjan@gmail.com
#
# Date:     2023-09-21
#
# Script Name:    ~/R-projects/reheat2_pilot/update_reheat/process_magnet.R
#
# Script Description:
# We will load the meta heart object and expand it with new data sets. 

library(tidyverse)
library(biomaRt)
library(limma)
library(edgeR)
library(annotables)


METAheart= readRDS("~/R-projects/HF_meta-analysis/data/METAheart.rds")
METAcounts <- readRDS("output/reheat1/Metaheart_counts2023.rds")

# add magnet --------------------------------------------------------------


## load MAGENET

gc= read.csv("~/R-projects/Collaborations/cheerio/raw_data/magnet/gene_count_matrix.csv")
#tc= read.csv("data/raw_data/magnet/transcript_count_matrix.csv")
meta= read.csv("~/R-projects/Collaborations/cheerio/raw_data/magnet/MAGE_metadata.txt")
pheno= read.csv("~/R-projects/Collaborations/cheerio/raw_data/magnet/phenoData.csv")

## translate gene IDs
gex.df= as_tibble(gc)%>% 
  left_join(annotables::grch38 %>% 
              rename(gene_id = ensgene,
                     gene= symbol )%>%
              dplyr::select(gene_id, gene))#%>% 
## sum counts for genes mapped to the same symbol
summed.gex= gex.df %>%
  dplyr::select(gene, everything(), -gene_id)%>%
  filter(gene != "")%>%
  group_by(gene)%>%
  summarise_all(sum)

# data conversions
df= column_to_rownames(summed.gex, "gene")

colnames(df)= substr(colnames(df), 1,11)

table(colnames(df) %in% meta$Run)

meta= meta %>% arrange(Sample.Name)

colnames(df)[match( colnames(df), meta$Run,)]

meta= meta %>% filter(Run %in% colnames(df), 
                      etiology %in% c("Non-Failing Donor", "Hypertrophic cardiomyopathy (HCM)" , "Dilated cardiomyopathy (DCM)")
)%>%
  as_tibble()

dim(df)
df= df[,meta$Run]
dim(df)

unique(meta$etiology)

dge<- DGEList(counts=df, group= meta$etiology)#, group=group)

#detect and remove low expressed gene
keep <- filterByExpr(dge, min.count	= 15, min.total.count= 20, min.prop = 0.75)
table(keep)

dge <- dge[keep,,keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge)

# use limma voom to transform count data to log2 counts per million
v <- voom(dge, plot=TRUE)
colnames(meta)

TARGET= meta %>% rename(Gender= sex,
                        Age= Age, 
                        Disease= etiology,
                        Sample= Run
                        )
TARGET = TARGET %>% 
  mutate(HeartFailure = ifelse(Disease == "Non-Failing Donor", "no", "yes"), 
         DCM  = ifelse(Disease == "Dilated cardiomyopathy (DCM)", "yes", "no"))%>%
  dplyr::select(Sample,HeartFailure, DCM, Age, Gender)


Magnet= list("GEX"= v$E,
             "TARGETS"= TARGET)


METAheart$Flam19= Magnet
METAcounts$Flam19= list("gex"= dge$counts, 
                        "target"= TARGET, 
                        "voom" = v$E)

#saveRDS(METAheart,"output/reheat1/METAheart2023_1.rds")

# add forte  --------------------------------------------------------------

#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE175764

df = read.csv("data/reheat1_data/forte2022.txt", sep = "\t")
#df2 = read.csv("update_reheat/add_data/forte2021.txt", sep = "\t")

#df= read.csv("data/GSE175764_Milena_Furtado_011718_HEART_FAILURE_vs_CONTROL-FINAL-counts.txt", sep ="\t")
target= read.csv("data/reheat1_data/forte2022_filereport_read_run_PRJNA733500_tsv.txt", sep = "\t")
df= df %>% as.data.frame() %>% column_to_rownames("X")
colnames(df) %in% target$run_accession
#colnames(df)[!colnames(df) %in% target$run_accession]= c("SRR14679384", "SRR14679390")

TARGET= target%>% 
  dplyr::rename(Sample = run_accession)       %>%
  mutate(HeartFailure= ifelse(grepl(sample_title, pattern = "Control"), "no", "yes"))%>%
  dplyr::select(Sample, HeartFailure)%>% mutate(DCM= "no")

df= df[, TARGET$Sample]

## run norm+DGE
dge<- DGEList(counts=df, group= TARGET$HeartFailure)#, group=group)

#detect and remove low expressed gene
keep <- filterByExpr(dge, min.count	= 10, min.total.count= 10, min.prop = 0.4)
table(keep)

dge <- dge[keep,,keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge)

# use limma voom to transform count data to log2 counts per million
v <- voom(dge, plot=TRUE)


PCA <- prcomp(t(v$E[,]) ,center = TRUE, scale. = F)

plot.pca = PCA$x %>%
  as.data.frame %>%
  rownames_to_column("Sample") %>%
  as_tibble()%>%
  left_join(TARGET, by= "Sample" )

p.pca = ggplot(plot.pca,aes(x= PC1, y= PC2, col= HeartFailure))+
  geom_point(size= 3)+
  theme_minimal()+
  labs(x= paste0("PC1 (",as.character(round(PCA$sdev[1]^2/sum(PCA$sdev^2)*100)),"%)"),
       y= paste("PC2 (",as.character(round(PCA$sdev[2]^2/sum(PCA$sdev^2)*100)),"%)"))+
  ggtitle(paste0(""))+
geom_text(aes(label= Sample),show.legend = FALSE)

p.pca


Forte2022= list("GEX"= v$E,
             "TARGETS"= TARGET)

#METAheart= readRDS("output/reheat1/METAheart2023_1.rds")
METAheart$Forte22= Forte2022
METAcounts$Forte22= list("gex"= dge$counts, 
                         "target"= TARGET, 
                         "voom" = v$E)

#saveRDS(METAheart,"output/reheat1/METAheart2023_1.rds")

# add wang22 / GSE203160 -----------------------------------------------------------
# load data
df = read.csv("data/reheat1_data/GSE203160.txt", sep = "\t")
df= df %>% as.data.frame() %>% column_to_rownames("X")

target = read.csv("data/reheat1_data/filereport_read_run_PRJNA838927_tsv.txt", sep = "\t")
target$run_accession %in% colnames(df)

target= target %>% mutate(HeartFailure = ifelse(grepl("Ctrl", sample_title), "no", "yes"))%>%
  rename(Sample = run_accession)%>%
  dplyr::select(Sample, HeartFailure)%>% 
  mutate(DCM = "no")

colnames(df) %in% target$Sample

df= df[, target$Sample]

## run norm+DGE
dge<- DGEList(counts=df, group= target$HeartFailure)#, group=group)

keep <- filterByExpr(dge, min.count	= 10, min.total.count= 10, min.prop = 0.4)
table(keep)

dge <- dge[keep,,keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge)

# use limma voom to transform count data to log2 counts per million
v <- voom(dge, plot=TRUE)


PCA <- prcomp(t(v$E[,]) ,center = TRUE, scale. = F)

plot.pca = PCA$x %>%
  as.data.frame %>%
  rownames_to_column("Sample") %>%
  as_tibble()%>%
  left_join(target, by= "Sample" )

p.pca = ggplot(plot.pca,aes(x= PC1, y= PC2, col= HeartFailure))+
  geom_point(size= 3)+
  theme_minimal()+
  labs(x= paste0("PC1 (",as.character(round(PCA$sdev[1]^2/sum(PCA$sdev^2)*100)),"%)"),
       y= paste("PC2 (",as.character(round(PCA$sdev[2]^2/sum(PCA$sdev^2)*100)),"%)"))+
  ggtitle(paste0(""))+
  geom_text(aes(label= Sample),show.legend = FALSE)

p.pca

Wang22= list("GEX"= v$E,
                "TARGETS"= target)

#METAheart= readRDS("output/reheat1/METAheart2023_1.rds")
METAheart$Wang22= Wang22

METAcounts$Wang22= list("gex"= dge$counts, 
                         "target"= target, 
                         "voom" = v$E)

#saveRDS(METAheart,"output/reheat1/METAheart2023_1.rds")

# add Rao22 ---------------------------------------------------------------
#data
file.path= "~/Downloads/GSE145154_RAW/"

files= list.files(file.path)

bulks= files[grepl("BRCM", files)]

data_list <- lapply(bulks, function(file) {
  file= paste0(file.path, file)
  print(file)
  #extracted_ids <- sub(".*/GSM\\d+_(BRCM-.).*", "\\2", file)
  #extracted_ids <- sub(".*/GSM\\d+_(BRCM-[^.]+)\\..*", "\\1", file)
  extracted_ids <- str_extract(file, "GSM.{7}")
  
  #extracted_ids
  data <- read.table(file, header = TRUE, sep = "\t", quote = "", comment.char = "", row.names = NULL)
  # You can customize the read.table parameters based on your data
  data<-
    data %>% as_tibble()%>%
    dplyr::select(gene_id, expected_count)%>%
    group_by(gene_id)%>%
    summarize(sum_expected_count = sum(expected_count, na.rm = TRUE))%>% 
    dplyr::rename(!!extracted_ids := sum_expected_count)
  
  return(data)
})

final_data <- data_list %>%
  reduce(full_join, by = "gene_id")

dat2 = final_data %>% 
  dplyr::rename(ensgene = gene_id)%>%
  left_join(grch38 %>% dplyr::select(ensgene, symbol), by= "ensgene")%>%
  dplyr::select(symbol, everything(),-ensgene)

dat2<-dat2%>% 
  group_by(symbol)%>%
  filter(symbol != "")%>% 
  summarise_all(sum)

count.m <- dat2%>%
  ungroup() %>% 
  mutate(across(where(is.numeric), as.integer))%>%
  as.data.frame()%>%
  column_to_rownames("symbol")%>%
  as.matrix()

meta= read.csv("/home/jan/R-projects/reheat2_pilot/data/reheat1_data/meta_rao.txt")
table(colnames(count.m) %in% meta$GEO_Accession..exp.)
##note, the blood samples are automatically removed by filtering for samples in the
## count df
target <- meta %>% 
  dplyr::select(GEO_Accession..exp., tissue_type, tissue)%>%
  filter(GEO_Accession..exp. %in% colnames(count.m))%>% 
  mutate(HeartFailure = ifelse(grepl("heart failure", tissue_type), "yes", "no"))%>%
  rename(Sample = GEO_Accession..exp.)%>%
  dplyr::select(Sample, HeartFailure)

# we found age ,sex and ef info, but it is not possible to map paitent IDs so will not 
# be used.
meta2 = read.csv("/home/jan/R-projects/reheat2_pilot/data/reheat1_data/Rao2021_bulk_meta.csv")

#order samples:
df <- count.m[, target$Sample]

## run norm+DGE
dge<- DGEList(counts=df, group= target$HeartFailure)#, group=group)

keep <- filterByExpr(dge, min.count	= 10, min.total.count= 10, min.prop = 0.4)
table(keep)

dge <- dge[keep,,keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge)

# use limma voom to transform count data to log2 counts per million
v <- voom(dge, plot=TRUE)


PCA <- prcomp(t(v$E[,]) ,center = TRUE, scale. = T)

plot.pca = PCA$x %>%
  as.data.frame %>%
  rownames_to_column("Sample") %>%
  as_tibble()%>%
  left_join(target, by= "Sample" )

p.pca = ggplot(plot.pca,aes(x= PC1, y= PC2, col= HeartFailure))+
  geom_point(size= 3)+
  theme_minimal()+
  labs(x= paste0("PC1 (",as.character(round(PCA$sdev[1]^2/sum(PCA$sdev^2)*100)),"%)"),
       y= paste("PC2 (",as.character(round(PCA$sdev[2]^2/sum(PCA$sdev^2)*100)),"%)"))+
  ggtitle(paste0(""))+
  geom_text(aes(label= Sample),show.legend = FALSE)

p.pca

## add to meta heart

Rao21= list("GEX"= v$E,
             "TARGETS"= target)

#METAheart = readRDS("output/reheat1/METAheart2023_1.rds")
METAheart$Rao21 = Rao21
METAcounts$Rao21= list("gex"= dge$counts, 
                        "target"= Rao21$TARGETS, 
                        "voom" = v$E)


# save the updated object -------------------------------------------------


saveRDS(METAheart,"output/reheat1/METAheart2023_1.rds")
saveRDS(METAcounts,"output/reheat1/Metaheart_counts2023.rds")



# add XU ------------------------------------------------------------------

x= read_delim("~/Downloads/GSE135055_raw_counts_GRCh38.p13_NCBI.tsv")

#fibrosis staining data
meta  = read.csv("data/reheat1_data/12916_2019_1469_MOESM6_ESM-.csv", skip = 1)


#sample ids
meta2 = read_tsv("data/reheat1_data/filereport_read_run_PRJNA557232_tsv.txt")
meta2 <- meta2%>% 
  mutate(Seq.ID = toupper(str_extract(sample_title, "[^-]+$")))%>%
  dplyr::select(Seq.ID, sample_alias)

target <- meta2 %>% left_join(meta, by= "Seq.ID")
#clinical data
meta3 <- read.csv("data/reheat1_data/12916_2019_1469_MOESM1_ESM.csv", skip = 2)
colnames(meta3)[1:6]<-c("pat_id", 
                        "Seq.ID", 
                        "HF", 
                        "sex", 
                        "age", 
                        "nyha")
meta3
target<- target %>% left_join(meta3, by= "Seq.ID")
target <- target %>% dplyr::select(Seq.ID, sample_alias, HF, sex, age, everything())


dat2 <- x %>%  dplyr::rename(entrez = GeneID)%>%
  left_join(grch38 %>% dplyr::select(entrez, symbol), by= "entrez")%>%
  dplyr::select(symbol, everything(),-entrez)
dat2<-dat2%>% 
  group_by(symbol)%>%
  filter(symbol != "")%>% 
  summarise_all(sum)


count.m <- dat2%>%
  ungroup() %>% 
  mutate(across(where(is.numeric), as.integer))%>%
  as.data.frame()%>%
  column_to_rownames("symbol")%>%
  as.matrix()

target.filt <- target %>% 
  filter(sample_alias %in% colnames(count.m))%>%
  mutate(HeartFailure= ifelse(HF =="DCM", "yes", "no"))%>%
  arrange(HeartFailure)

count.m <- count.m[, target.filt$sample_alias]
## run norm+DGE
dge<- DGEList(counts=count.m, group= target.filt$HF)#, group=group)

keep <- filterByExpr(dge, min.count	= 10, min.total.count= 10, min.prop = 0.4)
table(keep)

dge <- dge[keep,,keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge)

# use limma voom to transform count data to log2 counts per million
v <- voom(dge, plot=TRUE)

target_for_metaheart<- target.filt%>% 
  rename(Sample=sample_alias, 
         Gender = sex,
         Age = age)%>% 
  mutate(Gender = ifelse(Gender =="M","male",
                         ifelse(Gender=="F", "female",Gender)),
         DCM = ifelse(HF =="DCM", "yes", "no"), 
         HeartFailure= ifelse(HF =="Healthy donor", "no", "yes"))%>%
  dplyr::select(Sample, HeartFailure,DCM, Gender, Age)

METAheart= readRDS("output/reheat1/METAheart2023_1.rds")

METAheart$Hua19= list("TARGETS"= target_for_metaheart, 
                        "GEX" = v$E)
saveRDS(METAheart,"output/reheat1/METAheart2023_1.rds")

METAcounts <- readRDS("output/reheat1/Metaheart_counts2023.rds")
METAcounts$Hua2019= list("gex"= dge$counts, 
                        "target"= target, 
                        "voom" = v$E)
saveRDS(METAcounts,"output/reheat1/Metaheart_counts2023.rds")


