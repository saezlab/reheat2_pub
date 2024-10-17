# Copyright (c) [2023] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we quantify the degree of
#' overlap of cell-type markers

library(tidyverse)

setwd("/Users/ricardoramirez/Dropbox/PostDoc/Research/ReHeaT2/")
outplot_file <- "./results/qc/mrkr_overlap.pdf"

edgeR_df <- tibble(file = list.files("./data/mrkrs/")) %>%
  dplyr::mutate(study = gsub("_mrkrs.csv","",file)) %>%
  dplyr::filter(study != "Reichart2022_DCM") %>%
  dplyr::mutate(file = paste0("./data/mrkrs/",file)) %>%
  dplyr::mutate(edgeR = map(file, read_csv, show_col_types = FALSE)) %>%
  unnest(c(edgeR)) %>%
  dplyr::select(-file)

mrkr_list <- edgeR_df %>%
  dplyr::filter(FDR < 0.01, logFC > 2) %>%
  dplyr::rename("lineage" = name) %>%
  group_by(study,lineage) %>%
  nest() %>%
  dplyr::mutate(data = map(data, ~.x[[1]])) %>%
  dplyr::mutate(key = paste0(study,"..", lineage)) %>%
  dplyr::ungroup() %>%
  dplyr::select(key, data) %>%
  group_by(key) %>%
  deframe()

jaccard_ix_mat <- sapply(mrkr_list, function(x){
  sapply(mrkr_list, function(y){
    length(intersect(x,y))/length(union(x,y))
  })
})

#jaccard_ix_mat[upper.tri(jaccard_ix_mat,diag = T)] <- NA

jaccard_ix_df <- jaccard_ix_mat %>%
  as.data.frame() %>%
  rownames_to_column("StudyA") %>%
  pivot_longer(-StudyA, names_to = "StudyB") %>%
  #na.omit() %>%
  dplyr::mutate(ctA = strsplit(StudyA, "[.][.]") %>%
                  map_chr(., ~.x[[2]]),
                ctB = strsplit(StudyB, "[.][.]") %>%
                  map_chr(., ~.x[[2]])) %>%
  dplyr::filter(StudyA != StudyB)

jix_plt <- jaccard_ix_df %>%
  mutate(align = "al") %>%
  ggplot(., aes(y = align, x = value)) +
  geom_boxplot() +
  geom_point(size = 0.5) +
  facet_grid(ctB ~ ctA,switch = "y") +
  ylab("") +
  xlab("Jaccard Index") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(breaks=c(0,0.5,1)) +
  theme(strip.text.y.left = element_text(angle = 0))

pdf(outplot_file, height = 2.5, width = 5.5)
plot(jix_plt)
dev.off()


jaccard_ix_df %>%
  dplyr::mutate(type_int = ifelse(ctA == ctB, "reciprocal", "not_reciprocal")) %>%
  group_by(type_int) %>%
  summarize(median(value))




