# Copyright (c) [2024] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Check if the multicellular signature relates to
#' pediatric DCM
#'

library(MOFAcellulaR)
library(tidyverse)
library(cowplot)
source("./code/reheat2_pilot/aesthetics.R")

# Main ###############################################################
# outs
mds_out <- "./results/predict_project/PED/ped_fplot.pdf"
factor_bplots <- "./results/predict_project/PED/ped_factors.pdf"
proj_w_file <- "./results/predict_project/PED/ped_proj_w.csv"
ttest_file <- "./results/predict_project/PED/ped_ttest_res.csv"
aov_file <- "./results/predict_project/PED/ped_aov_res.csv"
factor_coords  <- "./results/predict_project/PED/ped_factor_coords.csv"

# Multi-cell model
mcell_reheat <- MOFA2::load_model("./results/meta_mofacell/metamodel_mofa.hdf5")

# New data
all_data <- readRDS("./data/val_studies/pediatric/pediatricDCM_pb_sample_celltype.rds")

# Counts
pb_data <- all_data$pb

# Col-data
coldat <- all_data$target  %>%
  column_to_rownames("pb_id") %>%
  dplyr::rename(ncells = "cell_count") %>%
  dplyr::filter(cell_type != "none") %>%
  dplyr::rename(sample = "orig.ident") %>%
  dplyr::mutate(group = ifelse(group == "CT", "NF", group))

pb_data <- pb_data[, rownames(coldat)]

# Meta data
meta_data <- coldat[, c("sample", "group")] %>%
  unique() %>%
  as_tibble()

# Processing the test data
# Defining cts
cts <- coldat$cell_type %>%
  unique() %>%
  set_names()

### Defining function to make projection and test it against meta

assoc_model <- function(MOFAcell_obj) {

  test_projection <- project_data(model = mcell_reheat,
                                  test_data = MOFAcell_obj)

  test_assocs_aov <- test_projection %>%
    as.data.frame(.) %>%
    rownames_to_column("sample") %>%
    pivot_longer(-sample, names_to = "Factor", values_to = "Score") %>%
    left_join(meta_data) %>%
    group_by(Factor) %>%
    nest() %>%
    mutate(aov_assocs = map(data, function(dat) {
      aov(Score ~ group, data = dat) %>%
        broom::tidy()
    })) %>%
    dplyr::select(Factor, aov_assocs) %>%
    unnest() %>%
    dplyr::filter(term == "group") %>%
    ungroup() %>%
    dplyr::filter(Factor %in% c("Factor1", "Factor2")) %>%
    dplyr::mutate(adj_p = p.adjust(p.value, "BH"))

  test_assocs_t <- test_projection %>%
    as.data.frame(.) %>%
    rownames_to_column("sample") %>%
    pivot_longer(-sample, names_to = "Factor", values_to = "Score") %>%
    left_join(meta_data) %>%
    dplyr::filter(Factor %in% c("Factor1", "Factor2")) %>%
    group_by(Factor) %>%
    nest() %>%
    mutate(pair_ttests = map(data, function(dat) {
      pairwise.t.test(dat$Score, dat$group,
                      p.adjust.method = "BH") %>%
        broom::tidy()
    })) %>%
    dplyr::select(Factor, pair_ttests) %>%
    unnest()

  f_interest <- c("Factor1", "Factor2")

  plots <- map(set_names(f_interest), function(fact) {

    useful_factors <- test_projection %>%
      as.data.frame(.) %>%
      rownames_to_column("sample") %>%
      pivot_longer(-sample, names_to = "Factor", values_to = "Score") %>%
      left_join(meta_data) %>%
      dplyr::filter(Factor %in% c(fact)) %>%
      ggplot(aes(x = group, y = Score, fill = group)) +
      geom_boxplot() +
      #geom_point(size = 3, alpha = 0.8) +
      theme_bw() +
      theme(axis.text = element_text(size =12),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      scale_fill_manual(values = etiology_colors[meta_data$group %>% unique]) +
      ylab(fact)
  })

  fplot <- test_projection %>%
    as.data.frame(.) %>%
    rownames_to_column("sample") %>%
    left_join(meta_data) %>%
    ggplot2::ggplot(aes(x = Factor1, y = Factor2, color = group)) +
    geom_point() +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 10),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank()) +
    scale_color_manual(values = etiology_colors[meta_data$group %>% unique])

  return(list("assocs_aov" = test_assocs_aov,
              "assocs_t" = test_assocs_t,
              "bplots" = plots,
              "fplots" = fplot))

}

# add complete model
# Creating MOFAcell object to make projection
MOFAcell_obj <- MOFAcellulaR::create_init_exp(pb_data, coldat) %>%
  MOFAcellulaR::filt_profiles(pb_dat = .,
                              cts = cts,
                              ncells = 10,
                              counts_col = "ncells",
                              ct_col = "cell_type") %>%
  MOFAcellulaR::tmm_trns(pb_dat_list = .,
                         scale_factor = 1000000) %>%
  MOFAcellulaR::center_views(pb_dat_list = .) %>%
  MOFAcellulaR::pb_dat2MOFA(pb_dat_list = .,
                            sample_column = "sample")

all_dat <- assoc_model(MOFAcell_obj)

write_csv(all_dat$assocs_t,ttest_file)
write_csv(all_dat$assocs_aov,aov_file)

# Make plots
pdf(mds_out, height = 2.5, width = 4.0)
plot(all_dat$fplots)
dev.off()

pdf(factor_bplots, width = 2.8, height = 2.5)
plot(all_dat$bplots$Factor1 + ylab("MCP1") + ylim(-0.7,.95))
plot(all_dat$bplots$Factor2)
dev.off()

# Projection

test_projection <- project_data(model = mcell_reheat,
                                test_data = MOFAcell_obj)

out_test <- test_projection[,c("Factor1", "Factor2")] %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  left_join(meta_data, by = c("sample_id"= "sample"))

write_csv(out_test, factor_coords)

