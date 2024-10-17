library(MOFAcellulaR)
library(tidyverse)
library(cowplot)
library(ComplexHeatmap)
library(biclust)
# 1. For a loading matrix, identify the number of genes that
# belong to 1 or more cell-types (by sign)

# 2. For genes with 1 cell-type, separate them for ORA

# 3. For genes with more than 1, separate them for ORA

# First define Hypergeometric analysis

GSE_analysis <- function(geneList,Annotation_DB){

  geneList = geneList[geneList %in% unique(unlist(Annotation_DB))]

  ResultsDF = matrix(0,nrow = length(Annotation_DB),ncol = 5)
  rownames(ResultsDF) = names(Annotation_DB)
  colnames(ResultsDF) = c("GenesInPathway","GenesInList","GeneNames",
                          "p_value","corr_p_value")

  DB_genecontent = length(unique(unlist(Annotation_DB)))
  GenesDB = DB_genecontent
  SelectedGenes = length(geneList)

  for(gset in rownames(ResultsDF)){
    GP = length(Annotation_DB[[gset]])
    GL = length(intersect(Annotation_DB[[gset]],geneList))

    ResultsDF[gset,"GenesInList"] = GL
    ResultsDF[gset,"GenesInPathway"] = GP
    ResultsDF[gset,"GeneNames"] = paste(intersect(Annotation_DB[[gset]],
                                                  geneList),
                                        collapse = ",")
    ResultsDF[gset,"p_value"] = phyper(q=GL - 1, m=GP, n=GenesDB-GP,
                                       k=SelectedGenes, lower.tail = FALSE,
                                       log.p = FALSE)
  }

  ResultsDF[,"corr_p_value"] = p.adjust(ResultsDF[,"p_value"],method = "BH")
  ResultsDF = data.frame(ResultsDF,stringsAsFactors = F)
  ResultsDF = ResultsDF[order(ResultsDF[,"p_value"]),]

  ResultsDF = ResultsDF %>%
    rownames_to_column("gset") %>%
    mutate_at(c("GenesInPathway","GenesInList",
                "p_value","corr_p_value"),
              as.numeric) %>%
    dplyr::arrange(corr_p_value,GenesInList)

  return(ResultsDF)
}

# Get genesets of interest

msigdb <- decoupleR::get_resource("MSigDB") %>%
  dplyr::filter(!grepl("COMPLEX", genesymbol))

msigdb_hallmarks <- msigdb %>%
  dplyr::filter(collection == "hallmark") %>%
  dplyr::select(genesymbol, geneset) %>%
  dplyr::mutate(weight = 1) %>%
  dplyr::rename("source" = geneset,
                "target" = genesymbol) %>%
  unique()

msigdb_canonical <- msigdb %>%
  dplyr::filter(collection %in% c("kegg_pathways", "reactome_pathways",
                                  "biocarta_pathways")) %>%
  dplyr::select(genesymbol, geneset) %>%
  dplyr::mutate(weight = 1) %>%
  dplyr::rename("source" = geneset,
                "target" = genesymbol) %>%
  unique()

gsets <- list("CM" = read_csv("./data/gsets/pathophysiol_processed/cm.csv"),
              "ECM" = read_csv("./data/gsets/pathophysiol_processed/mural.csv"),
              "VASC" = read_csv("./data/gsets/pathophysiol_processed/vasc.csv"),
              "HMARKS" = msigdb_hallmarks,
              "CAN" = msigdb_canonical)

gsets <- map(gsets, function(gset_col) {
  gset_col %>%
    dplyr::select(-weight) %>%
    group_by(source) %>%
    nest() %>%
    dplyr::mutate(data = map(data, ~.x[[1]])) %>%
    deframe()

})

# Check genes
loadings <- read_csv("./results/meta_mofacell/gene_loadings.csv") %>%
  dplyr::filter(Factor %in% c("Factor1", "Factor2"))


for(Fact in c("Factor1", "Factor2")) {

  unique_sets_file <- paste0("./results/meta_mofacell/", Fact,"_",
                             "single_cell_enrichment_v2.csv")

  multicell_sets_file <- paste0("./results/meta_mofacell/", Fact,"_",
                                "multicell_enrichment_v2.csv")

  split_loadings <- loadings %>%
    dplyr::filter(Factor == Fact) %>%
    dplyr::mutate(direction = ifelse(value > 0, "pos", "neg")) %>%
    group_by(direction) %>%
    nest()

  # Let's get cell-type unique dicts

  unique_sets <- map2(set_names(split_loadings$direction),
                      split_loadings$data,
                      function(dir, dat) {

                        filt_dat <- dat %>%
                          ungroup() %>%
                          dplyr::filter(abs(value) >= 0.1) %>%
                          group_by(feature) %>%
                          dplyr::mutate(n_reps = n()) %>%
                          dplyr::filter(n_reps == 1) %>%
                          dplyr::select(feature, ctype) %>%
                          group_by(ctype) %>%
                          nest() %>%
                          dplyr::mutate(data = map(data, ~.x[[1]])) %>%
                          deframe()

                      })

  # Enrichments by hand

  unique_acts <- map(unique_sets, function(dir_sets) {

    # Cardios
    cardio_gsets <- GSE_analysis(geneList = dir_sets[["CM"]],
                                 Annotation_DB = gsets$CM) %>%
      #dplyr::filter(corr_p_value <= 0.1) %>%
      dplyr::mutate(ct = "CM") %>%
      mutate(collection = "muscular")

    # ECM
    ECMcts <- c("Fib", "PC", "Endo", "vSMCs")
    ECMcts <- ECMcts[ECMcts %in% names(dir_sets)]

    ecm_gsets <- map(ECMcts, function(ct) {

      ECM <- GSE_analysis(geneList = dir_sets[[ct]],
                          Annotation_DB = gsets$ECM) %>%
        #dplyr::filter(corr_p_value <= 0.1) %>%
        dplyr::mutate(ct = ct)


    }) %>%
      bind_rows() %>%
      mutate(collection = "ecm")

    # VASC
    VASCcts <- c("PC", "Endo", "vSMCs")
    VASCcts <- VASCcts[VASCcts %in% names(dir_sets)]

    vasc_gsets <- map(VASCcts, function(ct) {


      ECM <- GSE_analysis(geneList = dir_sets[[ct]],
                          Annotation_DB = gsets$VASC) %>%
        #dplyr::filter(corr_p_value <= 0.1) %>%
        dplyr::mutate(ct = ct)


    }) %>%
      bind_rows() %>%
      mutate(collection = "vasculature")

    # All cells

    cts <- names(dir_sets)

    hmark_gsets <- map(cts, function(ct) {

      ECM <- GSE_analysis(geneList = dir_sets[[ct]],
                          Annotation_DB = gsets$HMARKS) %>%
        #dplyr::filter(corr_p_value <= 0.1) %>%
        dplyr::mutate(ct = ct)


    }) %>%
      bind_rows() %>%
      mutate(collection = "hallmark")

    all_sets <- bind_rows(hmark_gsets, vasc_gsets, ecm_gsets, cardio_gsets) %>%
      dplyr::arrange(ct) %>%
      dplyr::select(gset, corr_p_value, ct,GenesInPathway, GenesInList, collection)
  })


  enframe(unique_acts, "direction") %>%
    unnest() %>%
    dplyr::filter(collection != "canonical") %>%
    write_csv(unique_sets_file)

  # Multicellular sets

  multicell_sets <- map2(set_names(split_loadings$direction),
                      split_loadings$data,
                      function(dir, dat) {

                        filt_dat <- dat %>%
                          ungroup() %>%
                          dplyr::filter(abs(value) >= 0.1) %>%
                          group_by(feature) %>%
                          dplyr::mutate(n_reps = n()) %>%
                          dplyr::filter(n_reps > 1) %>%
                          ungroup() %>%
                          dplyr::select(feature) %>%
                          unique() %>%
                          deframe()

                        return(filt_dat)

                      })

  mcell_acts <- map(multicell_sets, function(gset) {

    can_gsets <- GSE_analysis(geneList = gset,
                              Annotation_DB = gsets$CAN) %>%
      #dplyr::filter(corr_p_value <= 0.1) %>%
      mutate(collection = "canonical")

    hmark_gsets <- GSE_analysis(geneList = gset,
                                Annotation_DB = gsets$HMARKS) %>%
      #dplyr::filter(corr_p_value <= 0.1) %>%
      mutate(collection = "hmarks")

    vasc_gsets <- GSE_analysis(geneList = gset,
                               Annotation_DB = gsets$VASC) %>%
      #dplyr::filter(corr_p_value <= 0.1) %>%
      mutate(collection = "vasculature")

    muscle_gsets <- GSE_analysis(geneList = gset,
                                 Annotation_DB = gsets$CM) %>%
      #dplyr::filter(corr_p_value <= 0.1) %>%
      mutate(collection = "muscular")

    ecm_gsets <- GSE_analysis(geneList = gset,
                              Annotation_DB = gsets$ECM) %>%
      #dplyr::filter(corr_p_value <= 0.1) %>%
      mutate(collection = "ecm")

    all_sets <- bind_rows(hmark_gsets, can_gsets,
                          vasc_gsets, muscle_gsets,
                          ecm_gsets) %>%
      dplyr::select(gset, corr_p_value, GenesInPathway, GenesInList, collection)

    return(all_sets)

  })


  enframe(mcell_acts, "direction") %>%
    unnest() %>%
    dplyr::filter(collection != "canonical") %>%
    write_csv(multicell_sets_file)


}





