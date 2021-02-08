library(plyr)
library(tidyverse)
library(fuzzyjoin)
library(ggplot2)
library(readxl)
library(MPStats)
library(arrow)
library(viridis)

library(caret)
library(gbm)

cell_features_20XX <- c(
    '2006A', '2007A', '2008A', '2009A',
    '2010A', '2010A',          '2012A',
    '2013A', '2014A', '2015A', '2016A',
    '2017A',          '2019A') %>%
    plyr::ldply(function(plate_id) {
        cat("Reading features from plate ", plate_id, " ...\n", sep="")
        arrow::read_parquet(
            paste0("product/covid19cq1_SARS_", plate_id, "_plate_scaled_Cell_MasterDataTable.parquet"),
            col_select = c("plate_id", "Compound", "dose_nM"))
    })


cell_features_lf_rem <- arrow::read_parquet(
    "product/lf_rem_plate_scaled_1999B_2020A_Cell_MasterDataTable.parquet",
    col_select=c("plate_id", "Condition", "Compound", "Lactoferrin_Concentration", "Remdesivir_Concentration"))


cell_features <- dplyr::bind_rows(
    cell_features_20XX,
    cell_features_lf_rem)

cluster_membership <- arrow:::read_parquet(
     "~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/top_hits_plate_scaled_200522a_umap2_2M_15_0.0/fig3_ROI_membership.parquet")


cluster_labels <- cluster_membership %>%
    dplyr::mutate(cell_index = dplyr::row_number()) %>%
    tidyr::pivot_longer(
        cols = tidyselect::starts_with("roi_"),
        names_to  = "cluster_label") %>%
    dplyr::arrange(desc(value)) %>%
    dplyr::distinct(cell_index, .keep_all = TRUE) %>%
    dplyr::arrange(cell_index) %>%
    dplyr::mutate(
        cluster_label =
            ifelse(!value, "ROI -1", cluster_label) %>%
            stringr::str_replace("roi_", "ROI "),
        cluster_index = cluster_label %>%
            stringr::str_replace("ROI ", "") %>%
            as.numeric(),
        cluster_index = dplyr::case_when(
            cluster_index == 0 ~ 1,
            cluster_index == 1 ~ 2,
            cluster_index == 2 ~ 3,
            cluster_index == 3 ~ 4,
            cluster_index == 4 ~ -1,
            TRUE ~ cluster_index),
        cluster_index = ifelse(cluster_index == -1, Inf, cluster_index),
        cluster_label = paste0("ROI ", cluster_index) %>% factor() %>% reorder(cluster_index)) %>%
    dplyr::select(-value)


cell_features <- cell_features %>%
    dplyr::bind_cols(cell_features, cluster_labels)


roi_membership_summary <- cell_features %>%
    dplyr::mutate(overall_total = dplyr::n()) %>%
    dplyr::group_by(cluster_label) %>%
    dplyr::mutate(overall_roi = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(cluster_label, Treatment, Compound) %>%

roi_membership_summary <- cell_features %>%
    plyr::ddply("cluster_label", function(df){
        data.frame(
            overall_total = nrow(cell_features),
            overall_roi = nrow(df),
            nc_total = cell_features %>%
                dplyr::filter(Compound == "NC" | Compound == "Negative Control") %>%
                dplyr::nrow(),
            nc_roi = df %>%
                dplyr::filter(Compound == "NC" | Compound == "Negative Control") %>%
                dplyr::nrow()) %>%
            dplyr::group_by(
                



cat("Number of cells: ", nrow(cell_features), "\n", sep="")
cat("Number of cells by ROI:\n")
cell_features %>% dplyr::count(cluster_label)


treatment_by_cluster <- dplyr::bind_cols(
    cell_features,
    cluster_labels) %>%
    dplyr::count(plate_id, Condition, Compound, dose_nM, Lactoferrin_Concentration, Remdesivir_Concentration, cluster_label) %>%
    dplyr::group_by(plate_id, Compound, dose_nM) %>%
    dplyr::mutate(fraction = n/sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
        dose_nM_label = dose_nM %>% signif(2) %>% as.character,
        cluster_label = factor(cluster_label),
        Compound = Compound %>% stringr::str_extract("^[a-zA-Z0-9-]+"))

treatment_by_cluster %>%
    readr::write_tsv(
        "product/figures/umap_features/top_hits_plate_scaled_200522a_no_virus_umap2_2M_15_0.0/fig3_treatment_by_roi.tsv")
