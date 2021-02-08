
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

feature_columns <- readr::read_tsv(
    "raw_data/cell_feature_columns.tsv")

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
    dplyr::bind_cols(
        cell_features,
        cluster_labels %>% dplyr::select(cluster_label))


clofazamine_test_data <- cell_features %>%
    dplyr::filter(Compound %in% c(
        "PC",
        "NC",
        "Clofazamine",
        "Niclosamide",
        "Amiodarone (hydrochloride)",
        "Lomitapide",
        "Z-FA-FMK",
        "S1RA",
        "Eliglustat")) %>%
    dplyr::mutate(
        test_region = dplyr::case_when(
            cluster_label == "ROI 2" ~ "ROI 2",
            cluster_label %in% c("ROI 1", "ROI 3", "ROI 4") ~ "ROI 1,3,4",
            TRUE ~ NA_character_)) %>%
    dplyr::filter(!is.na(test_region)) %>%
    dplyr::count(Compound, dose_nM, test_region) %>%
    tidyr::pivot_wider(
        id_cols = c("Compound", "dose_nM"),
        names_from = "test_region",
        values_from = "n") %>%
    dplyr::mutate(
        roi_2_fraction = `ROI 2` / (`ROI 1,3,4` + `ROI 2`),
        roi_2_fraction_sd = sqrt((`ROI 1,3,4` + `ROI 2`) * roi_2_fraction * (1 - roi_2_fraction)) / (`ROI 1,3,4` + `ROI 2`))
                      


plot <- ggplot2::ggplot(
    data = clofazamine_test_data %>%
        dplyr::filter(`ROI 1,3,4` > 20) %>%
        dplyr::filter(Compound != "PC", Compound != "NC")) +
    ggplot2::theme_bw() +
    ggplot2::geom_hline(
        yintercept = 0.14286413,
        color = "blue") +
    ggplot2::geom_text(
        data = data.frame(label = c("Negative Control"), y = 0.14286413 + .01, x = 8),
        mapping = ggplot2::aes(
            label = label,
            x = log10(x),
            y = y)) +
    ggplot2::geom_line(
        mapping = ggplot2::aes(
            x = log10(dose_nM),
            y = roi_2_fraction,
            color = Compound),
        size = 1.5) +
    ggplot2::geom_errorbar(
        mapping = ggplot2::aes(
            x = log10(dose_nM),
            ymin = roi_2_fraction - roi_2_fraction_sd,
            ymax = roi_2_fraction + roi_2_fraction_sd,
            color = Compound)) +
    ggplot2::geom_point(
        mapping = ggplot2::aes(
            x = log10(dose_nM),
            y = roi_2_fraction,
            color = Compound)) +
    ggplot2::scale_x_continuous(
        "Dose nM",
        breaks = log10(c(3.996004, 7.992008, 15.984016, 31.168831, 62.337662, 124.675325, 255.744256, 495.504496,  991.008991, 1998.001998)),
        labels = c("4", "8", "16", "31", "63", "125", "256", "495", "991", "1998")) +
    ggplot2::scale_y_continuous(
        "Fraction of infected cells in ROI 2",
        limits = c(0, 0.5))

ggsave(filename="product/figures/umap_features/clofazamine_roi2_20200625_2.pdf", width=6, height=6)
