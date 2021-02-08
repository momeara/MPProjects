library(plyr)
library(tidyverse)
library(fuzzyjoin)
library(ggplot2)
library(readxl)
library(MPStats)
library(arrow)
library(viridis)

cell_features <- arrow::read_parquet(
    "product/lf_rem_plate_scaled_1999B_2020A_Cell_MasterDataTable.parquet",
    col_select = c(
        "plate_id",
        "row",
        "column",
        "Condition",
        "Remdesivir_Concentration",
        "Lactoferrin_Concentration"))

cluster_membership <- arrow:::read_parquet(
    "~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/covid19cq1_SARS_1999B_200523_umap2_into_top_hits_plate_scaled_200522a_15_0.0/lf_10a_membership.parquet")

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
        cluster_index = ifelse(cluster_index == -1, Inf, cluster_index),
        cluster_label =
            paste0("ROI ", cluster_index) %>%
            factor() %>%
            reorder(cluster_index)) %>%
    dplyr::select(-value)


treatment_by_cluster <- dplyr::bind_cols(
    cell_features,
    cluster_labels) %>%
    dplyr::count(
        plate_id,
        row,
        column,
        Condition,
        Remdesivir_Concentration,
        Lactoferrin_Concentration,
        cluster_label,
        name = "well_count") %>%
    dplyr::group_by(
        plate_id,
        Condition,
        Remdesivir_Concentration,
        Lactoferrin_Concentration,
        cluster_label) %>%
    dplyr::summarize(
        well_count_mean = mean(well_count)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
        rem_dose_label = signif(Remdesivir_Concentration, 2) %>%
            factor() %>%
            reorder(Remdesivir_Concentration),
        # for Lf: 1 ug/mL = 14.49 nM
        lf_dose_label = signif(Lactoferrin_Concentration * 11.49, 2) %>%
            factor() %>%
            reorder(Lactoferrin_Concentration))


p <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.background = ggplot2::element_rect(
            fill = "#00274C", colour = "#00274C"),
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom") +
    ggplot2::geom_tile(
        data = treatment_by_cluster %>% dplyr::filter(cluster_label == "ROI 0"),
        mapping = ggplot2::aes(
            x = lf_dose_label,
            y = rem_dose_label,
            fill = well_count_mean ^ (1 / 2.7))) +
    ggplot2::facet_wrap(~plate_id) +
    ggplot2::scale_x_discrete("Latoferrin dose (nM)") +
    ggplot2::scale_y_discrete("Remdesivir dose (nM)") +
    viridis::scale_fill_viridis(
        "Cell count",
        breaks = c(0, 10, 25, 50) ^ (1 / 2.7),
        labels = c("0", "10", "25", "50"))

ggplot2::ggsave(
    filename = "product/figures/umap_features/lf_rem_plate_scaled_1999B_2020A/lf_10a_classification_by_plate.pdf",
    plot = p,
    height = 4,
    width = 7)



cell_features <- arrow::read_parquet(
    "product/lf_rem_plate_scaled_1999B_2020A_Cell_MasterDataTable.parquet",
    col_select = c(
        "ImageNumber",
        "Image_Metadata_WellID",
        "Image_Metadata_FieldID",
        "plate_id",
        "Condition",
        "Remdesivir_Concentration",
        "Lactoferrin_Concentration",
        "Cells_Number_Object_Number"))


cluster_cell_examples <- dplyr::bind_cols(
    cell_features,
    cluster_labels) %>%
    dplyr::group_by(cluster_label) %>%
    dplyr::sample_n(40) %>%
    dplyr::mutate(
        rem_dose_label = signif(Remdesivir_Concentration, 2) %>%
            factor() %>%
            reorder(Remdesivir_Concentration),
        # for Lf: 1 ug/mL = 14.49 nM
        lf_dose_label = signif(Lactoferrin_Concentration * 11.49, 2) %>%
            factor() %>%
            reorder(Lactoferrin_Concentration))


cluster_cell_examples %>%
    dplyr::filter(cluster_label == "ROI 0") %>%
    readr::write_tsv("product/figures/umap_features/lf_rem_plate_scaled_1999B_2020A/lf_10a_example_cells.tsv")
