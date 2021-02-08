
library(plyr)
library(tidyverse)
library(ggplot2)
library(readxl)
library(MPStats)
library(arrow)

plate_ids <- readr::read_tsv("raw_data/plate_ids.tsv")


### CX5
infectivity_score_CX5_features <- c(
    "Cells_Intensity_IntegratedIntensityEdge_Virus",
    "Cells_Intensity_MeanIntensityEdge_Virus",
    "Cells_Intensity_MaxIntensityEdge_Virus",
    "Cells_Intensity_MaxIntensity_Virus")
add_infectivity_score_CX5 <- function(cell_features){
    cell_features %>%
        dplyr::mutate(
            infectivity_score = -5.064328 +
                Cells_Intensity_IntegratedIntensityEdge_Virus * 1.487025e-01 +
                Cells_Intensity_MeanIntensityEdge_Virus * -3.840196e+01 +
                Cells_Intensity_MaxIntensityEdge_Virus * 4.270269e+01 +
                Cells_Intensity_MaxIntensity_Virus * 4.254849e+01)
}


infectivity_score_field <- plate_ids %>%
    dplyr::filter(
        schema == 'covid19primary',
        plate_id %>% stringr::str_detect("^100.")) %>%
    plyr::adply(1, function(df) {
        plate_id <- df$plate_id[1]
        cell_features <- arrow::read_parquet(
            file = paste0(
                "product/covid19primary_SARS_",
                plate_id,
                "_Cell_MasterDataTable.parquet"),
            col_select = c(
                "Compound",
                "row",
                "column",
                "Image_Metadata_Field",
                tidyselect::all_of(infectivity_score_CX5_features)))
    }) %>%
    add_infectivity_score_CX5() %>%
    dplyr::group_by(plate_id) %>%
    dplyr::mutate(
        infectivity_score_infected_mean = ifelse(
                !(Compound %in% c("Positive Control")),
                infectivity_score,
                NA) %>%
            mean(na.rm = TRUE),
        infectivity_score_pc_mean = ifelse(
                Compound %in% c("Positive Control"),
                infectivity_score,
                NA) %>%
            mean(na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
        infectivity_score_cell =
            (infectivity_score - infectivity_score_pc_mean) /
            (infectivity_score_infected_mean - infectivity_score_pc_mean)) %>%
    dplyr::group_by(
        plate_id,
        dose_nM,
        Compound,
        row,
        column,
        Image_Metadata_Field) %>%
    dplyr::summarize(
        infectivity_score_field = mean(infectivity_score_cell, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(plate_id, dose_nM, Compound, row, column) %>%
    dplyr::arrange(infectivity_score_field) %>%
    dplyr::mutate(
        infectivity_score_field_rank = dplyr::row_number() / dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
        keep_field = infectivity_score_field_rank >= .66)

infectivity_score_field %>% arrow::write_parquet(
    "product/infectivity_score_field_10XX_200519.parquet")

infectivity_score_well <- infectivity_score_field %>%
    dplyr::filter(infectivity_score_field_rank >= .66) %>%
    dplyr::group_by(plate_id, row, column, dose_nM, Compound) %>%
    dplyr::summarize(
        infectivity_score_well_mean = mean(infectivity_score_field),
        infectivity_score_well_sem =
            sd(infectivity_score_field) / sqrt(dplyr::n())) %>%
    dplyr::ungroup()
infectivity_score_well %>% arrow::write_parquet(
    "product/infectivity_score_well_10XX_200519.parquet")



infectivity_score_treatment <- infectivity_score_field %>%
    dplyr::filter(infectivity_score_field_rank >= .66) %>%
    dplyr::group_by(plate_id, dose_nM, Compound) %>%
    dplyr::summarize(
        infectivity_score_treatment_mean = mean(infectivity_score_field),
        infectivity_score_treatment_sem =
            sd(infectivity_score_field) / sqrt(dplyr::n())) %>%
    dplyr::ungroup()

infectivity_score_treatment %>% arrow::write_parquet(
    "product/infectivity_score_treatment_10XX_200519.parquet")
