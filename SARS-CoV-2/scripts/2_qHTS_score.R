
library(plyr)
library(tidyverse)
library(ggplot2)
library(readxl)
library(MPStats)
library(arrow)


infectivity_cell_scores <- readr::read_tsv("raw_data/plate_ids.tsv") %>%
    plyr::adply(1, function(df){
        cat("reading plate '", df$plate_id[1], "'\n", sep="")
        arrow::read_parquet(
            file=paste0("product/SARS_", df$plate_id[1], "_Cell_MasterDataTable.parquet"),
            col_select=c(
                "plate_id",
                "dose_nM",
                "Compound",
                "row",
                "column",
                "Image_Metadata_Field",
                "Cells_Intensity_IntegratedIntensityEdge_Virus",
                "Cells_Intensity_MeanIntensityEdge_Virus",
                "Cells_Intensity_MaxIntensityEdge_Virus",
                "Cells_Intensity_MaxIntensity_Virus"))
    }) %>%
    dplyr::mutate(
        infectivity_score =
            -5.064328 +
            Cells_Intensity_IntegratedIntensityEdge_Virus * 1.487025e-01 +
            Cells_Intensity_MeanIntensityEdge_Virus * -3.840196e+01 +
            Cells_Intensity_MaxIntensityEdge_Virus * 4.270269e+01 +
            Cells_Intensity_MaxIntensity_Virus * 4.254849e+01)

infectivity_cell_scores %>%
    arrow::write_parquet("intermediate_data/S25_infectivity_cell_scores.parquet")

infectivity_well_scores <- infectivity_cell_scores %>%
    dplyr::group_by(plate_id) %>%
    dplyr::mutate(infectivity_plate_score = mean(infectivity_score)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(master_plate_id, plate_id, dose_nM, Compound, row, column) %>%
    dplyr::summarize(
        infectivity_well_score = mean(infectivity_score/infectivity_plate_score)) %>%
    dplyr::ungroup()

infectivity_well_scores %>%
    arrow::write_parquet("intermediate_data/S25_infectivity_well_scores.parquet")


#' Does the qHTS screen show a consistent trend in inhibition?
#' Input well scores with
#'     well_features %>%
#'        tidyr::pivot_wider(
#'             names_from=dose_nM,
#'             names_prefix="dose_",
#'             values_from=score)
#'
qHTS_has_trend <- function(well_features, score){
        with(well_features,
             (dose_50 > dose_250)  & (dose_250 > dose_500) |
             (dose_250 > dose_500) & (dose_500 > dose_1000) |
             (dose_500 > dose_1000) & (dose_1000 > dose_2000))
}

#' 
qHTS_significant_inhibition <- function(well_features, threshold=.5){
    with(well_features,
        (dose_50 < threshold) |
        (dose_250 < threshold) |       
        (dose_500 < threshold) |
        (dose_1000 < threshold) |
        (dose_2000 < threshold))
}


