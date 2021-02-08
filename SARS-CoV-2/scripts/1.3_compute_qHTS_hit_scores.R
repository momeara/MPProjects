
library(plyr)
library(tidyverse)
library(ggplot2)
library(readxl)
library(MPStats)
library(arrow)


infectivity_cell_scores <- readr::read_tsv("raw_data/plate_ids.tsv") %>%
    plyr::adply(1, function(df) {
        cat("reading plate '", df$plate_id[1], "'\n", sep = "")
        arrow::read_parquet(
            file = paste0("product/SARS_", df$plate_id[1], "_Cell_MasterDataTable.parquet"),
            col_select = c(
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
        infectivity_well_score = mean(infectivity_score / infectivity_plate_score)) %>%
    dplyr::ungroup()

infectivity_well_scores %>%
    arrow::write_parquet("intermediate_data/S25_infectivity_well_scores.parquet")


###################################################3


umap_well_scores <- dplyr::bind_cols(
    dplyr::bind_rows(
        arrow::read_parquet(
            file = "product/SARS_10030250A_Cell_MasterDataTable.parquet",
            col_select = c("Image_Metadata_WellID", "Compound")) %>%
            dplyr::mutate(dose_nM = 250),
        arrow::read_parquet(
            file = "product/SARS_10030500A_Cell_MasterDataTable.parquet",
            col_select = c("Image_Metadata_WellID", "Compound")) %>%
            dplyr::mutate(dose_nM = 500),
        arrow::read_parquet(
            file = "product/SARS_10031000A_Cell_MasterDataTable.parquet",
            col_select = c("Image_Metadata_WellID", "Compound")) %>%
            dplyr::mutate(dose_nM = 1000),
        arrow::read_parquet(
            file = "product/SARS_10032000A_Cell_MasterDataTable.parquet",
            col_select = c("Image_Metadata_WellID", "Compound")) %>%
            dplyr::mutate(dose_nM = 2000)))

    
