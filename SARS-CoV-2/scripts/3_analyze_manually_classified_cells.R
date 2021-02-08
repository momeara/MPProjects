
library(plyr)
library(tidyverse)
library(arrow)


load("intermediate_data/image_scores.Rdata")

manual_cells <- readr::read_csv("raw_data/manual_classification_200427.csv")
outlier_fields <- readr::read_csv("raw_data/1003_OutlierFields.csv") %>%
    dplyr::mutate(Image_Metadata_Field = as.numeric(Image_Metadata_Field)) %>%
    dplyr::mutate(
        row = FDA_384_Well_ID %>%
            stringr::str_extract("^[A-Z]") %>%
            purrr::map_int(~which(LETTERS == ., arr.ind = T)),
        column = FDA_384_Well_ID %>%
            stringr::str_extract("[0-9]+$") %>%
            as.integer()) %>%
    dplyr::mutate(master_plate_id = 1003)

z <- outlier_fields %>%
    dplyr::inner_join(
       cell_features %>%
           dplyr::mutate(
               Image_Metadata_Field = as.numeric(Image_Metadata_Field)),
       by = c(
           "Metadata_PlateID" = "Image_Metadata_PlateID",
           "row", "column", "Compound", "Image_Metadata_Field"))

z %>%
    dplyr::mutate(master_plate_id = 1003) %>%
    dplyr::select(
        master_plate_id,
        plate_id = Barcode,
        Metadata_PlateID,
        row,
        column,
        Compound,
        Image_Metadata_Field,
        cell_index) %>%
    readr::write_tsv("product/S25_plate_1003_outliers_200428.tsv")


master_plate_id <- 1003
cell_features <- dplyr::bind_rows(
    arrow::read_parquet(
        file = paste0("product/SARS_", master_plate_id, "0050A_Cell_MasterDataTable.parquet")) %>%
        dplyr::mutate(dose_nM = 50),
    arrow::read_parquet(
        file = paste0("product/SARS_", master_plate_id, "0250A_Cell_MasterDataTable.parquet")) %>%
        dplyr::mutate(dose_nM = 250),
    arrow::read_parquet(
        file = paste0("product/SARS_", master_plate_id, "0500A_Cell_MasterDataTable.parquet")) %>%
        dplyr::mutate(dose_nM = 500),
    arrow::read_parquet(
        file = paste0("product/SARS_", master_plate_id, "1000A_Cell_MasterDataTable.parquet")) %>%
        dplyr::mutate(dose_nM = 1000),
    arrow::read_parquet(
        file = paste0("product/SARS_", master_plate_id, "2000A_Cell_MasterDataTable.parquet")) %>%
        dplyr::mutate(dose_nM = 2000)) %>%
    dplyr::mutate(cell_index = dplyr::row_number())
