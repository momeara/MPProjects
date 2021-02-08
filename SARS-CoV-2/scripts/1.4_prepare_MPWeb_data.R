


library(plyr)
library(tidyverse)
library(readxl)
library(RMySQL)
library(magrittr)
library(tictoc)
library(arrow)


load("intermediate_data/image_scores_CX5_100X.Rdata")


z <- image_scores_CX5_100X %>%
    dplyr::mutate(
        well_id = paste0(Image_Metadata_PlateID, "_", Image_Metadata_WellID),
        plate_alias = plate_id,
        `Smiles String`= `Smiles String` %>%
             stringr::str_replace("\r", "") %>%
             stringr::str_replace("\n", "")) %>%
    dplyr::group_by(well_id) %>%
    dplyr::summarize(
        project_id = "SARS-CoV-2_S25",
        project_description = "SARS-CoV-2 Drug Repurposing Collection 5-plate QHTS",
        plate_id = Image_Metadata_PlateID[1],
        plate_alias = plate_alias[1],
        row = Image_Metadata_WellID[1] %>% stringr::str_extract("^."),
        column=column[1],
        treatment = Compound[1],
        treatment_description = paste0("Therapeutic Class:", `Therapeutic Class`[1], ",smiles:", `Smiles String`[1],",CDR:",CDR[1]),
        dose=dose_nM[1],
        dose_units="nM",
        replicate=1,
        image_paths=paste0(
            "https://sextonpublicimages.s3.amazonaws.com/SARS-COV-2/", Image_Metadata_PlateID, "/pics/",
            Image_Metadata_PlateID, "_",
            Image_Metadata_WellID, "_",
            Image_Metadata_Field, ".jpeg",
            collapse=";"))
z %>% readr::write_tsv("product/MPWeb_SARS-CoV-2_CX5_100X_well_map_200509.tsv")


image_scores_1999B <- arrow::read_parquet("product/image_scores_1999B.parquet")
z <- image_scores_1999B %>%
    dplyr::mutate(
        well_id = paste0(Image_Metadata_PlateID, "_", Image_Metadata_WellID),
        plate_alias = plate_id) %>%
    dplyr::group_by(well_id) %>%
    dplyr::summarize(
        project_id = "SARS-CoV-2_S25",
        project_description = "SARS-CoV-2 Drug Repurposing Collection Lactoferrin vs. Remdesivir",
        plate_id = Image_Metadata_PlateID[1],
        plate_alias = plate_alias[1],
        row = Well_ID[1] %>% stringr::str_extract("^."),
        column=column[1],
        treatment = paste0(rem_label[1], lf_label[1]),
        treatment_description = paste0("Condition:", Condition[1], ",Cell Line:Huh-7"),
        dose=lf_label[1],
        dose_units="Lactoferrin µg/mL",
        replicate=1,
        image_paths=paste0(
            "https://sextonpublicimages.s3.amazonaws.com/SARS-COV-2/", Image_Metadata_PlateID, "/pics/",
            Image_Metadata_PlateID, "_",
            Image_Metadata_WellID, "_",
            Image_Metadata_FieldID, ".jpeg",
            collapse=";"))
z %>% readr::write_tsv("product/MPWeb_SARS-CoV-2_1999B_well_map_200509.tsv")

image_scores_20XX <- arrow::read_parquet("product/image_scores_20XX.parquet")
z <- image_scores_20XX %>%
    dplyr::mutate(
        well_id = paste0(Image_Metadata_PlateID, "_", Image_Metadata_WellID),
        plate_alias = plate_id) %>%
    dplyr::group_by(well_id) %>%
    dplyr::summarize(
        project_id = "SARS-CoV-2_S25",
        project_description = "SARS-CoV-2 Drug Repurposing Collection 10x Dose Response",
        plate_id = Image_Metadata_PlateID[1],
        plate_alias = plate_alias[1],
        row = Image_Metadata_WellID[1] %>% stringr::str_extract("^."),
        column=column[1],
        treatment = Compound_Name[1],
        treatment_description = paste0("Microscpe:CX5,Condition:", Condition[1], ",Cell Line:Huh-7"),
        dose=dose_nM[1],
        dose_units="nM",
        replicate=1,
        image_paths=paste0(
            "https://sextonpublicimages.s3.amazonaws.com/SARS-COV-2/", Image_Metadata_PlateID, "/pics/",
            Image_Metadata_PlateID, "_",
            Image_Metadata_WellID, "_",
            Image_Metadata_Field, ".jpeg",
            collapse=";"))
z %>% readr::write_tsv("product/MPWeb_SARS-CoV-2_20XX_well_map_200509.tsv")


