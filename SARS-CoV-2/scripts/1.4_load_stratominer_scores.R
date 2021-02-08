


library(plyr)
library(tidyverse)
library(ggplot2)
library(readxl)
library(MPStats)
library(arrow)


stratominer_well_scores <- readr::read_tsv(
    file="raw_data/Stratominer_2016A-2019A_all_hits_200514.txt") %>%
    dplyr::mutate(
       plate_id = Barcode,
       Plate_Name = paste0("SARS_", plate_id),
    row = WellID %>%
        stringr::str_extract("^[A-Z]") %>%
        purrr::map_int(~which(LETTERS==., arr.ind=T)),
    column = WellID %>%
        stringr::str_extract("[0-9]+$") %>%
        as.integer(),
    dose_nM = Concentration * 1000)

stratominer_well_scores %>% arrow::write_parquet(
    file="intermediate_data/stratominer_well_scores_2016A-2019A_200514.parquet")


# these are just the well scores with the feature level averages?
statominer_results <- readr::read_tsv("raw_data/Stratominer_2006A-2009A_all_results_200514.txt")    


########################################3

stratominer_well_scores_1999B_2020A <- readr::read_csv(
    file="raw_data/Stratominer_1999B_2020A_200515.csv") %>%
    dplyr::mutate(
        Barcode = ifelse(Barcode == "1999A", "1999B", Barcode),
        plate_id = Barcode,
        Plate_Name = paste0("SARS_", plate_id),
    row = wellLocation %>%
        stringr::str_extract("^[a-zA-Z]") %>%
        stringr::str_to_upper() %>%
        purrr::map_int(~which(LETTERS==., arr.ind=T)),
    column = wellLocation %>%
        stringr::str_extract("[0-9]+$") %>%
        as.integer()) %>%
    dplyr::rename(Remdesivir_Concentration = RemdesivirConcentration) %>%
    dplyr::mutate(Remdesivir_Concentration = ifelse(is.na(Remdesivir_Concentration), 0, Remdesivir_Concentration))

stratominer_well_scores_1999B_2020A %>%
    arrow::write_parquet("intermediate_data/stratominer_well_scores_1999B_2020A_200515.parquet")
