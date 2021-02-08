

library(plyr)
library(tidyverse)
library(MPStats)

metadata <- readr::read_csv(
  "raw_data/iAEC2_Combo_Metadata_20201001.csv",
  col_types = readr::cols(
    .default = readr::col_logical(),         
    Plate = readr::col_double(),
    `Dispensed Well` = readr::col_character(),
    `Fluid name` = readr::col_character(),
    Concentration = readr::col_double(),
    Lomitapide = readr::col_double(),
    Remdesivir = readr::col_double(),
    S1RA = readr::col_double(),
    Silmitasertib = readr::col_double(),
    `Z-FA-FMK` = readr::col_double(),
    Amiodarone = readr::col_double(),
    Apilimod = readr::col_double(),
    Camostat = readr::col_double(),
    Clofazimine = readr::col_double(),
    Hydroxychloroquine = readr::col_double(),
    `Ipratropium Bromide` = readr::col_double(),
    Dynasore = readr::col_double(),
    Fluvoxamine = readr::col_double(),
    Terfernadine = readr::col_double(),
    Amiloride = readr::col_double(),
    Sertraline = readr::col_double())) %>%
  dplyr::select(
    -`Plate ID`,               
    -tidyselect::starts_with("X")) %>%
  dplyr::rename(
    well_id = `Dispensed Well`,
    plate_id = Plate)


lf_metadata <- metadata %>%
  dplyr::filter(
    !is.na(`Fluid name`),
    `Fluid name` != "2 Fluids") %>%
  tidyr::pivot_longer(
     cols = c(5:20),
     names_to="drug_2",
     values_to="drug_concentration_2",
     values_drop_na = TRUE) %>%
  dplyr::transmute(
     plate_id,
     well_id,
     drug_1 = "Lactoferrin",
     drug_concentration_1 = Concentration,
     drug_2,
     drug_concentration_2)
    
two_fluids_metadata <- metadata %>%
  dplyr::filter(
    !is.na(`Fluid name`),
    `Fluid name` == "2 Fluids") %>%
  dplyr::select(-Concentration) %>%
  tidyr::pivot_longer(
    cols = c(4:19),
    names_to="drug",
    values_to="drug_concentration",
    values_drop_na = TRUE) %>%
  dplyr::group_by(plate_id, well_id) %>%
  dplyr::mutate(drug_index = dplyr::row_number()) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(
    id_cols = c(plate_id, well_id),
    names_from = drug_index,         
    values_from = c(drug, drug_concentration))         


tidy_metadata <- dplyr::bind_rows(
  lf_metadata,
  two_fluids_metadata)

tidy_metadata %>%
  readr::write_tsv(
    "product/iAEC2_Combo_Metadata_tidy_20201001.tsv")        
