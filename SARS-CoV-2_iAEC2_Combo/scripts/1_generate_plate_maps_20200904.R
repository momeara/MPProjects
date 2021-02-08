

library(plyr)
library(tidyverse)
library(MPStats)



treatments <- readr::read_tsv("raw_data/iAEC2_combination_treatments_20200904.tsv") %>%
    dplyr::mutate(IC50 = log10(IC50_nM) - 9)

wells <- treatments %>%
    MPStats::generate_combination_plate_maps(
        n_rays_per_combination = 1,
        n_doses_per_ray = 8,
        n_replicas_per_dose = 3)
    
# which plates are which treatments on?
dplyr::bind_rows(
    wells %>% dplyr::select(plate_index, treatment = treatment_1),
    wells %>% dplyr::select(plate_index, treatment = treatment_2)) %>%
    dplyr::count(treatment, plate_index) %>%
    dplyr::filter(!is.na(treatment)) %>%
    data.frame

wells %>% readr::write_tsv(
    "product/iAC2_plate_maps_20200904.tsv")


