
library(plyr)
library(tidyverse)
library(googlesheets4)
#library(dm)


source("parameters.R")


# check substance_name collisions
duplicate_substance_names <- activities_summary %>%
    dplyr::distinct(substance_name, .keep_all = TRUE) %>%
    dplyr::count(substance_smiles) %>%
    dplyr::filter(n > 1) %>%
    dplyr::left_join(
        activities_summary %>% dplyr::select(
            substance_name,
            substance_smiles,
            by = "substance_name")) %>%
    data.frame()
assertthat::assert_that(
    nrow(duplicate_substance_names) == 0,
    msg = "Substances with the same smiles have distinct substance names.")


# check dock_ids are distinct
duplicate_dock_id_substances <- substances %>%
    dplyr::count(substance_dock_id) %>%
    dplyr::filter(n > 1)
assertthat::assert_that(
    nrow(duplicate_dock_id_substances) == 0,
    msg = "Substances with distinct smiles have duplicate dock_ids")


# check input/rdkit smiles
duplicate_smiles_rdkit <- substances %>%
    dplyr::distinct(substance_smiles, .keep_all = TRUE) %>%
    dplyr::select(substance_name, substance_smiles, substance_smiles_rdkit) %>%
        dplyr::inner_join(
            substances %>%
                dplyr::distinct(substance_smiles, .keep_all = TRUE) %>%
                dplyr::count(substance_smiles_rdkit) %>%
                dplyr::filter(n>1),
        by = c("substance_smiles_rdkit"))
