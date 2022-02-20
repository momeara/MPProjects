
library(plyr)
library(tidyverse)

source("parameters.R")

activities_summary <- readr::read_tsv("raw_data/activities_summary_20210323.tsv")
substances <- readr::read_tsv("intermediate_data/substances_sanitized_20210323.tsv")
activities <- readr::read_tsv("raw_data/activities_20210323.tsv")

# check for substances with alternate names
alternate_substance_names <- activities_summary %>%
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
    nrow(alternate_substance_names) == 0,
    msg = "Substances with the same smiles have distinct substance names.")

# check input/rdkit smiles
duplicate_smiles_rdkit <- substances %>%
    dplyr::distinct(substance_smiles, .keep_all = TRUE) %>%
    dplyr::select(substance_name, substance_smiles, substance_smiles_rdkit) %>%
    dplyr::inner_join(
        substances %>%
            dplyr::distinct(substance_smiles, .keep_all = TRUE) %>%
            dplyr::count(substance_smiles_rdkit) %>%
            dplyr::filter(n > 1),
        by = c("substance_smiles_rdkit")) %>%
    dplyr::filter(!is.na(substance_smiles_rdkit)) %>%
    dplyr::arrange(substance_smiles_rdkit)
assertthat::assert_that(
    nrow(duplicate_smiles_rdkit) == 0,
    msg = "Substances with different input smiles have the same sanitized smiles.")


# check dock_ids are distinct
duplicate_dock_id_substances <- substances %>%
    dplyr::count(substance_dock_id) %>%
    dplyr::filter(n > 1)
assertthat::assert_that(
    nrow(duplicate_dock_id_substances) == 0,
    msg = "Substances with distinct smiles have duplicate dock_ids")



unidentified_activities <- activities %>%
    dplyr::anti_join(activities_summary, by = "substance_name") %>%
    dplyr::distinct(substance_name, .keep_all = TRUE) %>%
    dplyr::arrange(substance_name)

