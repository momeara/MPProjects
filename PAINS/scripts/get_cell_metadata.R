

library(plyr)
library(tidyverse)

plate_tracking <- arrow::read_parquet("intermediate_data/plate_tracking.parquet")
plate_map <- arrow::read_parquet("intermediate_data/plate_map.parquet")
compound_map <- arrow::read_parquet("intermediate_data/compound_map.parquet")
dataset_ids <- readr::read_tsv("raw_data/dataset_ids.tsv")


get_cell_metadata <- function(dataset_id, verbose = TRUE) {
    dataset_path <- paste0("raw_data/", dataset_id, "/cpdata.h5")
    if (verbose) {
        cat("Getting cell metadata for dataset from path '", dataset_path, "' ...\n", sep = "")
    }
    dataset <- hdf5r::H5File$new(dataset_path, "r")

    cell_metadata <- tibble::tibble(
        object_number = dataset[['object/rObjectNumber']][],
        image_number = dataset[['object/rImageNumber']][])
    if (verbose) {
        cat("n cells: ", nrow(cell_metadata), "\n", sep = "")
    }

    if (verbose) {
        cat("Joining cell -> image\n")
    }
    image_metadata <- dataset[['image/rMetadata']][] %>%
        dplyr::select(
            row = Metadata_rownum,
            column = Metadata_colnum,
            site = Metadata_site) %>% # I think this is the field?
        dplyr::mutate(
            image_number = dataset[['image/rImageNumber']][])
    cell_metadata <- cell_metadata %>%
        dplyr::left_join(image_metadata, by = c("image_number"))

    if (verbose) {
        cat("Joining image -> dataset\n")
    }
    cell_metadata <- cell_metadata %>%
        dplyr::mutate(dataset_id = dataset_id) %>%
        dplyr::left_join(dataset_ids, by = c("dataset_id"))

    if (verbose) {
        cat("Joining dataset -> plate map\n")
    }
    cell_metadata <- cell_metadata %>%
        dplyr::left_join(
            plate_tracking %>%
            dplyr::select(
                compound_plate = `Compound Plate`,                   # BR00109091
                compound_plate_nickname = `Compound Plate Nickname`, # Ono20
                plate_map_name = `Plate Map for Mathias`,            # JLD012020-Ono20
                plate_barcode = `Assay Plate`,                       # BR00110363
                time_point = `Treatment Duration`),                  # 24 hours
            by = c("plate_barcode"))

    if (verbose) {
        cat("Joining plate map -> compound\n")
    }
    cell_metadata <- cell_metadata %>%
        dplyr::left_join(
            plate_map %>%
            dplyr::select(
                plate_map_name = `Plate Map Name`,   # JLD012020-Ono20
                well_position = `Well Position`,     # A01
                compound_id = `Broad Sample`,        # BRD-K00259736-001-17-2
                dose_mg_per_ml = `mg Per ml`,        # 0.133884969
                dose_mM = `mmoles Per Liter`,        # 5
                solvent = Solvent) %>%               # DMSO
            dplyr::mutate(
                row = well_position %>%
                    stringr::str_extract("^[A-Z]") %>%
                    purrr::map_int(~which(LETTERS == ., arr.ind = T)),
                column = well_position %>%
                    stringr::str_extract("[0-9]+$") %>%
                    as.integer()) %>%
            dplyr::select(-well_position),
            by = c("plate_map_name", "row", "column"))

    if (verbose) {
        cat("Joining compound -> compound info\n")
    }
    cell_metadata <- cell_metadata %>%
        dplyr::left_join(
            compound_map %>%
            dplyr::group_by(compound_id) %>%
            dplyr::summarize(
                compound_set = compound_set %>% na.omit() %>% unique() %>% paste0(collapse = "|"),
                class = class %>% na.omit() %>% unique() %>% paste0(collapse = "|"),
                alias = alias %>% na.omit() %>% unique() %>% paste0(collapse = "|"),
                smiles = na.omit(smiles)[1]),
            by = c("compound_id"))

    cell_metadata
}
