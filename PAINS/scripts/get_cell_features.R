
library(plyr)
library(tidyverse)

source("scripts/get_cell_metadata.R")


get_cell_features <- function(
        dataset_id,
        cell_feature_columns,
        verbose = TRUE) {
    cell_metadata <- get_cell_metadata(dataset_id)
    metadata_columns <- names(cell_metadata)

    dataset_path <- paste0("raw_data/", dataset_id, "/cpdata.h5")
    if (verbose) {
        cat("Getting cell features for dataset from path '", dataset_path, "' ...\n", sep = "")
    }

    dataset <- hdf5r::H5File$new(dataset_path, "r")
    object_columns <- tibble::tibble(feature = dataset[["object/cFeatureName"]][])
    object_features <- dataset[["object/M"]][, ]
    colnames(object_features) <- object_columns$feature
    object_features <- object_features %>%
        tibble::as_tibble()

    if (verbose) {
        cat("Number of feature columns: ", cell_feature_columns %>% nrow(), "\n")
        n_total_cells <- object_features %>% nrow()
        cat("Total number of cells: ", n_total_cells, "\n", sep = "")
    }
    cell_features <-
        dplyr::bind_cols(
            cell_metadata,
            object_features) %>%
        dplyr::filter(
            Cell_Children_Cytoplasm_Count > 0,
            Nucleus_Children_Cells_Count > 0,
            Nucleus_Children_Cytoplasm_Count > 0) %>%
        dplyr::select(
            tidyselect::one_of(c(
                metadata_columns,
                cell_feature_columns$feature))) %>%
        dplyr::filter_at(
            .vars = vars(tidyselect::one_of(cell_feature_columns$feature)),
            ~ !is.na(.))

    if (verbose) {
        n_filtered_cells <- cell_features %>% nrow()
        cat("Filtered number of cells ", n_filtered_cells, "\n", sep = "")
    }
    cell_features
}
