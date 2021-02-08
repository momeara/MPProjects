

library(plyr)
library(tidyverse)
library(MPStats)
library(arrow)

dataset_ids <- readr::read_tsv("raw_data/dataset_ids.tsv")

cell_feature_columns <- readr::read_tsv("raw_data/cell_feature_columns.tsv")
cell_metadata_columns <- readr::read_tsv("raw_data/cell_metadata_columns.tsv")


embedding <- dataset_ids %>%
    dplyr::select(dataset_id) %>%
    plyr::adply(1, function(dataset) {
        dataset_id <- dataset$dataset_id[1]
        embedding_path <- paste0("intermediate_data/per_plate_into_48h_2M/cell_metadata/", dataset_id, "_cell_metadata.parquet")
        arrow::read_parquet(embedding_path)
    }) 



cell_feature_columns <- tibble::tibble(feature = c("UMAP_1", "UMAP_2"))

embedding <- embedding %>% dplyr::sample_n(5000000)


cell_dataset <- MPStats::populate_cds(
    cell_features = embedding,
    cell_feature_columns = cell_feature_columns,
    cell_metadata_columns = cell_metadata_columns,
    embedding_type = c("UMAP"),
    embedding = embedding %>% dplyr::select(UMAP_1, UMAP_2),
    verbose = TRUE)

cell_dataset <- cell_dataset %>%
    monocle3::cluster_cells(
        k = 200,
        resolution = .00001,
        num_iter = 10,
        verbose = TRUE)

# this can't handle all the cells at once, so look into how the
# clustering can be done on a subset of the cells and then extended to
# the rest
