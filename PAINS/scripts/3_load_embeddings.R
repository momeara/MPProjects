
library(plyr)
library(tidyverse)
library(arrow)

source("scripts/get_cell_features.R")

source("scripts/mount_S3_bucket.R")

mount_S3_bucket()

dataset_ids <- readr::read_tsv(
    "raw_data/dataset_ids.tsv")

command <- paste0("
sudo cp -r ~/bucket_cellprofilerdata/PAINS/UMAP_embeddings/per_plate_into_48h_2M intermediate_data/
sudo chmod a+r intermediate_data/per_plate_into_48h_2M
sudo chmod a+w intermediate_data/per_plate_into_48h_2M
mkdir intermediate_data/per_plate_into_48h_2M/cell_metadata
")
system(command)

cell_feature_columns <- readr::read_tsv("raw_data/cell_feature_columns.tsv")
cell_metadata_columns <- readr::read_tsv("raw_data/cell_metadata_columns.tsv")

# collect the metadata for all the datasets
dataset_ids %>%
    dplyr::select(dataset_id) %>%
    plyr::a_ply(1, function(dataset) {
        dataset_id <- dataset$dataset_id[1]
        output_path <- paste0("intermediate_data/per_plate_into_48h_2M/cell_metadata/", dataset_id, "_cell_metadata.parquet")
        if (file.exists(output_path)) {
            cat("Output path: '", output_path, "' already exists, skipping...\n", sep = "")
            return(NULL)
        }
        cat(
            "Collecting metadata for dataset for dataset with dataset_id: '", dataset_id, "'\n",
            sep = "")
        if (!file.exists(paste0("raw_data/", dataset_id, "/cpdata.h5"))) {
            command <- paste0(
                "  mkdir raw_data/", dataset_id, "\n",
                "  sudo cp ~/bucket_cellprofilerdata/PAINS/", dataset_id, "/cpdata.h5 raw_data/", dataset_id, "/cpdata.h5\n",
                "  sudo chmod a+r raw_data/", dataset_id, "/cpdata.h5\n",
                "  sudo rm -rf ~/tmp/cellprofilerdata/PAINS/*")
            cat(paste0(command, "\n"))
            system(command)
        }
        cell_features <- get_cell_features(dataset_id, cell_feature_columns) %>%
            dplyr::select(
                tidyselect::one_of(cell_metadata_columns$column),
                tidyselect::matches("_Location_"),
                tidyselect::matches("ObjectNumber"))
        command <- paste0("sudo rm -rf raw_data/", dataset_id)
        cat("  ", command, "\n", sep = "")
        system(command)
        embedding <- arrow::read_parquet(paste0("intermediate_data/per_plate_into_48h_2M/UMAP_embedding_", dataset_id, "_into_48h_2M/umap_embedding.parquet"))
        full_embedding <- dplyr::bind_cols(
            cell_features,
            embedding)
        full_embedding %>% arrow::write_parquet(
            paste0("intermediate_data/per_plate_into_48h_2M/cell_metadata/", dataset_id, "_cell_metadata.parquet"))
    })


embedding <- dataset_ids %>%
    dplyr::select(dataset_id) %>%
    plyr::adply(1, function(dataset) {
        dataset_id <- dataset$dataset_id[1]
        embedding_path <- paste0("intermediate_data/per_plate_into_48h_2M/cell_metadata/", dataset_id, "_cell_metadata.parquet")
        arrow::read_parquet(embedding_path)
    })

embedding %>% arrow::write_parquet(
    paste0("intermediate_data/per_plate_into_48h_2M/cell_metadata.parquet"))
