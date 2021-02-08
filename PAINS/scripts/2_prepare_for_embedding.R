library(plyr)
library(tidyverse)
library(arrow)

source("scripts/get_cell_metadata.R")
source("scripts/get_cell_features.R")

plate_tracking <- arrow::read_parquet("intermediate_data/plate_tracking.parquet")
plate_map <- arrow::read_parquet("intermediate_data/plate_map.parquet")
compound_map <- arrow::read_parquet("intermediate_data/compound_map.parquet")
dataset_ids <- readr::read_tsv("raw_data/dataset_ids.tsv")
cell_feature_columns <- readr::read_tsv("raw_data/cell_feature_columns.tsv")

if (!file.exists("raw_data/cell_metadata_columns.tsv")) {
    dataset_id <- "3b581cf7-b04d-4a02-a647-e3d9a6f1bc09"
    cell_metadata <- get_cell_metadata(dataset_id)
    metadata_columns <- names(cell_metadata)
    metadata_columns <- data.frame(column = metadata_columns)
    metadata_columns %>%
        readr::write_tsv("raw_data/cell_metadata_columns.tsv")
} else {
    metadata_columns <- readr::read_tsv("raw_data/cell_metadata_columns.tsv")
}

# collet 2M cells from the 48 hour time point
cell_features <- dataset_ids %>%
    dplyr::filter(
        dataset_id != "b0197b57-5e0c-4b7e-b59a-f4d68f4d768a", # Ono23@48h has 5 datasets, exclude 1
        `Treatment Duration` == "48 hours") %>%
    dplyr::select(dataset_id) %>%
    plyr::adply(1, function(dataset) {
        dataset_id <- dataset$dataset_id[1]
        cat(
            "Collecting features for compound plate: dataset_id: '", dataset_id, "'\n",
            sep = "")
        command <- paste0(
            "  mkdir raw_data/", dataset_id, "\n",
            "  sudo cp ~/bucket_cellprofilerdata/PAINS/", dataset_id, "/cpdata.h5 raw_data/", dataset_id, "/cpdata.h5\n",
            "  sudo chmod a+r raw_data/", dataset_id, "/cpdata.h5")
        cat(paste0(command, "\n"))
        system(command)
        cell_features <- get_cell_features(dataset_id, cell_feature_columns) %>%
            dplyr::sample_n(62500)
        command <- paste0("sudo rm -rf raw_data/", dataset_id)
        cat(paste0("             ", command, "\n"))
        system(command)
        return(cell_features)
    })

cell_features %>%
    arrow::write_parquet(
        sink = "raw_data/48h_2M_Cell_MasterDataTable.parquet")

# copy to S3 bucket
system("sudo cp raw_data/48h_2M_Cell_MasterDataTable.parquet ~/bucket_cellprofilerdata/PAINS/UMAP_embeddings/")

# retrieve from S3 bucket
system("sudo cp ~/bucket_cellprofilerdata/PAINS/UMAP_embeddings/48h_2M_Cell_MasterDataTable.parquet raw_data/")

# embed the 2M cells collected from the 48h hour timepoint
system("
        cd ~/opt/MPLearn/vignettes/PAINS &&
        /home/ubuntu/anaconda3/envs/sextonlab/bin/python \\
            ~/anaconda3/envs/sextonlab/bin/embed_umap \\
            --dataset /home/ubuntu/opt/MPStats/vignettes/PAINS/raw_data/48h_2M_Cell_MasterDataTable.parquet \\
            --tag UMAP_embedding_48h_2M \\
            --feature_columns ~/opt/MPStats/vignettes/PAINS/raw_data/cell_feature_columns.tsv \\
            --random_subset 2000000 \\
            --pca_batch_size 50000 \\
	    --umap_low_memory \\
	    --verbose
")

# embed each dataset into the 48h_2M embedding
dataset_ids %>%
    dplyr::select(dataset_id) %>%
    plyr::a_ply(1, function(dataset) {
        dataset_id <- dataset$dataset_id[1]
        cat(
            "embedding compounds for dataset with dataset_id: '", dataset_id, "'\n",
            sep = "")
        command <- paste0(
            "  mkdir raw_data/", dataset_id, "\n",
            "  sudo cp ~/bucket_cellprofilerdata/PAINS/", dataset_id, "/cpdata.h5 raw_data/", dataset_id, "/cpdata.h5\n",
            "  sudo chmod a+r raw_data/", dataset_id, "/cpdata.h5\n",
            "  sudo rm -rf ~/tmp/cellprofilerdata/PAINS/*")
        cat(paste0(command, "\n"))
        system(command)
        dataset_path <- paste0("/home/ubuntu/opt/MPStats/vignettes/PAINS/raw_data/", dataset_id, "_Cell_MasterDataTable.parquet")
        get_cell_features(dataset_id, cell_feature_columns) %>%
            arrow::write_parquet(dataset_path)
        
        command <- paste0("
        cd ~/opt/MPLearn/vignettes/PAINS &&
        /home/ubuntu/anaconda3/envs/sextonlab/bin/python \\
            ~/anaconda3/envs/sextonlab/bin/embed_umap \\
            --ref_embed_dir /home/ubuntu/opt/MPStats/vignettes/PAINS/intermediate_data/UMAP_embedding_48h_2M \\
            --dataset ", dataset_path, " \\
            --tag UMAP_embedding_", dataset_id, "_into_48h_2M \\
            --feature_columns ~/opt/MPStats/vignettes/PAINS/raw_data/cell_feature_columns.tsv \\
	    --verbose")
        cat(command, "\n", sep = "")
        system(command)
        
        command <- paste0(
            "sudo rm -rf raw_data/", dataset_id, "\n",
            "rm -rf ", dataset_path)
        cat("  ", command, "\n", sep = "")
        system(command)
    })


# collect the metadata for all the datasets
dataset_ids %>%
    dplyr::select(dataset_id) %>%
    dplyr::filter(dplyr::row_number() >= 15) %>%
    plyr::a_ply(1, function(dataset) {
        dataset_id <- dataset$dataset_id[1]
        cat(
            "Collecting metadata for dataset for dataset with dataset_id: '", dataset_id, "'\n",
            sep = "")
        command <- paste0(
            "  mkdir raw_data/", dataset_id, "\n",
            "  sudo cp ~/bucket_cellprofilerdata/PAINS/", dataset_id, "/cpdata.h5 raw_data/", dataset_id, "/cpdata.h5\n",
            "  sudo chmod a+r raw_data/", dataset_id, "/cpdata.h5\n",
            "  sudo rm -rf ~/tmp/cellprofilerdata/PAINS/*")
        cat(paste0(command, "\n"))
        system(command)
        dataset_path <- paste0("/home/ubuntu/opt/MPStats/vignettes/PAINS/raw_data/", dataset_id, "_Cell_MasterDataTable.parquet")
        cell_metadata <- get_cell_metadata(dataset_id) %>%
            arrow::write_parquet(dataset_path)

        command <- paste0("sudo rm -rf raw_data/", dataset_id)
        cat("  ", command, "\n", sep = "")
        system(command)
    })




##############################################################################
# Re-do embedding with greater negative sampling rate and better convergence #
##############################################################################

# embed the 2M cells collected from the 48h hour timepoint
system("
        cd ~/opt/MPLearn/vignettes/PAINS &&
        /home/ubuntu/anaconda3/envs/sextonlab/bin/python \\
            ~/anaconda3/envs/sextonlab/bin/embed_umap \\
            --dataset /home/ubuntu/opt/MPStats/vignettes/PAINS/raw_data/48h_2M_Cell_MasterDataTable.parquet \\
            --tag UMAP_embedding_48h_2M_20200816 \\
            --feature_columns ~/opt/MPStats/vignettes/PAINS/raw_data/cell_feature_columns_20200814.tsv \\
            --random_subset 2000000 \\
            --pca_batch_size 50000 \\
            --umap_negative_sample_rate 20 \\
            --umap_n_epochs 2000 \\
	    --umap_low_memory \\
	    --verbose
")
# note that this requires somewhere between 347-547G and ~10 hours
# and more than 225G of space free space to save the transformation
# as such I failed to have enough space to keep the transformation
