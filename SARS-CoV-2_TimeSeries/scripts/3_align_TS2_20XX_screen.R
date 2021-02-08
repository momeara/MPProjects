library(plyr)
library(tidyverse)
library(arrow)
library(monocle3)
library(BiocNeighbors)
library(MPStats)

# gathered in 3_pseudo_time_202006_UMAP_embedding
cell_features <- arrow::read_parquet(
    "raw_data/covid19cq1_TS_scaled_Cell_MasterDataTable.parquet")

embedding <- arrow::read_parquet(
    "~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/UMAP_embedding_TS_scaled_full/umap_embedding.parquet")


time_points <- cell_features$time_point %>% unique()


    


cds <- MPStats::populate_cds(
    cell_features = cell_features,
    cell_feature_columns = cell_feature_columns,
    cell_metadata_columns = cell_metadata_columns,
    embedding = embedding)







