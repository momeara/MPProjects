library(plyr)
library(tidyverse)
library(arrow)
library(monocle3)
library(BiocNeighbors)
library(MPStats)
library(miloR)

# gathered in 3_pseudo_time_202006_UMAP_embedding
cell_features <- arrow::read_parquet(
    "raw_data/covid19cq1_TS_scaled_Cell_MasterDataTable.parquet")

embedding <- arrow::read_parquet(
    "~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/UMAP_embedding_TS_scaled_full/umap_embedding.parquet")


time_points <- cell_features$time_point %>% unique()


milo_design <- cell_metadata %>%
    stats::xtabs(formula = ~ time_point + plate_id, .) %>%
    as.data.frame() %>%
    dplyr::filter(Freq > 0)


milo_test <- MPStats::populate_csd(
    cell_features = cell_features,
    cell_feature_columns = cell_feature_columns,
    cell_metadata_columns = cell_metadata_columns)                      
    miloR::Milo() %>%
    miloR::buildGraph(k = 20, d = 30) %>%
    miloR::makeNhoods(k = 20, d = 30, refined = TRUE, prop = 0.2) %>%
    miloR::countCells(samples = "Sample", meta.data = cell_meta) %>%
    miloR::testNhoods(design = ~ Condition, design.df = milo_design)




 



cds <- MPStats::populate_cds(
    cell_features = cell_features,
    cell_feature_columns = cell_feature_columns,
    cell_metadata_columns = cell_metadata_columns,
    embedding = embedding)







