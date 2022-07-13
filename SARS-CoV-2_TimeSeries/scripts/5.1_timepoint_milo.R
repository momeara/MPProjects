



library(tidyverse)
library(arrow)
library(brms)
library(miloR)

cell_feature_columns <- readr::read_tsv(
    "product/cell_feature_columns_TS_202008.tsv",
    show_col_types = FALSE)

cell_metadata_columns <- readr::read_tsv(
    "product/cell_metadata_columns_TS_202008.tsv",
    show_col_types = FALSE)


k <- 100
resolution <- 1e-4
num_iter <- 10
cell_metadata <- dplyr::bind_cols(
    arrow::read_parquet(
        file = "product/TS2_plate_scaled_Cell_MasterDataTable.parquet",
        col_select = cell_metadata_columns$feature),
    arrow::read_parquet(
        file = paste0(
        "intermediate_data/UMAP_embedding_into_TS2_2M_epochs=2000_20200901/",
        "clusters_leiden_k=", k, "_res=", resolution, "_num_iter=", num_iter, ".parquet"))) %>%
    dplyr::mutate(
        time_point = factor(
            x = time_point,
            levels = c(
                "Uninfected",
                "8 hours", "18 hours", "24 hours", "30 hours", "36 hours", "48 hours"),
            labels = c(
                "Uninfected",
                "8 hours", "18 hours", "24 hours", "30 hours", "36 hours", "48 hours")),
        well_id = paste0(toupper(letters[column]), row)) %>%
    dplyr::arrange(cluster_label) %>%
    dplyr::mutate(
        cluster_code = cluster_label,
        cluster_label = paste("cluster ", cluster_label) %>% as.factor() %>% forcats::fct_inorder())


cell_embedding <- arrow::read_parquet(
    file = "intermediate_data/UMAP_embedding_into_TS2_2M_epochs=2000_20200901/umap_embedding.parquet"),

cell_milo <- cell_embedding %>%
    dplyr::select(UMAP1, UMAP2) %>%
    as.matrix() %>%
    t() %>%
    miloR::buildGraph(k = 1000, d = 10) %>%
    miloR::countCells(samples = "samples", meta.data = cell_metadata) %>%
    miloR::calcNhoodDistance(d = 10) %>%

milo_design <- cell_metadata %>%
    stats::xtabs(formula = ~ time_point + well_id, .) %>%
    as.data.frame() %>%
    dplyr::filter(Freq > 0)
    
da_milo <- milo %>%
    miloR::testNhoods(
        design = ~ time_point,
        design.df = milo_design)
