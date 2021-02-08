
library(plyr)
library(tidyverse)
library(arrow)
library(monocle3)
library(MPStats)
library(miloR)
library(scater)

n_cells <- 100000

cell_feature_columns <- readr::read_tsv(
    "product/cell_feature_columns_TS_202008.tsv")

cell_metadata_columns <- readr::read_tsv(
    "product/cell_metadata_columns_TS_202008.tsv")


cell_features <- dplyr::bind_cols(
    arrow::read_parquet(
        file = "product/TS2_2M_Cell_MasterDataTable.parquet"),
    arrow::read_parquet(
        file = "intermediate_data/UMAP_embedding_TS2_2M_epochs=2000_20200901/umap_embedding.parquet"))
cell_features <- cell_features %>% dplyr::sample_n(n_cells)

cell_metadata <- cell_features %>%
    dplyr::select(tidyselect::one_of(cell_metadata_columns$feature))


#######################
# Plate Batch Effects #
#######################

plate_ids <- cell_features$plate_id %>% unique()

milo <- cell_features %>%
    dplyr::select(tidyselect::one_of(cell_feature_columns$feature)) %>%
    as.matrix() %>%
    t() %>%
    miloR::buildGraph(k = 2000, d = 10)    

milo <- milo %>% scater::runPCA(ncomponents=50)
milo <- milo %>% scater::runUMAP()

plot <- milo %>% scater::plotUMAP()
ggplot2::ggsave(
    filename = paste0("product/figures/TS2/milo/umap_", n_cells, "_cells_20210128.pdf"),
    plot = plot) 


a <- milo %>% miloR::plotNhoodSizeHist()
ggplot2::ggsave(
    filename = paste0("product/figures/TS2/milo/Nhood_size_hist", n_cells, "_20210128.pdf"))

cell_metadata <- cell_metadata %>% dplyr::mutate(
    samples = paste(time_point, plate_id, sep = "_"))

milo <- milo %>% miloR::countCells(samples = "samples", meta.data = cell_metadata)

milo <- milo %>% miloR::calcNhoodDistance(d=10)

milo_design <- cell_metadata %>%
    stats::xtabs(formula = ~ time_point + plate_id, .) %>%
    as.data.frame() %>%
    dplyr::filter(Freq > 0)
    
da_milo <- milo %>%
    miloR::testNhoods(
        design = ~ time_point,
        design.df = milo_design)

milo <- milo %>% miloR::buildNhoodGraph()

z <- scater::plotUMAP(milo) +
    miloR::plotNhoodGraphDA(milo, da_milo, alpha=1, layout="UMAP") +
    patchwork::plot_layout(guides="collect")

ggplot2::ggsave(
    paste0("product/figures/TS2/milo/umap_nhoods_", n_cells, "_20210128.pdf"),
    width = 10,
    height = 5)
