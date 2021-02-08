



library(plyr)
library(tidyverse)
library(fuzzyjoin)
library(ggplot2)
library(MPStats)
library(arrow)

library(monocle3)
library(EGAD)


cell_feature_columns <- readr::read_tsv(
    "product/cell_feature_columns_TS_202008.tsv") %>%
    dplyr::mutate(feature = feature %>% stringr::str_replace("Zernike_", "Zernike__")) %>%
    tidyr::separate(
        col = feature,
        into = c("object", "type", "measurement", "dye", "p1", "p2"),
        sep = "_",
        remove = FALSE,
        fill = "right") %>%
    dplyr::mutate(dye = ifelse(dye == "", NA, dye))


cell_metadata_columns <- readr::read_tsv(
    "product/cell_metadata_columns_TS_202008.tsv")


cell_features <- arrow::read_parquet(
    "product/TS2_2M_Cell_MasterDataTable.parquet") %>%
    dplyr::sample_n(500000)

cell_metadata <- cell_features %>%
    dplyr::select(tidyselect::one_of(cell_metadata_columns$feature))

cell_features <- cell_features %>%
    dplyr::select(tidyselect::one_of(cell_feature_columns$feature))

feature_covar <- EGAD::build_coexp_network(
	exprs = cell_features,
	gene.list = cell_feature_columns$feature)
               


###########################
### Monocle3 wants the data like this
###
###     
###             cell
###        ______________
###        |            |
###  gene  | expression |
###        |    data    | 
###        |            |
###        |------------|
###
### and then UMAP embeds the cells

# so to embed the features
# "gene" <- cell_metadata
# "cell" <- cell_feature_columns


cell_index <- 1:500000
row.names(cell_metadata) <- cell_index
cell_metadata$gene_short_name <- cell_index
cell_features <- cell_features %>% as.matrix()
row.names(cell_features) <- cell_index
colnames(cell_features) <- cell_feature_columns$feature
rownames(cell_feature_columns) <- cell_feature_columns$feature


cds <- monocle3::new_cell_data_set(
    expression_data = cell_features,
    gene_metadata = cell_metadata,
    cell_metadata = cell_feature_columns)

cds <- cds %>% monocle3::preprocess_cds(
    num_dims = 704,
    norm_method = "none",
    verbose = TRUE)

cds <- cds %>% monocle3::reduce_dimension(
    preprocess_method = "PCA",
    a = .6,
    b = .9,
    n_epochs = 10000,
    verbose = TRUE)

cds <- cds %>% monocle3::cluster_cells(
    resolution = 1e-2,
    num_iter = 10,
    verbose = TRUE)

plot <- cds %>% monocle3::plot_cells(
    cell_size = 1,
    alpha = .8,
    color_cells_by = "object",
    label_cell_groups = FALSE) +
    ggplot2::scale_color_discrete("CellProfiler Object") +
    ggplot2::theme_bw() +
    ggplot2::coord_fixed()

ggplot2::ggsave(
    filename = "product/figures/TS2/embed_features/UMAP_500k_by_object_20210125.pdf",
    plot = plot,
    width=7,
    height=7)

plot <- cds %>% monocle3::plot_cells(
    cell_size = 1,
    alpha = .8,
    color_cells_by = "type",
    label_cell_groups = FALSE) +
    ggplot2::scale_color_discrete("CellProfiler Type") +
    ggplot2::theme_bw() +
    ggplot2::coord_fixed()

ggplot2::ggsave(
    filename = "product/figures/TS2/embed_features/UMAP_500k_by_type_20210125.pdf",
    plot = plot,
    width=7,
    height=7)


plot <- cds %>% monocle3::plot_cells(
    cell_size = 1,
    alpha = .8,
    color_cells_by = "measurement",
    label_cell_groups = FALSE) +
    ggplot2::scale_color_discrete("CellProfiler Measurement") +
    ggplot2::theme_bw() +
    ggplot2::coord_fixed()

ggplot2::ggsave(
    filename = "product/figures/TS2/embed_features/UMAP_500k_by_measurement_20210125.pdf",
    plot = plot,
    width=9,
    height=7)


plot <- cds %>% monocle3::plot_cells(
    cell_size = 1,
    alpha = .8,
    color_cells_by = "dye",
    label_cell_groups = FALSE) +
    ggplot2::scale_color_discrete("CellProfiler Dye") +
    ggplot2::theme_bw() +
    ggplot2::coord_fixed()

ggplot2::ggsave(
    filename = "product/figures/TS2/embed_features/UMAP_500k_by_dye_20210125.pdf",
    plot = plot,
    width=7,
    height=7)

