

library(plyr)
library(tidyverse)
library(arrow)
library(monocle3)

cell_metadata_columns <- readr::read_tsv("raw_data/cell_metadata_columns_TS.tsv")

cell_feature_columns <- readr::read_tsv("raw_data/cell_feature_columns_TS.tsv")

plate_ids <- readr::read_tsv("raw_data/plate_ids.tsv") %>%
    dplyr::filter(plate_id %>% stringr::str_detect("^TS")) %>%
    magrittr::extract2("plate_id")


cell_metadata_dataframe <- plate_ids %>%
    plyr::ldply(function(plate_id){
        cat("Loading metadata for plate '", plate_id, "'\n", sep = "")
        arrow::read_parquet(
            paste0("product/covid19cq1_SARS_", plate_id, "_Cell_MasterDataTable.parquet"),
            col_select = c(
                cell_metadata_columns$feature,
                Nuclei_Number_Object_Number,
                Nuclei_Location_Center_X,
                Nuclei_Location_Center_Y))
    }) %>%
    dplyr::mutate(
        plate_id = factor(
            plate_id,
            levels=plate_ids,
            labels=plate_ids %>%
                stringr::str_replace("TS", "")))

cell_metadata_dataframe %>%
    dplyr::count(plate_id) %>%
    dplyr::arrange(plate_id)

# load cell features as a matrix with dimensions [feature, cell]
cell_features_dataframe <- plate_ids %>%
    plyr::ldply(function(plate_id){
        cat("Loading features for plate '", plate_id, "'\n", sep = "")
        arrow::read_parquet(
            paste0("product/covid19cq1_SARS_", plate_id, "_Cell_MasterDataTable.parquet"),
            col_select = cell_feature_columns$feature)
    })

cell_sample <- cell_metadata_dataframe %>%
    nrow() %>%
    sample(1000000)

cell_metadata <- cell_metadata_dataframe %>%
    dplyr::slice(cell_sample)





# shift features to be non-negative to prevent nans in log-normalization
cell_features <- cell_features_dataframe %>%
    dplyr::slice(cell_sample) %>%
    dplyr::mutate_at(vars(tidyselect::matches("_EulerNumber$")), ~ . + 2) %>%
    dplyr::mutate_at(vars(tidyselect::matches("_Orientation$")), ~ . + 90)

cell_features <- cell_features %>%
    as.matrix() %>%
    t()

row.names(cell_features) <- cell_feature_columns$feature
row.names(cell_feature_columns) <- cell_feature_columns$feature
names(cell_features) <- 1:ncol(cell_features)
row.names(cell_metadata) <- 1:ncol(cell_features)

cell_data_set <- monocle3::new_cell_data_set(
    expression_data = cell_features,
    cell_metadata = cell_metadata,
    gene_metadata = cell_feature_columns)

# not sure where this is suppose to be set, but it is stripped in new_cell_data_set
# and needed in cluster_cells
row.names(colData(cell_data_set)) <- 1:ncol(cell_features)

cell_data_set <- cell_data_set %>% monocle3::preprocess_cds(
    num_dims = 100,
    verbose = TRUE)

#cell_data_set <- cell_data_set %>% monocle3::align_cds(cell_data_set, preprocess_method="PCA")

cell_data_set <- cell_data_set %>% monocle3::reduce_dimension(
    preprocess_method = "PCA",
    cores = 30,
    umap.min_dist = 0,
    umap.n_neighbors = 30,
    verbose = TRUE)

cell_data_set <- cell_data_set %>% monocle3::cluster_cells(
    resolution = 1e-8,
    num_iter = 10,
    verbose = TRUE)

cell_data_set <- cell_data_set %>% monocle3::learn_graph(
    verbose = TRUE)

cell_data_set <- cell_data_set %>% monocle3::order_cells(
    verbose = TRUE)


plot <- cell_data_set %>%
    monocle3::plot_cells(
        color_cells_by = "plate_id",
        label_cell_groups = FALSE,
        cell_size = .05,
        alpha = 2) +
    ggplot2::facet_wrap(~plate_id)
        

ggplot2::ggsave(
    filename = "product/figures/TS/monocle3/cell_1M_embedded_min_dist0_n_neighbors30_Plate_Name_facets_200708.png",
    plot = plot,
    width=10,
    height=10)


plot <- cell_data_set %>%
    monocle3::plot_cells(
        color_cells_by = "Compound",
        label_cell_groups = FALSE,
        cell_size = .05,
        alpha = 2) +
    ggplot2::facet_grid(
        rows = vars(Compound),
        cols = vars(plate_id))
        

ggplot2::ggsave(
    filename = "product/figures/TS/monocle3/cell_1M_embedded_min_dist0_n_neighbors30_plate_id_Compound_facets_200708.png",
    plot = plot,
    width = 15,
    height = 15)



ggplot2::ggsave(
    filename = "product/figures/TS/monocle3/cell_1M_embedded_min_dist0_n_neighbors50_Plate_Name_200708.pdf",
    plot = plot)



ggplot2::ggsave(
    filename = "product/figures/TS/monocle3/vanilla_10k_n_neighbors50_trajectory_200708.pdf",
    plot = plot)

