

library(plyr)
library(tidyverse)
library(arrow)


plate_ids <- readr::read_tsv("raw_data/plate_ids.tsv")
feature_columns <- readr::read_tsv(
    "raw_data/cell_feature_columns.tsv")

compounds_of_interest <- c(
    "DGD1202",
    "MCTI253",
    "MCTI178",
    "CIRC825")

cell_features <- plate_ids %>%
    dplyr::filter(
        master_plate_id == "1005",
        plate_id %>% stringr::str_detect("C$")) %>%
    plyr::adply(1, function(plate){
        schema <- plate$schema[1]
        plate_id <- plate$plate_id[1]
        input_fname <- paste0(
            "product/", schema, "_SARS_", plate_id, "_Cell_MasterDataTable.parquet")
        cat("Loading cell_features from '", input_fname, "' ...\n", sep = "")
        arrow::read_parquet(input_fname)
    })


cell_features_sample <- cell_features %>%
    dplyr::filter(
        Compound %in% c(
            compounds_of_interest,
            "Positive Control",
            "Negative Control")) %>%
    cell_features %>% dplyr::filter(Compound %in% compounds_of_interest),
    cell_features %>%
    dplyr::filter(Compound %in% "Negative Control") %>%
    dplyr::sample_n(1500000),
    cell_features %>%
    dplyr::filter(Compound %in% "Positive Control") %>%
    dplyr::sample_n(250000))

cell_features_sample %>%
    arrow::write_parquet(
        "raw_data/covid19primary_1005_hits_Cell_MasterDataTable.parquet")


######
#library(uwot)
#
#cell_features <- arrow::read_parquet(
#    "raw_data/covid19primary_1005_hits_Cell_MasterDataTable.parquet")
#
#    
#
#embedding <- uwot::umap(
#    cell_features,
#    n_components = 2,
#    metric = "euclidean",
#    min_dist = 0,
#    n_neighbors = 15,
#    fast_sgd = TRUE,
#    n_threads = 30,
#    verbose = TRUE,
#    nn_method = "annoy",
#    pca = ncol(cell_features),
#    ret_model = TRUE)

system("
        cd ~/opt/MPLearn/vignettes/SARS-CoV-2/S25 &&
        ~/anaconda3/envs/sextonlab/bin/embed_umap \
            --dataset /home/ubuntu/opt/MPStats/vignettes/SARS-CoV-2/raw_data/covid19primary_1005_hits_Cell_MasterDataTable.parquet \
            --tag UMAP_embedding_1005_hits_full \
            --feature_columns ~/opt/MPStats/vignettes/SARS-CoV-2/raw_data/cell_feature_columns.tsv \
	    --umap_low_memory \
	    --verbose
")


cell_metadata <- cell_features_sample %>%
    dplyr::select(
        -tidyselect::one_of(feature_columns$feature))
        
expression_data <- cell_features_sample %>%
    dplyr::select(
        tidyselect::one_of(feature_columns$feature)) %>%
    as.matrix() %>%
    t()

row.names(expression_data) <- feature_columns$feature
row.names(feature_columns) <- feature_columns$feature
names(expression_data) <- 1:ncol(expression_data)
row.names(cell_metadata) <- 1:ncol(expression_data)


cds <- monocle3::new_cell_data_set(
    expression_data = expression_data,
    cell_metadata = cell_metadata,
    gene_metadata = feature_columns)

reducedDims(cds)[['UMAP']] <- arrow::read_parquet(
    "~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/UMAP_embedding_1005_hits_full/umap_embedding.parquet")

row.names(colData(cds)) <- 1:ncol(expression_data)
cds <- cds %>% monocle3::cluster_cells(
    resolution = 1e-8,
    num_iter = 10,
    verbose = TRUE)


cds %>% monocle3::plot_cells(
    reduction_method = "UMAP",
    cell_size = 0.1,
    color_cells_by="cluster",
    alpha=.8) +
    #ggplot2::scale_color_manual(values="black") +
    ggplot2::facet_grid(
        rows = dplyr::vars(Compound),
        cols = dplyr::vars(dose_nM)) +
    ggplot2::theme(
        panel.border =
            ggplot2::element_rect(color = "black", fill = NA, size = .67),
        panel.grid.major.x =
            ggplot2::element_line(size = 0.05, color = "grey20"),
#        panel.grid.minor.x =
#            ggplot2::element_line(size = 0.05, color = "grey50"),
        panel.grid.major.y =
            ggplot2::element_line(size = 0.05, color = "grey20"))
#        panel.grid.minor.y =
#            ggplot2::element_line(size = 0.25, color = "grey50"))
ggplot2::ggsave(
    filename = "product/figures/plate_1005_hits/umap_embedding_compound_by_dose_200709.png",
    height = 10,
    width = 10)
    


