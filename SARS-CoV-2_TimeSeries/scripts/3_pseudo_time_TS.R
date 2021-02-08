library(plyr)
library(tidyverse)
library(arrow)
library(monocle3)
library(BiocNeighbors)
library(MPStats)

cell_metadata_columns <- readr::read_tsv("raw_data/cell_metadata_columns_TS.tsv")

cell_feature_columns <- readr::read_tsv("raw_data/cell_feature_columns_TS.tsv")

plate_ids <- readr::read_tsv("raw_data/plate_ids.tsv") %>%
    dplyr::filter(plate_id %>% stringr::str_detect("^TS")) %>%
    magrittr::extract2("plate_id")

########################
# Gather cell features #
########################

# Gather the features from the database
# source("1.2_load_cell_features.R")

# or load them back from S3
plate_ids %>%
    plyr::l_ply(function(plate_id) {
        cat("Copying features for plate '", plate_id, "'\n", sep = "")
        command <- paste0("cp ~/bucket_umich-insitro/CQ1/pseudo_time_202006_cell_profiler_features/covid19cq1_SARS_", plate_id, "_Cell_MasterDataTable.parquet product/")
        cat(command, "\n", sep = "")
        system(command)
    })
# the s3fuse fills up the ~/tmp folder after copying lots of data so clear it out
system("sudo rm -rf ~/tmp/*")

cell_features <- plate_ids %>%
    plyr::ldply(function(plate_id) {
        cat("Loading features for plate '", plate_id, "'\n", sep = "")
        arrow::read_parquet(
            paste0("product/covid19cq1_SARS_", plate_id, "_Cell_MasterDataTable.parquet"),
            col_select = c(
                cell_metadata_columns$feature,
                Nuclei_Number_Object_Number,
                Nuclei_Location_Center_X,
                Nuclei_Location_Center_Y,
                cell_feature_columns$feature)) %>%
            dplyr::mutate_at(cell_feature_columns$feature, ~ scale(.)[, 1])
    }) %>%
    dplyr::mutate(
        plate_id = factor(
            plate_id,
            levels = plate_ids,
            labels = plate_ids %>%
                stringr::str_replace("TS", "")))

cell_features %>%
    arrow::write_parquet(
       "raw_data/covid19cq1_TS_scaled_Cell_MasterDataTable.parquet")

########################
# Embed 2M cell subset #
########################

system("
        cd ~/opt/MPLearn/vignettes/SARS-CoV-2/S25 &&
        /home/ubuntu/anaconda3/envs/sextonlab/bin/python \
            ~/anaconda3/envs/sextonlab/bin/embed_umap \
            --dataset /home/ubuntu/opt/MPStats/vignettes/SARS-CoV-2/raw_data/covid19cq1_TS_scaled_Cell_MasterDataTable.parquet \
            --tag UMAP_embedding_TS_scaled_2M \
            --feature_columns ~/opt/MPStats/vignettes/SARS-CoV-2/raw_data/cell_feature_columns_TS.tsv \
            --no_standardize_features \
            --random_subset 2000000 \
	    --umap_low_memory \
	    --verbose
")

######################
# Re-embed all cells #
######################

system("
        cd ~/opt/MPLearn/vignettes/SARS-CoV-2/S25 &&
        /home/ubuntu/anaconda3/envs/sextonlab/bin/python \
            ~/anaconda3/envs/sextonlab/bin/embed_umap \
            --ref_embed_dir intermediate_data/UMAP_embedding_TS_scaled_2M \
            --dataset raw_data/covid19cq1_TS_scaled_Cell_MasterDataTable.parquet \
            --tag UMAP_embedding_TS_scaled_full \
            --feature_columns ~/opt/MPStats/vignettes/SARS-CoV-2/raw_data/cell_feature_columns_TS.tsv \
            --no_standardize_features \
	    --verbose
")


#############################
# Identify infected cluster #
#############################

# identify infected cells using
# MPLearn/vignettes/SARS-CoV-2/notebooks/pseudo_time_202006.ipynb

cell_features <- arrow::read_parquet(
    "raw_data/covid19cq1_TS_scaled_Cell_MasterDataTable.parquet")
# 5090202

infected_ROI_membership <- arrow::read_parquet(
    "~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/UMAP_embedding_TS_scaled_full/infected_ROI_membreship.parquet")
# 5090202

infected_cells <- cell_features %>%
    dplyr::bind_cols(
        infected_ROI_membership %>%
        dplyr::mutate(is_infected = roi_0)) %>%
    dplyr::filter(is_infected)

infected_cells %>%
    arrow::write_parquet(
        "~/opt/MPStats/vignettes/SARS-CoV-2/raw_data/covid19cq1_TS_scaled_infected_Cell_MasterDataTable.parquet")

########################
# Embed infected cells #
########################

system("
        cd ~/opt/MPLearn/vignettes/SARS-CoV-2/S25 &&
        /home/ubuntu/anaconda3/envs/sextonlab/bin/python \\
            ~/anaconda3/envs/sextonlab/bin/embed_umap \\
            --dataset /home/ubuntu/opt/MPStats/vignettes/SARS-CoV-2/raw_data/covid19cq1_TS_scaled_infected_Cell_MasterDataTable.parquet \\
            --tag UMAP_embedding_TS_scaled_infected \\
            --feature_columns ~/opt/MPStats/vignettes/SARS-CoV-2/raw_data/cell_feature_columns_TS.tsv \\
            --no_standardize_features \\
	    --umap_low_memory \\
	    --verbose
")

###################################################
# Explore embedding parameters for infected cells #
###################################################

umap_params <- expand.grid(
    umap_a = 1,
    umap_b = c(.5),
    umap_n_neighbors = c(500),
    umap_n_epochs = 2000,
    umap_negative_sample_rate = c(20))


umap_params %>% plyr::a_ply(1, function(params) {
    dataset_tag <- paste0("UMAP_embedding_TS_scaled_infected_a=", params$umap_a[1], "_b=", params$umap_b[1], "_n_neighbors=", params$umap_n_neighbors[1], "_negative_sample_rate=", params$umap_negative_sample_rate[1], "_n_epochs=", params$umap_n_epochs[1])
    if(dir.exists(
        paste0("~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/", dataset_tag))){
        cat("Embedding '", dataset_tag, "' exists, skipping...\n")
        return(NULL)
    }
    command <- paste0("
        cd ~/opt/MPLearn/vignettes/SARS-CoV-2/S25 &&
        /home/ubuntu/anaconda3/envs/sextonlab/bin/python \\
            ~/anaconda3/envs/sextonlab/bin/embed_umap \\
            --dataset /home/ubuntu/opt/MPStats/vignettes/SARS-CoV-2/raw_data/covid19cq1_TS_scaled_infected_Cell_MasterDataTable.parquet \\
            --tag ", dataset_tag, " \\
            --feature_columns ~/opt/MPStats/vignettes/SARS-CoV-2/raw_data/cell_feature_columns_TS.tsv \\
            --umap_a ", params$umap_a[1], " \\
            --umap_b ", params$umap_b[1], " \\
            --umap_n_neighbors ", params$umap_n_neighbors[1], " \\
            --umap_negative_sample_rate ", params$umap_negative_sample_rate[1], " \\
            --no_standardize_features \\
            --umap_low_memory \\
            --verbose")
    cat(command, "\n", sep = "")
    system(command)

    cat("Loading embedding ...\n")
    infected_embedding <- arrow::read_parquet(
        paste0("~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/", dataset_tag, "/umap_embedding.parquet"))

    infected_cells_and_embedding <- infected_cells %>%
        dplyr::bind_cols(infected_embedding)
    
    plot <- ggplot2::ggplot(
        data = infected_cells_and_embedding) +
        theme_bw() +
        ggplot2::geom_point(
            mapping = ggplot2::aes(
                x = UMAP_1,
                y = UMAP_2),
            size = .2,
            alpha = 1) +
        ggplot2::coord_fixed() +
        ggplot2::facet_wrap(
            facets = vars(plate_id))

    ggplot2::ggsave(
        filename = paste0("product/figures/pseudo_time_202006/", dataset_tag, "_plate_id_facets_200722.pdf"),
        plot = plot,
        width = 15,
        height = 10)
    
    
    })



runs <- list.files(
    path = "/home/ubuntu/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data",
    pattern = "UMAP_embedding_TS_scaled_infected_") %>%
    tibble::tibble(embedding_tag = .) %>%
    dplyr::mutate(
        a = stringr::str_match(embedding_tag, "a=([^_]+)")[, 2] %>%
            as.numeric(),
        b = stringr::str_match(embedding_tag, "b=([^_]+)")[, 2] %>%
            as.numeric(),
        n_neighbors = stringr::str_match(embedding_tag, "n_neighbors=([^_]+)")[, 2] %>%
            as.numeric(),
        negative_sample_rate = stringr::str_match(embedding_tag, "negative_sample_rate=([^_]+)")[, 2] %>%
            as.numeric(),
        negative_sample_rate = dplyr::case_when(
            is.na(negative_sample_rate) ~ 5,
            TRUE ~ negative_sample_rate),
        n_epochs = stringr::str_match(embedding_tag, "n_epochs=([^_]+)")[, 2] %>%
            as.numeric(),
        n_epochs = dplyr::case_when(
            is.na(n_epochs) ~ 200,
            TRUE ~ n_epochs))

# evaluate how clustering 'clumpiness' changes with the parameters"
# where clumpiness is the expected L
source("scripts/embedding_statistics.R")
runs <- runs %>% plyr::adply(1, function(run) {
    cat("Collecting statistics for embedding '", run$embedding_tag[1], "' ...\n", sep = "")
    embedding_fname <- paste0("~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/", run$embedding_tag[1], "/umap_embedding.parquet")
    if (!file.exists(embedding_fname)) {
        cat("   Embedding doesn't exist, skipping ...\n")
        return(NULL)
    }
    embedding <- arrow::read_parquet(embedding_fname)
    point_process <- populate_point_process(
        data = embedding,
        coordinate_columns = c("UMAP_1", "UMAP_2"))
    expected_l <- clustering_density(point_process = point_process)
    data.frame(expected_l = expected_l)
})


plot <- runs %>%
    tidyr::pivot_longer(
        cols = c(a, b, n_neighbors, negative_sample_rate, n_epochs),
        names_to = "parameter") %>%
    ggplot2::ggplot(data = .) +
    ggplot2::geom_point(
        mapping = ggplot2::aes(
            x = value,
            y = expected_l)) +
    ggplot2::facet_wrap(
        facets = dplyr::vars(parameter),
        scales = "free_x")

ggplot2::ggsave(
    filename = "product/figures/pseudo_time_202006/infected_UMAP_expected_l_by_params.pdf",
    width = 6,
    height = 4)



##########################################
# Plot compound by time-point embedding #
##########################################

infected_embedding <- arrow::read_parquet(
    "~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/UMAP_embedding_TS_scaled_infected_a=1_b=0.5_n_neighbors=500_negative_sample_rate=20_n_epochs=2000/umap_embedding.parquet")

infected_cells_and_embedding <- infected_cells %>%
    dplyr::bind_cols(infected_embedding)

plot <- ggplot2::ggplot(
    data = infected_cells_and_embedding) +
    theme_bw() +
    ggplot2::geom_point(
        mapping = ggplot2::aes(
            x = UMAP_1,
            y = UMAP_2),
        size = .01,
        shape = 16,
        alpha = .4) +
    ggplot2::coord_fixed() +
    ggplot2::facet_grid(
        rows = vars(Compound),
        cols = vars(plate_id))

ggplot2::ggsave(
    filename = "product/figures/pseudo_time_202006/UMAP_embedding_TS_scaled_infected_a=1_b=0.5_n_neighbors=500_negative_sample_rate=20_n_epochs=2000_embedding_plate_id_Compound_facets_200722.pdf",
    plot = plot,
    width = 15,
    height = 15)

####################################################
# Use Monocle3 to cluster and estimate pseudo-time #
####################################################

infected_cds <- MPStats::populate_cds(
    cell_features = infected_cells_and_embedding,
    cell_feature_columns = cell_feature_columns,
    cell_metadata_columns = cell_metadata_columns,
    embedding_type = c("UMAP"),
    embedding = infected_cells_and_embedding %>% dplyr::select(UMAP_1, UMAP_2),
    verbose = TRUE)

# as resolution gets bigger --> more clusters
infected_cds <- infected_cds %>%
    monocle3::cluster_cells(
        resolution = .0000000001,
        num_iter = 10,
        verbose = TRUE)

infected_cds %>% MPStats::serialize_clusters(
    output_fname = "~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/UMAP_embedding_TS_scaled_infected_a=1_b=0.5_n_neighbors=500_negative_sample_rate=20_n_epochs=2000/clusters_leiden_res=1e-10.parquet")

infected_cds <- infected_cds %>%
    monocle3::learn_graph(
        close_loop = FALSE,
        verbose = TRUE)

infected_cds <- infected_cds %>%
    monocle3::order_cells(verbose = TRUE)

save(infected_cds, file="~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/UMAP_embedding_TS_scaled_infected/infected_cds.Rdata")

plot <- infected_cds %>%
    monocle3::plot_cells(
        color_cells_by = "cluster",
        label_cell_groups = FALSE,
        label_leaves = FALSE,
        label_branch_points = FALSE,
        cell_size = .2,
        alpha = 1) +
    ggplot2::coord_fixed() +
    ggplot2::facet_wrap(
        facets = vars(plate_id))

ggplot2::ggsave(
    filename = "product/figures/pseudo_time_202006/UMAP_embedding_TS_scaled_infected_a=1_b=0.5_n_neighbors=500_negative_sample_rate=20_n_epochs=2000_trajectory_plate_id_facets_200722.pdf",
    plot = plot,
    width = 15,
    height = 10)

#####################################################
# plot clusters distribution by drug and time-point #
#####################################################

infected_cells_and_embedding <- infected_cells_and_embedding %>%
    dplyr::bind_cols(
        arrow::read_parquet("~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/UMAP_embedding_TS_scaled_infected_a=1_b=0.5_n_neighbors=500_negative_sample_rate=20_n_epochs=2000/clusters_leiden_res=1e-10.parquet"))

infected_cells_and_embedding <- infected_cells_and_embedding %>%
    dplyr::mutate(
        time_point_hours = Time_Point %>%
            stringr::str_extract("^[0-9]+") %>%
            as.numeric(),
        cluster_label = cluster_label %>%
            factor(
                levels = c(11, 14, 13, 12, 8, 6, 9, 4, 7, 5, 2, 1, 10, 3, 15),
                labels = c(11, 14, 13, 12, 8, 6, 9, 4, 7, 5, 2, 1, 10, 3, 15)))

data <- infected_cells_and_embedding %>%
    dplyr::group_by(cluster_label) %>%
        dplyr::mutate(cluster_count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(plate_id, cluster_label) %>%
        dplyr::mutate(plate_cluster_count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(
        time_point_hours,
        plate_id,
        cluster_label,
        cluster_count,
        plate_cluster_count)

plot <- ggplot2::ggplot(data = data) +
    ggplot2::theme_bw() +
    ggplot2::geom_line(
        mapping = ggplot2::aes(
            x = time_point_hours,
            y = plate_cluster_count / cluster_count,
            color = cluster_label),
        size = .8) +
    ggplot2::geom_point(
        mapping = ggplot2::aes(
            x = time_point_hours,
            y = plate_cluster_count / cluster_count,
            color = cluster_label),
        size = 2) +
    ggplot2::scale_x_continuous("Time Point (h)") +
    ggplot2::scale_y_continuous("Fraction of Cells in Cluster")

ggplot2::ggsave(
    filename = "product/figures/pseudo_time_202006/UMAP_embedding_TS_scaled_infected_a=1_b=0.5_n_neighbors=500_negative_sample_rate=20_n_epochs=2000/clusters_by_time_point.pdf",
    plot = plot,
    width = 6,
    height = 6)

plot <- ggplot2::ggplot(data = data) +
    ggplot2::theme_bw() +
    ggplot2::geom_line(
        mapping = ggplot2::aes(
            x = time_point_hours,
            y = plate_cluster_count / cluster_count),
        size = .8) +
    ggplot2::geom_point(
        mapping = ggplot2::aes(
            x = time_point_hours,
            y = plate_cluster_count / cluster_count),
        size = 2) +
    ggplot2::facet_wrap(facets = dplyr::vars(cluster_label)) +
    ggplot2::scale_x_continuous("Time Point (h)") +
    ggplot2::scale_y_continuous("Fraction of Cells in Cluster")

ggplot2::ggsave(
    filename = "product/figures/pseudo_time_202006/UMAP_embedding_TS_scaled_infected_a=1_b=0.5_n_neighbors=500_negative_sample_rate=20_n_epochs=2000/clusters_by_time_point_facets.pdf",
    plot = plot,
    width = 6,
    height = 6)


data <- infected_cells_and_embedding %>%
    dplyr::group_by(cluster_label) %>%
        dplyr::mutate(cluster_count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Compound, cluster_label) %>%
        dplyr::mutate(Compound_cluster_count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(
        Compound,
        cluster_label,
        cluster_count,
        Compound_cluster_count)

plot <- ggplot2::ggplot(data = data) +
    ggplot2::theme_bw() +
    ggplot2::geom_bar(
        mapping = ggplot2::aes(
            x = cluster_label,
            y = Compound_cluster_count / cluster_count),
        stat = "identity") +
    ggplot2::facet_wrap(
        facets = dplyr::vars(Compound),
        scales = "free_y") +
    ggplot2::scale_x_discrete("Cluster") +
    ggplot2::scale_y_continuous("Fraction of Cells in Cluster")

ggplot2::ggsave(
    filename = "product/figures/pseudo_time_202006/UMAP_embedding_TS_scaled_infected_a=1_b=0.5_n_neighbors=500_negative_sample_rate=20_n_epochs=2000/clusters_by_Compound_facets.pdf",
    plot = plot,
    width = 9,
    height = 6)


data <- infected_cells_and_embedding %>%
    dplyr::group_by(cluster_label) %>%
        dplyr::mutate(cluster_count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Compound, plate_id, cluster_label) %>%
        dplyr::mutate(Compound_time_point_cluster_count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(
        Compound,
        time_point_hours,
        plate_id,
        cluster_label,
        cluster_count,
        Compound_time_point_cluster_count)

plot <- ggplot2::ggplot(data = data) +
    ggplot2::theme_bw() +
    ggplot2::geom_line(
        mapping = ggplot2::aes(
            x = cluster_label,
            y = Compound_time_point_cluster_count / cluster_count,
            group = plate_id,
            color = plate_id)) +
    ggplot2::facet_wrap(
        facets = dplyr::vars(Compound),
        scales = "free_y") +
    ggplot2::scale_x_discrete("Cluster") +
    ggplot2::scale_y_continuous("Fraction of Cells in Cluster")

ggplot2::ggsave(
    filename = "product/figures/pseudo_time_202006/UMAP_embedding_TS_scaled_infected_a=1_b=0.5_n_neighbors=500_negative_sample_rate=20_n_epochs=2000/clusters_by_Compound_time_point_facets.pdf",
    plot = plot,
    width = 9,
    height = 6)




##########################################
# PLS enriched features for each cluster #
##########################################

feature_importance <- data.frame(cluster_label = 1:15) %>%
    plyr::adply(1, function(df) {
        cat("Computing feature importance for cluster '", df$cluster_label[1], "'\n", sep = "")
        data_matrix <- infected_cells_and_embedding  %>%
            dplyr::mutate(label = dplyr::case_when(
                as.character(cluster_label) != df$cluster_label[1] ~ "other",
                as.character(cluster_label) == df$cluster_label[1] ~ "cluster",
                TRUE ~ NA_character_) %>%
                    factor(
                        levels = c("cluster", "other"),
                        labels = c("cluster", "other"))) %>%
            dplyr::select(
                label,
                tidyselect::one_of(cell_feature_columns$feature)) %>%
            tidyr::drop_na()
        in_training <- caret::createDataPartition(
            data_matrix$label,
            p = .80,
            list = FALSE)
        data_train <- data_matrix %>% dplyr::slice(in_training)
        data_test  <- data_matrix %>% dplyr::slice(-in_training)
        cat("    N train: ", nrow(data_train), "\n", sep = "")
        cat("    N test:  ", nrow(data_test), "\n", sep = "")
        fit_control <- caret::trainControl(## 10-fold CV
            method = "repeatedcv",
            number = 2,
            classProbs = TRUE)
        pls_fit <- caret::train(
            label ~ .,
            data = data_train,
            method = "pls",
            trControl = fit_control,
            preProc = c("center", "scale"),
            tuneLength = 3,
            verbose = TRUE)
        pred_fit <- predict(pls_fit, newdata = data_test)
        feature_importance <- pls_fit %>% caret::varImp()
        feature_importance$importance %>%
            dplyr::bind_cols(
                cell_feature_columns %>%
                dplyr::select(feature)) %>%
            dplyr::rename(feature_score = Overall) %>%
            dplyr::arrange(desc(feature_score)) %>%
            head(20) %>%
            dplyr::mutate(cluster_label = df$cluster_label[1])
    })

p <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_point(
        data = feature_importance %>%
            dplyr::mutate(
                cluster_index = cluster_label %>%
                    stringr::str_extract("[0-9]+$") %>%
                    as.numeric(),
                cluster_label = factor(cluster_label) %>%
                    reorder(cluster_index)),
        mapping = ggplot2::aes(
            x = feature_score,
            y = feature)) +
    ggplot2::scale_x_continuous("PLS feature importance") +
    ggplot2::scale_y_discrete(
        "Feature",
        drop = FALSE) +
    ggplot2::facet_wrap(
        ~cluster_label,
        ncol = 3,
        scales = "free_y")



ggplot2::ggsave(
    filename = "product/figures/pseudo_time_202006/UMAP_embedding_TS_scaled_infected_a=1_b=0.5_n_neighbors=500_negative_sample_rate=20_n_epochs=2000/clusters_feature_importance.pdf",
    plot = p,
    width = 15,
    height = 25)

feature_importance %>%
    readr::write_tsv("product/figures/pseudo_time_202006/UMAP_embedding_TS_scaled_infected_a=1_b=0.5_n_neighbors=500_negative_sample_rate=20_n_epochs=2000/clusters_feature_importance.tsv")





##########################################
# Monocle3 based representative features #
##########################################

marker_test_res <- monocle3::top_markers(
    infected_cds,
    cores = 30,
    verbose = TRUE)

top_specific_markers <- marker_test_res %>%
    dplyr::group_by(cell_group) %>%
    dplyr::top_n(1, pseudo_R2)

top_specific_marker_ids <- top_specific_markers %>%
    dplyr::pull(gene_id) %>%
    unique()

plot <- monocle3::plot_genes_by_group(
    infected_cds,
    top_specific_marker_ids,
    ordering_type = "maximal_on_diag",
    max.size = 3)

ggplot2::ggsave(
    filename = "product/figures/pseudo_time_202006/UMAP_embedding_TS_scaled_infected_a=1_b=0.5_n_neighbors=500_negative_sample_rate=20_n_epochs=2000_top_markers_200722.pdf",
    plot = plot,
    width = 12,
    height = 12)
