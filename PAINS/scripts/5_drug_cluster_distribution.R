


library(plyr)
library(tidyverse)
library(MPStats)
library(arrow)

dataset_ids <- readr::read_tsv("raw_data/dataset_ids.tsv")

cell_feature_columns <- readr::read_tsv("raw_data/cell_feature_columns.tsv")
cell_metadata_columns <- readr::read_tsv("raw_data/cell_metadata_columns.tsv")


embedding <- dataset_ids %>%
    dplyr::select(dataset_id) %>%
    plyr::adply(1, function(dataset) {
        dataset_id <- dataset$dataset_id[1]
        embedding_path <- paste0("intermediate_data/per_plate_into_48h_2M/cell_metadata/", dataset_id, "_cell_metadata.parquet")
        arrow::read_parquet(embedding_path)
    }) 

cluster_membership <- arrow::read_parquet(
    "intermediate_data/per_plate_into_48h_2M/regions_of_interest_membership.parquet")

cluster_labels <- cluster_membership %>%
    dplyr::mutate(cell_index = dplyr::row_number()) %>%
    tidyr::pivot_longer(
        cols = tidyselect::starts_with("roi_"),
        names_to  = "cluster_label") %>%
    dplyr::arrange(desc(value)) %>%
    dplyr::distinct(cell_index, .keep_all = TRUE) %>%
    dplyr::arrange(cell_index) %>%
    dplyr::mutate(
        cluster_label =
            ifelse(!value, "ROI -1", cluster_label) %>%
            stringr::str_replace("roi_", "ROI "),
        cluster_index = cluster_label %>%
            stringr::str_replace("ROI ", "") %>%
            as.numeric(),
        cluster_index = ifelse(cluster_index == -1, Inf, cluster_index),
        cluster_label = factor(cluster_label) %>% reorder(cluster_index)) %>%
    dplyr::select(-value)


embedding <- dplyr::bind_cols(embedding, cluster_labels)

###########################

# drug class by cluster normalized by cluster
drug_class_by_cluster <- embedding %>%
    dplyr::count(class, cluster_label, cluster_index) %>%
    dplyr::group_by(cluster_label) %>%
    dplyr::mutate(n_normalized = n / sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(class), class != "")

plot <- ggplot2::ggplot(
    data = drug_class_by_cluster) +
    ggplot2::theme_bw() +
    ggplot2::theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey40", colour = "grey40")) +
    ggplot2::geom_tile(
        mapping = ggplot2::aes(
            x = cluster_label,
            y = class,
            fill = n_normalized)) +
    ggplot2::ggtitle(
        label = "Drug class by UMAP cluster",
        subtitle = "Normalized by cluster") +
    ggplot2::scale_x_discrete("UMAP Cluster") +
    ggplot2::scale_y_discrete("Drug Class")

ggplot2::ggsave(
    filename = paste0("product/figures/per_plate_into_48h_2M/drug_class_by_cluster_normalized_by_cluster_", MPStats::date_code(), ".pdf"),
    plot = plot,
    width = 8, height = 4)

ggplot2::ggsave(
    filename = paste0("product/figures/per_plate_into_48h_2M/drug_class_by_cluster_normalized_by_cluster_", MPStats::date_code(), ".png"),
    plot = plot,
    width = 8, height = 4)

################################################
# drug class by cluster normalized by cluster
drug_class_by_cluster <- embedding %>%
    dplyr::count(class, cluster_label, time_point, cluster_index) %>%
    dplyr::group_by(class) %>%
    dplyr::mutate(n_normalized = n / sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(class), class != "")

plot <- ggplot2::ggplot(
    data = drug_class_by_cluster) +
    ggplot2::theme_bw() +
    ggplot2::theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey40", colour = "grey40")) +
    ggplot2::geom_tile(
        mapping = ggplot2::aes(
            x = cluster_label,
            y = class,
            fill = log(n_normalized))) +
    ggplot2::facet_wrap(facets = dplyr::vars(time_point)) +
    ggplot2::scale_fill_continuous("log(class normalized cell count)") +
    ggplot2::ggtitle(
        label = "Drug class by UMAP cluster",
        subtitle = "Normalized by cluster") +
    ggplot2::scale_x_discrete("UMAP Cluster") +
    ggplot2::scale_y_discrete("Drug Class")

ggplot2::ggsave(
    filename = paste0("product/figures/per_plate_into_48h_2M/drug_class_by_cluster_by_time_point_normalized_by_drug_class_", MPStats::date_code(), ".pdf"),
    plot = plot,
    width = 10, height = 4)

ggplot2::ggsave(
    filename = paste0("product/figures/per_plate_into_48h_2M/drug_class_by_cluster_by_time_point_normalized_by_drug_class_", MPStats::date_code(), ".png"),
    plot = plot,
    width = 10, height = 4)


linear_fit <- lm(
    formula = log(n) ~ class + cluster_label,
    data= drug_class_by_cluster)
    
drug_class_by_cluster <- drug_class_by_cluster %>%
    dplyr::mutate(
        baseline_log_n = linear_fit %>% predict(),
        baseline_residual = linear_fit$residuals)


plot <- ggplot2::ggplot(
    data = drug_class_by_cluster) +
    ggplot2::theme_bw() +
    ggplot2::theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey40", colour = "grey40")) +
    ggplot2::geom_tile(
        mapping = ggplot2::aes(
            x = cluster_label,
            y = class,
            fill = baseline_residual)) +
    ggplot2::facet_wrap(facets = dplyr::vars(time_point)) +
    ggplot2::scale_fill_continuous("Linear log(n) residual") +
    ggplot2::ggtitle(
        label = "Drug class by UMAP cluster",
        subtitle = "Normalized by cluster") +
    ggplot2::scale_x_discrete("UMAP Cluster") +
    ggplot2::scale_y_discrete("Drug Class")

ggplot2::ggsave(
    filename = paste0("product/figures/per_plate_into_48h_2M/drug_class_by_cluster_by_time_point_linear_residual_", MPStats::date_code(), ".pdf"),
    plot = plot,
    width = 10, height = 4)

ggplot2::ggsave(
    filename = paste0("product/figures/per_plate_into_48h_2M/drug_class_by_cluster_by_time_point_linear_residual_", MPStats::date_code(), ".png"),
    plot = plot,
    width = 10, height = 4)
