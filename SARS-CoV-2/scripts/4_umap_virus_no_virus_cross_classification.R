library(plyr)
library(tidyverse)
library(fuzzyjoin)
library(ggplot2)
library(readxl)
library(MPStats)
library(arrow)
library(viridis)

library(caret)
library(gbm)

feature_columns <- readr::read_tsv(
    "raw_data/cell_feature_columns.tsv")

cell_features <- arrow::read_parquet(
    "product/top_hit_cells_plate_scaled_200522a.parquet",
    col_select = c("plate_id", "Compound", "dose_nM"))


cluster_membership <- arrow:::read_parquet(
    "~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/top_hits_plate_scaled_200522a_umap2_2M_15_0.0/roi_membership.parquet")


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
        cluster_index = dplyr::case_when(
            cluster_index == 0 ~ 1,
            cluster_index == 1 ~ 2,
            cluster_index == 2 ~ 3,
            cluster_index == 3 ~ 4,
            cluster_index == 4 ~ -1,
            TRUE ~ cluster_index),
        cluster_index = ifelse(cluster_index == -1, Inf, cluster_index),
        cluster_label = paste0("ROI ", cluster_index) %>% factor() %>% reorder(cluster_index)) %>%
    dplyr::select(-value)


cluster_membership_no_virus <- arrow:::read_parquet(
    "~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/top_hits_plate_scaled_200522a_no_virus_umap2_2M_15_0.0/roi_membership.parquet")


cluster_labels_no_virus <- cluster_membership_no_virus %>%
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
        cluster_label = paste0("ROI ", cluster_index) %>% factor() %>% reorder(cluster_index)) %>%
    dplyr::select(-value)


virus_no_virus_cross_clusters <- data.frame(
    cluster_label_with_virus = cluster_labels$cluster_label %>%
        factor(
            level  = c("ROI Inf", "ROI 14", "ROI 4", "ROI 3", "ROI 2", "ROI 1", "ROI 12", "ROI 15", "ROI 10", "ROI 15", "ROI 9", "ROI 8", "ROI 7", "ROI 13", "ROI 11", "ROI 6", "ROI 5"),
            labels = c("ROI Inf", "ROI 14", "ROI 4", "ROI 3", "ROI 2", "ROI 1", "ROI 12", "ROI 15", "ROI 10", "ROI 15", "ROI 9", "ROI 8", "ROI 7", "ROI 13", "ROI 11", "ROI 6", "ROI 5")),
    cluster_label_no_virus = cluster_labels_no_virus$cluster_label) %>%
    dplyr::count(cluster_label_with_virus, cluster_label_no_virus)


# quick counts
z <- virus_no_virus_cross_clusters %>%
    tidyr::pivot_wider(
        id_cols = cluster_label_no_virus,
        names_from = cluster_label_with_virus,
        values_from = "n") %>%
    as.matrix() %>%
    blockcluster::coclusterContinuous(nbcocluster = 10)

p <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.background = ggplot2::element_rect(
            fill = "#00274C", colour = "#00274C"),
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom") +
    ggplot2::geom_tile(
        data = virus_no_virus_cross_clusters,
        mapping = ggplot2::aes(
            x = cluster_label_no_virus,
            y = cluster_label_with_virus,
            fill = log(n))) +
    ggplot2::scale_x_discrete("Clusters with out virus channel") +
    ggplot2::scale_y_discrete(
        "Clusters with virus channel") +
    viridis::scale_fill_viridis(
        "Cell count")

ggplot2::ggsave(
    filename = "product/figures/umap_features/top_hits_plate_scaled_200522a_no_virus_umap2_2M_15_0.0/virus_no_virus_cross_clusters.pdf",
    plot = p)
