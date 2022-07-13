
library(tidyverse)
library(arrow)
library(seriation)

cell_feature_columns <- readr::read_tsv(
    "product/cell_feature_columns_TS_202008.tsv",
    show_col_types = FALSE)

cell_metadata_columns <- readr::read_tsv(
    "product/cell_metadata_columns_TS_202008.tsv",
    show_col_types = FALSE)



feature_list <- readr::read_csv(
    "/home/maom/turbo/weak_supervision/pseudo_time_202008_cell_profiler_features/feature_list.txt",
    show_col_types = FALSE)

df_att <- readr::read_csv(
    "/home/maom/turbo/weak_supervision/pseudo_time_202008_cell_profiler_features/df_att.csv",
    show_col_types = FALSE)


k <- 100
resolution <- 1e-4
num_iter <- 10
cell_clusters <- dplyr::bind_cols(
    arrow::read_parquet(
        file = "product/TS2_plate_scaled_Cell_MasterDataTable.parquet",
        col_select = c(
            cell_metadata_columns$feature,
            "Cells_Number_Object_Number")),
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
        cluster_label = paste("cluster", cluster_label) %>% as.factor() %>% forcats::fct_inorder())

cell_att_cluster <- df_att %>%
    dplyr::select(
        plate_id, row, column, Image_Metadata_FieldID, Cells_Number_Object_Number,
        bag, att_score, normalize_att_score) %>%
    dplyr::left_join(
        cell_clusters %>%
            dplyr::select(
                plate_id, time_point, row, column, Image_Metadata_FieldID, Cells_Number_Object_Number,
                cluster_code, cluster_label),
        by = c("plate_id", "row", "column", "Image_Metadata_FieldID", "Cells_Number_Object_Number"))

#############
plot_data <- cell_att_cluster %>%
    dplyr::group_by(time_point, cluster_label) %>%
    dplyr::summarize(
        mean_normalized_attention_score = mean(normalize_att_score),
        .groups = "drop")

plot_data <- plot_data %>%
    dplyr::inner_join(
        plot_data %>%
            dplyr::group_by(cluster_label) %>%
            dplyr::summarize(
                max_mean_normalized_attention_score = max(mean_normalized_attention_score),
                .groups = "drop") %>%
            dplyr::filter(
                max_mean_normalized_attention_score > .1,
                !(as.character(cluster_label) %in% c("cluster 236", "cluster 237"))) %>%
            dplyr::mutate(cluster_index = dplyr::row_number()),
        by = "cluster_label")


cluster_order <- data.frame(
    cluster_index = plot_data %>%
        dplyr::select(cluster_index, time_point, mean_normalized_attention_score) %>%
        tidyr::pivot_wider(
            id_cols = "cluster_index",
            names_from = "time_point",
            values_from = "mean_normalized_attention_score") %>%
        dplyr::select(-cluster_index) %>%
        as.matrix() %>%
        dist() %>%
        seriation::seriate(method = "OLO") %>%
        seriation::get_order()) %>%
    dplyr::mutate(
        cluster_ordered_index = dplyr::row_number())

plot_data <- plot_data %>%
    dplyr::left_join(
        cluster_order,
        by = "cluster_index") %>%
    dplyr::arrange(cluster_ordered_index) %>%
    dplyr::mutate(
        cluster_label = cluster_label %>%
            as.character() %>%
            as.factor() %>%
            forcats::fct_inorder())


plot <- ggplot2::ggplot(
    data = plot_data) +
    ggplot2::theme_bw() +
    ggplot2::geom_tile(
        mapping = ggplot2::aes(
            x = time_point,
            y = cluster_label,
            fill = mean_normalized_attention_score)) +
    ggplot2::scale_x_discrete(
        "Time Point",
        expand = c(0, 0)) +
    ggplot2::scale_y_discrete(
        "Cluster",
        expand = c(0, 0)) +
    ggplot2::theme(legend.position = "bottom") +
    viridis::scale_fill_viridis(
        "Mean\nAtt Score",
        option = "C") +
    ggplot2::theme(
        axis.text.x = ggplot2::element_text(
            angle = 90,
            vjust = 0.5,
            hjust = 1)) +
    ggplot2::theme(
        legend.key.width = ggplot2::unit(.6, "cm"))
        

ggplot2::ggsave(
    filename = "product/figures/TS2/attention_scores/mean_attention_score_by_cluster_time_point_heatmap_signif_20220512.pdf",
    plot = plot,
    width = 2,
    height = 10,
    useDingbats = FALSE)


######


cluster_density_vs_mean_attention <- cell_att_cluster %>%
    dplyr::group_by(time_point, plate_id, row, column, cluster_label) %>%
    dplyr::summarize(
        mean_normalized_attention_score = mean(normalize_att_score),
        cluster_count = dplyr::n(),
        .groups = "drop") %>%
    dplyr::group_by(time_point, plate_id, row, column) %>%
        dplyr::mutate(well_count = sum(cluster_count)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
        log_rate = ifelse(
            cluster_count == 0,
            -10,
            log(cluster_count / well_count)))


plot <- ggplot2::ggplot(
    data = cluster_density_vs_mean_attention) +
    ggplot2::theme_bw() +
    ggplot2::geom_point(
        mapping = ggplot2::aes(
            x = log_rate,
            y = log(mean_normalized_attention_score)),
        alpha = .7,
        shape = 16,
        size = 1) +
    ggplot2::geom_smooth(
        mapping = ggplot2::aes(
            x = log_rate,
            y = log(mean_normalized_attention_score)),
        method = "lm",
        formula = "y ~ x") +
    ggplot2::scale_x_continuous(
        "Cells in cluster per thousand cells in well",
        breaks = c(-10, -9, -8, -7, -6, -5, -4),
        labels = c(0, signif(exp(c(-9, -8, -7, -6, -5, -4)) * 1e4, 2))) +
    ggplot2::scale_y_continuous(
        "Mean Att Score in Cluster",
        breaks = c(-8, -6, -4, -2, -0),
        labels = signif(exp(c(-8, -6, -4, -2, -0)), 2)) +
#    viridis::scale_color_viridis(
#        "Time Point",
#        option = "B",
#        begin = 0.2,
#        end = 0.8,
#        discrete = TRUE,
#        guide = ggplot2::guide_legend(
#            nrow = 1,
#            title.hjust = 0.5)) +
    ggplot2::facet_wrap(facets = dplyr::vars(time_point)) +
    ggplot2::theme(
        legend.position = "bottom")

ggplot2::ggsave(
    filename = "product/figures/TS2/attention_scores/attention_score_by_cluster_density_20220512.pdf",
    plot = plot,
    width = 8,
    height = 6,
    useDingbats = FALSE)

    


