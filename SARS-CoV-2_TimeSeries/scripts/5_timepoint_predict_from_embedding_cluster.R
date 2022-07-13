


library(tidyverse)
library(arrow)
library(brms)
library(fastknn)

cell_feature_columns <- readr::read_tsv(
    "product/cell_feature_columns_TS_202008.tsv",
    show_col_types = FALSE)

cell_metadata_columns <- readr::read_tsv(
    "product/cell_metadata_columns_TS_202008.tsv",
    show_col_types = FALSE)


############################################3

k <- 100
resolution <- 1e-4
num_iter <- 10
cell_clusters <- dplyr::bind_cols(
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
    


cluster_count <- cell_clusters %>%
    dplyr::count(
        time_point, plate_id, well_id, cluster_label,
        name = "cluster_count",
        .drop = FALSE) %>%
    dplyr::group_by(time_point, plate_id, well_id) %>%
    dplyr::mutate(well_count = sum(cluster_count)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(log_rate = ifelse(
        cluster_count == 0,
        -10,
        log(cluster_count / well_count)))


model_8 <- lm(
    formula = log_rate ~ 0 + time_point + cluster_label + time_point:cluster_label,
    data = cluster_count %>% dplyr::filter(
        time_point %in% c("Uninfected", "8 hours")) %>%
        dplyr::slice_sample(n = 2000))

modelz <- brms::brm(
    formula = log_rate ~  time_point * cluster_label,
    data = cluster_count %>% dplyr::filter(
        time_point %in% c("Uninfected", "8 hours")) %>%
        dplyr::slice_sample(n = 2000))

coefs_8 <- data.frame(
    coefficient = model_8 %>% coefficients() %>% names(),
    coefficient_value = model_8 %>% coefficients()) %>%
    dplyr::mutate(
        interaction = coefficient %>%
            stringr::str_detect(":"),
        cluster_label = coefficient %>%
            stringr::str_extract("cluster_label.*$") %>%
            stringr::str_replace("cluster_label", "")) %>%
    dplyr::arrange(
        dplyr::desc(coefficient_value))


plot_data <- cluster_count %>%
    dplyr::semi_join(
        coefs %>%
        dplyr::filter(
            interaction,
            (coefficient_value > 4.8)),
        by = "cluster_label") %>%
    dplyr::bind_rows(
        cluster_count %>%
        dplyr::filter(
            cluster_label %in% paste("cluster ", 1:5)))


plot <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_jitter(
        data = plot_data,
        mapping = ggplot2::aes(
            x = time_point,
            y = log_rate,
            color = cluster_label),
        width = .15,
        size = .8,
        alpha = .3,
        shape = 16) +
    ggplot2::geom_line(
        data = plot_data %>%
            dplyr::group_by(time_point, cluster_label) %>%
            dplyr::summarize(log_rate = mean(log_rate), .groups = "drop"),
        mapping = ggplot2::aes(
            x = time_point,
            y = log_rate,
            color = cluster_label,
            group = cluster_label),
        size = 1.3,
        alpha = .8) +
    ggplot2::scale_x_discrete("Hours Post Infection") +
    ggplot2::scale_y_continuous(
        "Cells per thousand in well",
        breaks = c(-10, -9, -8, -7, -6, -5, -4),
        labels = c(0, signif(exp(c(-9, -8, -7, -6, -5, -4)) * 1e4, 2))) +
    ggplot2::scale_color_discrete(
        "Cluster Label")

        
ggplot2::ggsave(
    filename = "product/figures/TS2/time_point_by_sig_clusters_20220501.pdf",
    plot = plot,
    width = 7,
    height = 6,
    useDingbats = FALSE)


####
model_data <- cluster_count %>%
    dplyr::semi_join(
        coefs %>%
        dplyr::filter(
            interaction,
            (coefficient_value > 4.8)),
        by = "cluster_label") %>%

model <- mgcv::gam(
    formula = time_point ~ s(log_rate),
    data = model_data)


###################################
# Test/train split for clustering #
###################################

k <- 100
resolution <- 1e-4
num_iter <- 10
cell_features <- dplyr::bind_cols(
    arrow::read_parquet(
        file = "product/TS2_plate_scaled_Cell_MasterDataTable.parquet",
        col_select = cell_metadata_columns$feature),
    arrow::read_parquet(
        file = "intermediate_data/UMAP_embedding_into_TS2_2M_epochs=2000_20200901/umap_embedding.parquet")) %>%
    dplyr::mutate(
        time_point = factor(
            x = time_point,
            levels = c(
                "Uninfected",
                "8 hours", "18 hours", "24 hours", "30 hours", "36 hours", "48 hours"),
            labels = c(
                "Uninfected",
                "8 hours", "18 hours", "24 hours", "30 hours", "36 hours", "48 hours")),
        well_id = paste0(toupper(letters[column]), row))
    

cell_features_train <- dplyr::bind_cols(
    cell_features %>%
        dplyr::filter(row <= 8),
    arrow::read_parquet(
        file = paste0(
        "intermediate_data/UMAP_embedding_into_TS2_2M_epochs=2000_20200901/",
        "clusters_train_leiden_k=", k, "_res=", resolution, "_num_iter=", num_iter, ".parquet"))) %>%
    dplyr::arrange(cluster_label) %>%
    dplyr::mutate(
        cluster_code = cluster_label,
        cluster_label = paste("cluster ", cluster_label) %>% as.factor() %>% forcats::fct_inorder())

cell_features_test <- cell_features %>%
    dplyr::filter(row > 8)

knn_labeler <- fastknn::fastknn(
    xtr = cell_features_train %>%
        dplyr::select(UMAP_1, UMAP_2) %>%
        as.matrix(),
    ytr = cell_features_train %>%
        purrr::pluck("cluster_code"),
    xte = cell_features_test %>%
        dplyr::select(UMAP_1, UMAP_2) %>%
        as.matrix(),
    k = 1)

cell_features_test <- cell_features_test %>%
    dplyr::mutate(
        cluster_label_pred = colnames(knn_labeler$prob)[max.col(knn_labeler$prob)])

cluster_count <- cell_features_test %>%
    dplyr::count(
        time_point, plate_id, well_id, cluster_label_pred,
        name = "cluster_count",
        .drop = FALSE) %>%
    dplyr::group_by(time_point, plate_id, well_id) %>%
    dplyr::mutate(well_count = sum(cluster_count)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(log_rate = ifelse(
        cluster_count == 0,
        -10,
        log(cluster_count / well_count)))

model <- lm(
    formula = log_rate ~ 0 + time_point + cluster_label_pred + time_point:cluster_label_pred,
    data = cluster_count %>% dplyr::filter(
        time_point %in% c("Uninfected", "48 hours")))

model <- brms::brm(
    formula =  cluster_count | trials(well_count) ~ (well_id | plate_id) + time_point + cluster_label,
    data = cluster_count,
    family = binomial)

coefs <- data.frame(
    coefficient = model %>% coefficients() %>% names(),
    coefficient_value = model %>% coefficients()) %>%
    dplyr::mutate(
        interaction = coefficient %>%
            stringr::str_detect(":"),
        cluster_label_pred = coefficient %>%
            stringr::str_extract("cluster_label_pred.*$") %>%
            stringr::str_replace("cluster_label_pred", "")) %>%
    dplyr::arrange(
        dplyr::desc(coefficient_value))


plot_data <- cluster_count %>%
    dplyr::left_join(
        coefs %>%
        dplyr::filter(
            interaction),
        by = "cluster_label_pred")


plot <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_line(
        data = plot_data %>%
            dplyr::filter(coefficient_value <= -2) %>%            
            dplyr::group_by(time_point, cluster_label_pred) %>%
            dplyr::summarize(log_rate = mean(log_rate), .groups = "drop"),
        mapping = ggplot2::aes(
            x = time_point,
            y = log_rate,
            group = cluster_label_pred),
        size = .8,
        alpha = .3) +
    ggplot2::geom_line(
        data = plot_data %>%
            dplyr::filter(coefficient_value > -2) %>%
            dplyr::group_by(time_point, cluster_label_pred) %>%
            dplyr::summarize(log_rate = mean(log_rate), .groups = "drop"),
        mapping = ggplot2::aes(
            x = time_point,
            y = log_rate,
            color = cluster_label_pred,
            group = cluster_label_pred),
        size = 1.2,
        alpha = 1) +
    ggplot2::scale_x_discrete("Hours Post Infection") +
    ggplot2::scale_y_continuous(
        "Number of cells in cluster per thousand cells in well",
        breaks = c(-10, -9, -8, -7, -6, -5, -4),
        labels = c(0, signif(exp(c(-9, -8, -7, -6, -5, -4)) * 1e4, 2))) +
    ggplot2::scale_color_discrete(
        "Cluster Label")

        
ggplot2::ggsave(
    filename = "product/figures/TS2/time_point_test_by_clusters_20220503.pdf",
    plot = plot,
    width = 7,
    height = 6,
    useDingbats = FALSE)

ggplot2::ggsave(
    filename = "product/figures/TS2/time_point_test_by_clusters_20220503.png",
    plot = plot,
    width = 7,
    height = 6)


model_data <- plot_data %>%
    dplyr::filter(cluster_label_pred %in% c(20, 32, 99)) %>%
    dplyr::mutate(
        time_point_hours = dplyr::case_when(
            time_point == "Uninfected" ~ 0,
            time_point == "8 hours" ~ 8,
            time_point == "18 hours" ~ 18,
            time_point == "24 hours" ~ 24,
            time_point == "30 hours" ~ 30,
            time_point == "36 hours" ~ 36,
            time_point == "48 hours" ~ 48))

model <- mgcv::gam(
    formula = time_point_hours ~ s(log_rate),
    data = model_data)

model_data <- model_data %>%
    dplyr::mutate(
        time_point_pred = model %>% predict())

plot <- ggplot2::ggplot(model_data) +
    ggplot2::theme_bw() +
    ggplot2::geom_violin(
        mapping = ggplot2::aes(
            x = time_point,
            y = time_point_hours - time_point_pred),
        fill = "black") +
    ggplot2::scale_y_continuous(
        "Predicted Hours Post Infection Residual") +
    ggplot2::scale_x_discrete(
        "Actual Hours Post Infection") +
    ggplot2::coord_flip()



ggplot2::ggsave(
    filename = "product/figures/TS2/predict_time_point_from_cluster_log_rate_smooth_20220503.pdf",
    plot = plot,
    width = 4,
    height = 3,
    useDingbats = FALSE)

ggplot2::ggsave(
    filename = "product/figures/TS2/predict_time_point_from_cluster_log_rate_smooth_20220503.png",
    plot = plot,
    width = 4,
    height = 3)

