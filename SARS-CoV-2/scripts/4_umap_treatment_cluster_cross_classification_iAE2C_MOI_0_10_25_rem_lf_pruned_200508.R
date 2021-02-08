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
    "raw_data/iAEC2_cell_features_200508.tsv")

cell_features <- arrow::read_parquet(
    "product/features_iAE2C_MOI_0_10_25_rem_lf_pruned_200508.parquet",
    col_select = c("Treatment", "MOI"))

#cluster_labels <- arrow::read_parquet(
#    "~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/top_hits_plate_scaled_200519_umap2_2M_15_0.0/hdbscan_clustering_min100.parquet")
#
#treatment_by_cluster <- dplyr::bind_cols(
#    cell_features,
#    cluster_labels) %>%
#    dplyr::filter(plate_id != "0999A") %>%
#    dplyr::count(plate_id, Compound, dose_nM, cluster_label) %>%
#    dplyr::mutate(
#        dose_nM_label = dose_nM %>% signif(3) %>% as.character,
#        cluster_label = factor(cluster_label),
#        Compound = Compound %>% stringr::str_extract("^[a-zA-Z0-9-]+"))


cluster_membership <- arrow:::read_parquet(
    "~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/iAE2C_MOI_0_10_25_rem_lf_200508/roi_membership.parquet")


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

treatment_by_cluster <- dplyr::bind_cols(
    cell_features,
    cluster_labels) %>%
    dplyr::count(MOI, Treatment, cluster_label) %>%
    dplyr::group_by(MOI, Treatment) %>%
    dplyr::mutate(fraction = n/sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
        cluster_label = factor(cluster_label))

treatment_by_cluster %>%
    readr::write_tsv(
        "product/figures/umap_features/iAE2C_MOI_0_10_25_rem_lf_pruned_200508/treatment_by_roi.tsv")



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
        data = treatment_by_cluster,
        mapping = ggplot2::aes(
            x = cluster_label,
            y = Treatment,
            fill = log(n))) +
    ggplot2::facet_wrap(~MOI) +
    ggplot2::scale_x_discrete("Cluster Label") +
    ggplot2::scale_y_discrete("Treatment") +
    viridis::scale_fill_viridis(
        "Cell count")


ggplot2::ggsave(
    filename = "product/figures/umap_features/iAE2C_MOI_0_10_25_rem_lf_pruned_200508/treatment_by_roi.pdf",
    plot = p)

################################
make_treatment_by_compound_roi_plot <- function (
    data,
    tag,
    width = 8,
    height = 8,
    additional_layers=NULL) {
    cat("Making roi plot for tag ", tag, "\n", sep = "")
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
            data = data,
            mapping = ggplot2::aes(
                x = MOI,
                y = Treatment,
                fill = n)) +
        ggplot2::facet_wrap(~cluster_label, nrow = 6) +
        ggplot2::scale_x_discrete("MOI") +
        ggplot2::scale_y_discrete("Treatment") +
        viridis::scale_fill_viridis(
            "Cell count") +
        additional_layers
    ggplot2::ggsave(
        filename = glue::glue(
            "product/figures/umap_features/iAE2C_MOI_0_10_25_rem_lf_pruned_200508/",
            "treatment_by_compound_roi_{tag}.pdf"),
        plot = p,
        width = width,
        height = height)
}


make_treatment_by_compound_roi_plot(
    data = treatment_by_cluster,
    tag = "all",
    width = 6,
    height = 10,
    additional_layers = list(
        viridis::scale_fill_viridis(
            "Cell count")))


    



################################

cell_features <- arrow::read_parquet(
    "product/features_iAE2C_MOI_0_10_25_rem_lf_pruned_200508.parquet",
    col_select = c(
        "Metadata_WellID_Cells",
        "Metadata_FieldID_Cells",
        "Metadata_FileLocation_Cells",
        "Treatment",
        "MOI",
        "TreatmentxMOI",
        "ImageNumber_Cells",
        "ObjectNumber_Cells",
        "AreaShape_Center_X_Cells",
        "AreaShape_Center_Y_Cells"))


cluster_cell_examples <- dplyr::bind_cols(
    cell_features,
    cluster_labels) %>%
    dplyr::group_by(cluster_label) %>%
    dplyr::sample_n(50)


cluster_cell_examples %>%
    readr::write_tsv("product/figures/umap_features/iAE2C_MOI_0_10_25_rem_lf_pruned_200508/roi_example_cells.tsv")


################################################

cell_features <- arrow::read_parquet(
    "product/top_hit_cells_plate_scaled_200522a.parquet",
    col_select = c(
        "Compound",
        "plate_id",
        "dose_nM",
        "infectivity_score",
        "infectivity_score_bin",
        feature_columns$feature))

data <- dplyr::bind_cols(
    cell_features,
    cluster_labels) %>%
    dplyr::filter(as.character(cluster_label) != "ROI Inf")

feature_importance <- data %>%
    dplyr::filter(
        cluster_label != "ROI 14") %>%
    dplyr::distinct(cluster_label) %>%
    plyr::adply(1, function(df) {
        roi_label <- as.character(df$cluster_label[1])
        cat("Computing feature importance for ROI '", roi_label, "'\n", sep = "")
        data_matrix <- data  %>%
            dplyr::mutate(label = dplyr::case_when(
                as.character(cluster_label) == "ROI 14" ~ "normal",
                as.character(cluster_label) == roi_label ~ "ROI",
                TRUE ~ NA_character_) %>%
                    factor(
                        levels = c("ROI", "normal"),
                        labels = c("ROI", "normal"))) %>%
            dplyr::select(
                label,
                tidyselect::one_of(feature_columns$feature)) %>%
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
                feature_columns %>%
                dplyr::select(feature)) %>%
            dplyr::rename(feature_score = Overall) %>%
            dplyr::arrange(desc(feature_score)) %>%
            head(20) %>%
            dplyr::mutate(cluster_label = roi_label)
    })


#                                              Overall
# Cells_RadialDistribution_RadialCV_Virus_8of8  100.00
# Cells_RadialDistribution_MeanFrac_Virus_3of8   95.68
# Cells_RadialDistribution_MeanFrac_Virus_2of8   94.99
# Cytoplasm_Intensity_MassDisplacement_Virus     94.52
# Cytoplasm_Intensity_StdIntensity_Virus         93.76
# Cells_RadialDistribution_MeanFrac_Virus_4of8   93.44
# Cells_Intensity_StdIntensity_Virus             93.34
# Cells_Intensity_MaxIntensity_Virus             93.32
# Cytoplasm_Intensity_MaxIntensity_Virus         93.08
# Cells_Intensity_MassDisplacement_Virus         92.60
# Cytoplasm_Intensity_StdIntensityEdge_Virus     92.21
# Cells_RadialDistribution_MeanFrac_Virus_1of8   91.67
# Cytoplasm_Intensity_MaxIntensityEdge_Virus     91.64
# Cells_Intensity_StdIntensityEdge_Virus         91.48
# Cells_Intensity_MaxIntensityEdge_Virus         90.89
# Cells_RadialDistribution_RadialCV_Virus_7of8   88.99
# Cells_RadialDistribution_MeanFrac_Virus_5of8   86.99
# Cytoplasm_Intensity_MeanIntensityEdge_Virus    84.05
# Cells_Intensity_MeanIntensityEdge_Virus        83.53
# Cells_RadialDistribution_MeanFrac_Virus_8of8   81.74


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
    filename = glue::glue(
        "product/figures/umap_features/top_hits_plate_scaled_200522a/",
        "roi_importance_200325.pdf"),
    plot = p,
    width = 15,
    height = 20)

feature_importance %>%
    readr::write_tsv("product/figures/umap_features/top_hits_plate_scaled_200522a/roi_importance_200325_source_data.tsv")

