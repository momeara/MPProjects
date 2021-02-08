
library(plyr)
library(tidyverse)
library(arrow)
library(caret)


source("scripts/make_scatter_boards.R")


data_path <- "intermediate_data/infected_patch_1999B_2020A_2021A_20201017"

viral_features <- arrow::read_parquet(
    file = paste0(data_path, "/viral_features.parquet"))
viral_feature_columns <- readr::read_tsv(
    file = paste0(data_path, "/viral_feature_columns.tsv"))
viral_metadata_columns <- readr::read_tsv(
    file = paste0(data_path, "/viral_metadata_columns.tsv"))

nuclei_features <- arrow::read_parquet(
    file = paste0(data_path, "/nuclei_features.parquet"))
nuclei_feature_columns <- readr::read_tsv(
    file = paste0(data_path, "/nuclei_feature_columns.tsv"))
nuclei_metadata_columns <- readr::read_tsv(
    file = paste0(data_path, "/nuclei_metadata_columns.tsv"))

syn_nuc_features <- arrow::read_parquet(
    file = paste0(data_path, "/syn_nuc_features.parquet"))
syn_nuc_feature_columns <- readr::read_tsv(
    file = paste0(data_path, "/syn_nuc_feature_columns.tsv"))
syn_nuc_metadata_columns <- readr::write_tsv(
    file = paste0(data_path, "/syn_nuc_metadata_columns.tsv"))

puncta_features <- arrow::read_parquet(
    file = paste0(data_path, "/puncata_features.parquet"))
puncta_metadata_columns <- readr::read_tsv(
    file = paste0(data_path, "/puncta_metadata_columns.tsv"))


#########################
# object counts by dose #
#########################

viral_counts <- viral_features %>%
    dplyr::count(
        condition,
        drug_1,
        drug_1_units,
        drug_1_concentration,
        drug_1_label,
        drug_2,
        drug_2_units,
        drug_2_concentration,
        drug_2_label)



make_scatter_boards(
    features = viral_features,
    feature_x = "Children_Nuclei_Count",
    feature_y = "AreaShape_Area",
    feature_color = "AreaShape_Compactness",
    scales = NULL,
    output_dir = "product/infected_patch_1999B_2020A_2021A_20201017/viral_features")

make_scatter_boards(
    features = viral_features,
    feature_x = "Intensity_MaxIntensityEdge_NP",
    feature_y = "Intensity_MedianIntensity_NP",
    feature_color = "AreaShape_Area",
    scales = list(
        scale_x_log10("Max NP Edge Intensity"),
        scale_y_log10("Median NP Intensity")),
    output_dir = "product/infected_patch_1999B_2020A_2021A_20201017/viral_features")


make_scatter_boards(
    features = nuclei_features %>%
        dplyr::group_by(plate_id) %>%
        dplyr::sample_n(500000) %>%
        dplyr::ungroup(),
    feature_x = "Intensity_MedianIntensity_Hoe",
    feature_y = "Intensity_MedianIntensity_NP",
    feature_color = "AreaShape_Area",
    scales = list(
        scale_x_log10("Median Hoe Intensity"),
        scale_y_log10("Median NP Intensity")),
    output_dir = "product/infected_patch_1999B_2020A_2021A_20201017/nuclei_features")


features <- nuclei_features %>%
    dplyr::select(
        ImageNumber,
        plate_id,
        drug_1_concentration,
        drug_2_concentration,
        drug_1_label,
        drug_2_label,
        Parent_ViralObj,
        nuclei_max_edge_NP_intensity = Intensity_MaxIntensityEdge_NP,
        nuclei_median_NP_intensity = Intensity_MedianIntensity_NP) %>%
    dplyr::group_by(plate_id) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(
        viral_features %>%
        dplyr::transmute(
            ImageNumber,
            Parent_ViralObj = Number_Object_Number,
            viral_area = AreaShape_Area,
            viral_parimeter = AreaShape_Perimeter,
            log_num_nuclei = log10(Mean_Nuclei_Number_Object_Number)),
        by = c("ImageNumber", "Parent_ViralObj")) %>%
    dplyr::filter(!is.na(viral_area))

make_scatter_boards(
    features = features,
    feature_x = "nuclei_median_NP_intensity",
    feature_y = "viral_area",
    feature_color = "log_num_nuclei",
    scales = list(
        scale_x_log10("Median Nuclei NP Intensity")),
    output_dir = "product/infected_patch_1999B_2020A_2021A_20201017/nuclei_viral_features")


make_scatter_boards(
    features = features,
    feature_x = "nuclei_median_NP_intensity",
    feature_y = "log_num_nuclei",
    feature_color = "viral_area",
    scales = list(
        scale_x_log10("Median Nuclei NP Intensity")),
    output_dir = "product/infected_patch_1999B_2020A_2021A_20201017/nuclei_viral_features")

make_scatter_boards(
    features = features,
    feature_x = "nuclei_max_edge_NP_intensity",
    feature_y = "viral_parimeter",
    feature_color = "log_num_nuclei",
    scales = list(
        scale_x_log10("Nuclei Max Edge NP Intensity")),
    output_dir = "product/infected_patch_1999B_2020A_2021A_20201017/nuclei_viral_features")


















# top_discriminating_features_pls <- function(
object_features = viral_features
object_feature_columns = viral_feature_columns
stratify_features = c("plate_id")
#)


object_features <- object_features %>%
    dplyr::select(
        tidyselect::feature
    dplyr::filter(condition %in% c("PC", "NC")) %>%
    dplyr::mutate(label = factor(condition))

object_feature_coluns <- 

in_training <- caret::createDataPartition(
    object_features$label,
    p = .80,
    list = FALSE)
data_train <- object_Features %>% dplyr::slice(in_training)
data_test  <- object_Featurse %>% dplyr::slice(-in_training)
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


