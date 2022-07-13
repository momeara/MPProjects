
library(plyr)
library(tidyverse)
library(magrittr)
library(arrow)
library(MPStats)

source("parameters.R")

##########################
# Time Series 2 (202008) #
##########################

image_scores_CQ1_TS_202008 <- arrow::read_parquet(
    "product/image_scores_CQ1_TS_202008.parquet")

image_metadata_columns_TS_202008 <- tibble::tibble(
    feature = image_scores_CQ1_TS_202008 %>% names())

cell_features_TS2PL1 <- arrow::read_parquet(
    file = "product/covid19cq1_SARS_TS2PL1_Cell_MasterDataTable.parquet")

objects <- c("Cells", "Nuclei", "Cytoplasm", "InfectedCells")
dyes <- c("ConA", "Hoe", "NP", "Spike")
coordinates <- c("X", "Y", "Z")

# these are the features we want to exclude
cell_feature_columns_TS_202008 <- tibble::tibble(
    feature = cell_features_TS2PL1 %>% names(),
    transform = "identity") %>%
    dplyr::anti_join(image_metadata_columns_TS_202008, by = "feature") %>%
    dplyr::anti_join(
        expand.grid(child = objects, parent = objects) %>%
        dplyr::mutate(feature = paste0(child, "_Parent_", parent)),
        by = "feature") %>%
    dplyr::filter(!(feature %in% c("Nuclei_Distance_Centroid_InfectedCells"))) %>%
    dplyr::anti_join(
        expand.grid(
            object = objects,
            feature = c(
                "Location",
                "AreaShape"),
            coordinate = coordinates) %>%
        dplyr::mutate(
            feature = paste(object, feature, "Center", coordinate, sep = "_")),
        by = "feature") %>%
    dplyr::anti_join(
        expand.grid(
            object = objects,
            statistic = c(
                "MaxIntensity",
                "CenterMassIntensity"),
            coordinate = coordinates,
            dye = dyes) %>%
        dplyr::mutate(
            feature = paste(object, "Location", statistic, coordinate, dye, sep = "_")),
        by = "feature") %>%
    dplyr::anti_join(
        expand.grid(
            object = objects,
            class = c("Positive", "Negative")) %>%
        dplyr::mutate(
            feature = paste(object, "Classify", class, sep = "_")),
        by = "feature") %>%
    dplyr::anti_join(
        expand.grid(
            parent = objects,
            child = objects) %>%
        dplyr::mutate(
            feature = paste(parent, "Children", child, "Count", sep = "_")),
        by = "feature") %>%
    dplyr::anti_join(
        data.frame(objects) %>%
        dplyr::mutate(feature = paste0(objects, "_Number_Object_Number")),
        by = "feature")

cell_feature_columns_TS_202008 %>% summary()

cell_feature_columns_TS_202008 %>%
    readr::write_tsv("product/cell_feature_columns_TS_202008.tsv")


cell_metadata_columns_TS_202008 <- cell_features_TS2PL1 %>%
    names() %>%
    data.frame(feature = .) %>%
    dplyr::anti_join(
        cell_feature_columns_TS_202008,
        by = "feature")

image_metadata_columns_TS_202008 %>%
    readr::write_tsv("product/cell_metadata_columns_TS_202008.tsv")






