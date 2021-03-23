
library(plyr)
library(tidyverse)
library(magrittr)
library(arrow)
library(MPStats)

source("parameters.R")

##########################
# Time Series 1 (202006) #
##########################

image_scores_CQ1_TS <- arrow::read_parquet(
    file = "intermediate_data/image_scores_CQ1_TS.parquet")
image_metadata_columns_TS <- tibble::tibble(
    feature = image_scores_CQ1_TS %>% names())

cell_features_TS6h <- arrow::read_parquet(
    file = "intermediate_data/covid19cq1_SARS_TS6h_Cell_MasterDataTable.parquet")

objects <- c("Cytoplasm", "Cells", "Nuclei", "InfectedCells")
dyes <- c("NP", "Spike", "ConA", "Hoe")
coordinates <- c("X", "Y", "Z")

# these are the features we want to exclude
cell_feature_columns_TS <- tibble::tibble(
    feature = cell_features_TS6h %>% names(),
    transform = "identity") %>%
    dplyr::anti_join(image_metadata_columns_TS, by = "feature") %>%
    dplyr::filter(!(feature %in% c("schema", "dose_nM"))) %>%
    dplyr::filter(!(feature %in% c("Nuclei_Distance_Centroid_InfectedCells"))) %>%    
    dplyr::anti_join(
        expand.grid(child = objects, parent = objects) %>%
        dplyr::mutate(feature = paste0(child, "_Parent_", parent)),
        by = "feature") %>%
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

# antijoin against the ones we want to keep to make sure we haven't missed any
cell_feature_columns_TS %>%
    dplyr::anti_join(
        expand.grid(
            object = objects,
            feature = c(
                "Intensity"),
            statistic = c(
                "IntegratedIntensityEdge",
                "IntegratedIntensity",
                "LowerQuartileIntensity",
                "MADIntensity",
                "MassDisplacement",
                "MaxIntensityEdge",
                "MaxIntensity",
                "MeanIntensityEdge",
                "MeanIntensity",
                "MedianIntensityEdge",
                "MedianIntensity",
                "MinIntensityEdge",
                "MinIntensity",
                "StdIntensityEdge",
                "StdIntensity",
                "UpperQuartileIntensity"),
            dye = dyes) %>%
        dplyr::mutate(
            feature = paste(object, feature, statistic, dye, sep = "_")),
        by = "feature") %>%
    dplyr::anti_join(
        expand.grid(
            object = objects,
            feature = c(
                "AreaShape"),
            statistic = c(
                "Area",
                "Compactness",
                "Eccentricity",
                "EulerNumber",
                "Extent",
                "FormFactor",
                "MajorAxisLength",
                "MaxFeretDiameter",
                "MaximumRadius",
                "MeanRadius",
                "MedianRadius",
                "MinFeretDiameter",
                "MinorAxisLength",
                "Orientation",
                "Perimeter",
                "Solidity")) %>%
        dplyr::mutate(
            feature = paste(object, feature, statistic, sep = "_")),
        by = "feature") %>%
    dplyr::anti_join(
        expand.grid(
            object = objects,
            zernike_level_1 = 0:10,
            zernike_level_2 = 0:10) %>%
        dplyr::mutate(
            feature = paste(
                object, "AreaShape_Zernike", zernike_level_1, zernike_level_2, sep = "_")),
        by = "feature") %>%
    dplyr::anti_join(
        expand.grid(
            object = objects,
            statistic = c(
                "FracAtD",
                "MeanFrac",
                "RadialCV"),
            dye = dyes,
            fraction = 1:8) %>%
        dplyr::mutate(
            feature = paste0(
                object, "_RadialDistribution_", statistic, "_", dye, "_", fraction, "of8")),
        by = "feature") %>%
    dplyr::anti_join(
        expand.grid(
            object = objects,
            statistic = c("ZernikeMagnitude"),
            dye = dyes,
            zernike_level_1 = 0:9,
            zernike_level_2 = 0:9) %>%
        dplyr::mutate(
            feature = paste(
                object, "RadialDistribution",
                statistic, dye,
                zernike_level_1, zernike_level_2, sep = "_")),
        by = "feature")


#check that there are no missing values
cell_features_TS6h %>% nrow()
cell_features_TS6h %>%
    dplyr::select(
        tidyselect::one_of(cell_feature_columns_TS$feature)) %>%
    tidyr::drop_na() %>%
    nrow()

cell_features_TS6h %>% summary()

cell_feature_columns_TS %>%
    readr::write_tsv("product/cell_feature_columns_TS.tsv")

cell_metadata_columns_TS <- cell_features_TS6h %>%
    names() %>%
    data.frame(feature = .) %>%
    dplyr::anti_join(
        cell_feature_columns_TS6h,
        by = "feature")

cell_metadata_columns_TS %>%
    readr::write_tsv("product/cell_metadata_columns_TS.tsv")


##########################
# Time Series 2 (202008) #
##########################

image_scores_CQ1_TS_202008 <- arrow::read_parquet(
    "intermediate_data/image_scores_CQ1_TS_202008.parquet")

image_metadata_columns_TS_202008 <- tibble::tibble(
    feature = image_scores_CQ1_TS_202008 %>% names())

cell_features_TS2PL1 <- arrow::read_parquet(
    file = "intermediate_data/covid19cq1_SARS_TS2PL1_Cell_MasterDataTable.parquet")

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






