
library(plyr)
library(tidyverse)
library(RMySQL)
library(magrittr)
library(tictoc)
library(arrow)

source("scripts/database.R")
con <- get_primary_database_connection()

source("scripts/load_features_from_database.R")


plate_ids <- readr::read_tsv("raw_data/plate_ids.tsv")
feature_columns <- readr::read_tsv("raw_data/cell_feature_columns.tsv")

load("intermediate_data/image_scores_CX5_100X.Rdata")
load("intermediate_data/image_scores_CX5_20XX.Rdata")
load("intermediate_data/image_scores_CQ1_999A.Rdata")

load("intermediate_data/image_scores_CQ1_1999B.Rdata")
load("intermediate_data/image_scores_CQ1_2020A.Rdata")

load("intermediate_data/image_scores_CQ1_2021A.Rdata")
load("intermediate_data/image_scores_CQ1_20XX.Rdata")
load("intermediate_data/image_scores_CQ1_TS.Rdata")

load("intermediate_data/image_scores_CQ1_TS_202008.Rdata")


## 100X series: 5-dose qHTS for 5 plates of drugs
plate_ids %>%
    dplyr::filter(master_plate_id == "1005", plate_id %>% stringr::str_detect("C$")) %>%
    plyr::a_ply(1, function(df) {
        collect_cell_features(
            con = con,
            schema = df$schema[1],
            image_scores = image_scores_CX5_100X,
            plate_id = df$plate_id[1])
    })

# 20XX series: 10-point dose-response in triplicate for top 140 compounds
plate_ids %>%
    dplyr::filter(schema == "covid19primary") %>%
    dplyr::filter(plate_id %>% stringr::str_detect("^2")) %>%
    plyr::a_ply(1, function(df) {
        collect_cell_features(
            con = con,
            schema = df$schema[1],
            image_scores = image_scores_CX5_20XX,
            plate_id = df$plate_id[1])
    })

# 20XX serise on the CQ1
plate_ids %>%
    dplyr::filter(
        schema == "covid19cq1",
        master_plate_id %>% stringr::str_detect("^20..")) %>%
    plyr::a_ply(1, function(df) {
        collect_cell_features(
            con = con,
            schema = df$schema[1],
            image_scores = image_scores_CQ1_20XX,
            plate_id = df$plate_id[1])
    })

readr::read_tsv("raw_data/plate_ids.tsv") %>%
    dplyr::filter(schema == "covid19cq1") %>%
    plyr::a_ply(1, function(df) {
        collect_cell_features(
            con = con,
            schema = df$schema[1],
            image_scores = image_scores_CQ1_20XX,
            plate_id = df$plate_id[1])
    })


collect_cell_features(
    con = con,
    schema = 'covid19cq1',
    image_scores = image_scores_CQ1_999A,
    plate_id = "0999A")

collect_cell_features(
    con = con,
    schema = "covid19cq1",
    image_scores = image_scores_CQ1_1999B,
    plate_id = "1999B")

collect_cell_features(
    con = con,
    schema = "covid19cq1",
    image_scores = image_scores_CQ1_2020A,
    plate_id = "2020A")

collect_cell_features(
    con=con,
    schema = "covid19cq1",
    image_scores = image_scores_CQ1_2021A,
    plate_id = "2021A")


collect_cell_features(
    con = con,
    schema = "covid19primary",
    image_scores = image_scores_CX5_2021A,
    plate_id = "2021A")





readr::read_tsv("raw_data/plate_ids.tsv") %>%
    dplyr::filter(schema=="covid19cq1") %>%
    dplyr::filter(plate_id %in% c("2006A", "2007A", "2008A", "2009A")) %>%
    plyr::a_ply(1, function(df){
        collect_nucleoli_features(
            con=con,
            schema=df$schema[1],
            plate_id=df$plate_id[1],
            image_scores=image_scores_CQ1_20XX)
    })

collect_nucleoli_features(
    con=con,
    schema='covid19cq1',
    plate_id="1999B",
    image_scores=image_scores_1999B)

collect_nucleoli_features(
    con=con,
    schema='covid19cq1',
    plate_id="2020A",
    image_scores=image_scores_2020A)


readr::read_tsv("raw_data/plate_ids.tsv") %>%
    dplyr::filter(schema == "covid19primary") %>%    
    dplyr::filter(plate_id %>% stringr::str_detect("^10")) %>%
    plyr::a_ply(1, function(df){
        collect_syncytia_features(
            con=con,
            schema=df$schema[1],
            image_scores=image_scores_CX5_100X)
    })

readr::read_tsv("raw_data/plate_ids.tsv") %>%
    dplyr::filter(schema == "covid19primary") %>%
    dplyr::filter(plate_id %>% stringr::str_detect("^20")) %>%
    plyr::a_ply(1, function(df){
        collect_syncytia_features(
            con=con,
            schema=df$schema[1],
            image_scores=image_scores_CX5_20XX)
    })


collect_syncytia_features(
    con,
    schema='covid19cq1',
    plate_id='1999B',
    image_scores=image_scores_1999B)

collect_syncytia_features(
    con,
    schema='covid19cq1',
    plate_id='2020A',
    image_scores=image_scores_2020A)


readr::read_tsv("raw_data/plate_ids.tsv") %>%
    dplyr::filter(schema=="covid19cq1") %>%
    dplyr::filter(plate_id %in% c("2006A", "2007A", "2008A", "2009A")) %>%
    plyr::a_ply(1, function(df){
        collect_syncytia_features(
            con=con,
            schema=df$schema[1],
            plate_id=df$plate_id[1],
            image_scores=image_scores_CQ1_20XX)
    })

##########################
# Time Series 1 (202006) #
##########################

readr::read_tsv("raw_data/plate_ids.tsv") %>%
    dplyr::filter(schema == "covid19cq1") %>%
    dplyr::filter(plate_id %in% c("TS18h")) %>%
    plyr::a_ply(1, function(df) {
        collect_cell_features(
            con = con,
            schema = df$schema[1],
            plate_id = df$plate_id[1],
            image_scores = image_scores_CQ1_TS)
        collect_InfectedCells_features(
            con = con,
            schema = df$schema[1],
            plate_id = df$plate_id[1],
            image_scores = image_scores_CQ1_TS)
    })



readr::read_tsv("raw_data/plate_ids.tsv") %>%
    dplyr::filter(schema == "covid19cq1") %>%
    dplyr::filter(plate_id %in% c("TS3h", "TS6h", "TS12h", "TS18h", "TS24h", "TS48h")) %>%
    plyr::a_ply(1, function(df) {
        collect_cell_features(
            con = con,
            schema = df$schema[1],
            plate_id = df$plate_id[1],
            image_scores = image_scores_CQ1_TS)
        collect_InfectedCells_features(
            con = con,
            schema = df$schema[1],
            plate_id = df$plate_id[1],
            image_scores = image_scores_CQ1_TS)
    })

image_scores_CQ1_TS <- arrow::read_parquet(
    file = "product/image_scores_CQ1_TS.parquet")
image_metadata_columns_TS <- tibble::tibble(
    feature = image_scores_CQ1_TS %>% names())

cell_features_TS6 <- arrow::read_parquet(
    file = "product/covid19cq1_SARS_TS6h_Cell_MasterDataTable.parquet")


objects <- c("Cytoplasm", "Cells", "Nuclei", "InfectedCells")
dyes <- c("NP", "Spike", "ConA", "Hoe")
coordinates <- c("X", "Y", "Z")

# these are the features we want to exclude
cell_feature_columns_TS6h <- tibble::tibble(
    feature = cell_features %>% names(),
    transform = "identity") %>%
    dplyr::anti_join(image_metadata_columns_TS, by="feature") %>%
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

# these are the features we want to keep
cell_feature_columns_TS6h %>%
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
cell_features %>% nrow()
cell_features %>%
    dplyr::select(
        tidyselect::one_of(cell_feature_columns_TS6h$feature)) %>%
    tidyr::drop_na() %>%
    nrow()

cell_feature_columns_TS6h %>%
    readr::write_tsv("raw_data/cell_feature_columns_TS.tsv")
cell_metadata_columns_TS %>%
    readr::write_tsv("raw_data/cell_metadata_columns_TS.tsv")


cell_features %>% summary()


##########################
# Time Series 2 (202008) #
##########################

image_scores_CQ1_TS_202008 <- arrow::read_parquet(
    file = "product/image_scores_CQ1_TS_202008.parquet")

readr::read_tsv("raw_data/plate_ids.tsv") %>%
    dplyr::filter(schema == "covid19cq1") %>%
    dplyr::filter(plate_id %in% c("TS2PL1", "TS2PL2", "TS2PL3")) %>%
    plyr::a_ply(1, function(df) {
        collect_cell_features(
            con = con,
            schema = df$schema[1],
            plate_id = df$plate_id[1],
            image_scores = image_scores_CQ1_TS_202008)
        collect_InfectedCells_features(
            con = con,
            schema = df$schema[1],
            plate_id = df$plate_id[1],
            image_scores = image_scores_CQ1_TS_202008)
    })

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



cell_feature_columns_TS_202008 %>%
    readr::write_tsv("raw_data/cell_feature_columns_TS_202008.tsv")
image_metadata_columns_TS_202008 %>%
    readr::write_tsv("raw_data/cell_metadata_columns_TS_202008.tsv")


cell_features %>% summary()
