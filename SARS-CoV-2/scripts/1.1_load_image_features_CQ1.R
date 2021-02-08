library(plyr)
library(tidyverse)
library(readxl)
library(magrittr)
library(arrow)

source("scripts/load_features_from_database.R")

source("scripts/database.R")
con <- get_primary_database_connection("covid19cq1")



load("intermediate_data/plate_map_1999B.Rdata")
load("intermediate_data/plate_map_2020A.Rdata")
load("intermediate_data/plate_map_2021A.Rdata")

load("intermediate_data/plate_map_20XX.Rdata")

load("intermediate_data/plate_map_TS.Rdata")

load("intermediate_data/plate_map_TS_202008.Rdata")





schema <- "covid19cq1"
objects <- c("Cells", "Cytoplasm", "Droplets", "Nuclei", "Nucleoli", "syn_nucs", "syncytia")

image_scores_CQ1_999A <- load_image_features_CQ1(
    con, schema = schema, plate_map = plate_map_999A, objects = objects)
save(image_scores_CQ1_999A, file = "intermediate_data/image_scores_CQ1_999A.Rdata")
image_scores_CQ1_999A %>% arrow::write_parquet("product/image_scores_CQ1_999A.parquet")

image_scores_CQ1_1999B <- load_image_features_CQ1(
    con, schema = schema, plate_map = plate_map_1999B, objects = objects)
save(image_scores_CQ1_1999B, file = "intermediate_data/image_scores_CQ1_1999B.Rdata")
image_scores_CQ1_1999B %>% arrow::write_parquet("product/image_scores_CQ1_1999B.parquet")

image_scores_CQ1_2020A <- load_image_features_CQ1(
    con, schema = schema, plate_map = plate_map_2020A, objects = objects)
save(image_scores_CQ1_2020A, file = "intermediate_data/image_scores_CQ1_2020A.Rdata")
image_scores_CQ1_2020A %>% arrow::write_parquet("product/image_scores_CQ1_2020A.parquet")

image_scores_CQ1_2021A <- load_image_features_CQ1(
    con, schema = schema, plate_map = plate_map_2021A, objects = objects)
save(image_scores_CQ1_2021A, file = "intermediate_data/image_scores_CQ1_2021A.Rdata")
image_scores_CQ1_2021A %>% arrow::write_parquet("product/image_scores_CQ1_2021A.parquet")


# check that the WellID is mapped correctly
image_scores_CQ1_1999B %>%
    dplyr::distinct(Image_Metadata_WellID, .keep_all = T) %>%
    dplyr::arrange(Image_Metadata_WellID) %>%
    dplyr::mutate(Well_ID = paste0(LETTERS[row], column)) %>%
    dplyr::select(Image_Metadata_WellID, row, column, Well_ID, Condition, rem_label, lf_label) %>%
    data.frame

#################
## 20XX Series ##
#################
schema <- "covid19cq1"
objects <- c("Cells", "Cytoplasm", "Droplets", "Nuclei", "Nucleoli", "syn_nucs", "syncytia")


image_scores_CQ1_20XX <- plate_map_20XX %>%
    plyr::ddply(c("plate_id"), function(plate_map) {
        load_image_features_CQ1(
            con,
            schema = schema,
            plate_map = plate_map,
            objects = objects)
    })

image_scores_CQ1_20XX %>% save(file = "intermediate_data/image_scores_CQ1_20XX.Rdata")
image_scores_CQ1_20XX %>% arrow::write_parquet("product/image_scores_CQ1_20XX.parquet")


#######################
## Time Series 202006 #
#######################


# someone messed up along the way and the well ids are NA for TS2PL3
# fortunately can scrape them out of one of the image filenames
image_scores_CQ1_TS_202008 <- plate_map_TS_202008 %>%
    dplyr::filter(plate_id %in% c("TS2PL1", "TS2PL2")) %>%
    plyr::ddply(c("plate_id"), function(plate_map) {
        load_image_features_CQ1(
            con,
            schema = "covid19cq1",
            plate_map = plate_map,
            objects = c("Cells", "Nuclei", "InfectedCells"))
    })

load_image_features_TS2PL3 <- function(con, schema, plate_map, objects) {
    con %>% set_schema(schema)
    plate_id <- plate_map$plate_id[1]
    table_name <- paste0("SARS_", plate_id, "_Per_Image")
    available_tables <- con %>% DBI::dbListTables()
    if (!(table_name %in% available_tables)) {
        cat("ERROR: Unable to get image features beacuse table '", table_name, "' is not in schema '", schema, "'\n", sep = "")
        return(data.frame())
    }
    cat("Loading image data from table '", table_name, "'\n", sep = "")
    image_features <- con %>%
        dplyr::tbl(table_name) %>%
        dplyr::select(
           Image_Metadata_PlateID,
           Image_Metadata_FieldID,
           Image_Metadata_Frame,
           Image_Metadata_WellID,
           ImageNumber,
           Image_Classify_Negative_NumObjectsPerBin,
           Image_Classify_Negative_PctObjectsPerBin,
           Image_Classify_Positive_NumObjectsPerBin,
           Image_Classify_Positive_PctObjectsPerBin,
           Image_URL_ConA,
           tidyselect::one_of(paste0("Image_Count_", objects))) %>%
        dplyr::collect(n = Inf) %>%
        dplyr::mutate(
            Image_Metadata_FieldID = Image_URL_ConA %>%
                stringr::str_match("F([0-9][0-9][0-9][0-9])T") %>%
                magrittr::extract(, 2) %>%
                unlist(),
            Image_Metadata_WellID = Image_URL_ConA %>%
                stringr::str_match("W([0-9][0-9][0-9][0-9])F") %>%
                magrittr::extract(, 2) %>%
                unlist()) %>%
        dplyr::mutate(
            row = floor((as.numeric(Image_Metadata_WellID) - 1) / 24) + 1,
            column = mod((as.numeric(Image_Metadata_WellID) - 1), 24) + 1) %>%
        dplyr::select(-Image_URL_ConA) %>%
        dplyr::left_join(
            plate_map %>%
                dplyr::select(-Well_ID),
            by = c("row", "column"))
}

image_scores_CQ1_TS2PL3_202008 <- plate_map_TS_202008 %>%
    dplyr::filter(plate_id %in% c("TS2PL3")) %>%
    plyr::ddply(c("plate_id"), function(plate_map) {
        load_image_features_TS2PL3(
            con,
            schema = "covid19cq1",
            plate_map = plate_map,
            objects = c("Cells", "Nuclei", "InfectedCells"))
    })

image_scores_CQ1_TS_202008 <- dplyr::bind_rows(
    image_scores_CQ1_TS_202008,
    image_scores_CQ1_TS2PL3_202008)
    
image_scores_CQ1_TS_202008 %>% save(file = "intermediate_data/image_scores_CQ1_TS_202008.Rdata")
image_scores_CQ1_TS_202008 %>% arrow::write_parquet("product/image_scores_CQ1_TS_202008.parquet")
