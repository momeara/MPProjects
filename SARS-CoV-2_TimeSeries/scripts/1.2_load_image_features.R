library(plyr)
library(tidyverse)
library(readxl)
library(magrittr)
library(arrow)
library(MPStats)

source("scripts/load_features_from_database.R")

con <- MPStats::get_database_connection(
    database_parameters = parameters$databases$covid19cq1)

load("intermediate_data/plate_map_TS.Rdata")
load("intermediate_data/plate_map_TS_202008.Rdata")



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
            table_name = paste0("SARS_", plate_map$plate_id[1], "_Per_Image"),
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
