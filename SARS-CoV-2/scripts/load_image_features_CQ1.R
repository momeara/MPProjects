

library(plyr)
library(tidyverse)
library(magrittr)
library(argparser)

source("scripts/database.R")

load_image_features_CQ1 <- function(con, schema, plate_map, objects) {
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
           tidyselect::one_of(paste0("Image_Count_", objects))) %>%
        dplyr::collect(n = Inf) %>%
        dplyr::mutate(
            row = floor((as.numeric(Image_Metadata_WellID) - 1) / 24) + 1,
            column = mod((as.numeric(Image_Metadata_WellID) - 1), 24) + 1) %>%
        dplyr::left_join(
            plate_map %>%
                dplyr::select(-Well_ID),
            by = c("row", "column"))
}
