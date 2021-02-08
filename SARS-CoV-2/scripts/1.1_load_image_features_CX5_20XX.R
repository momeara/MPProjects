library(plyr)
library(tidyverse)
library(readxl)
library(RMySQL)
library(magrittr)
library(tictoc)
library(arrow)

load("intermediate_data/plate_map_20XX.Rdata")
load("intermediate_data/plate_map_2021A.Rdata")

source("scripts/database.R")
con <- get_primary_database_connection()
schema <- "covid19primary"

load_image_features_CX5 <- function(con, schema, plate_map){
    con %>% set_schema(schema)
    table_name <- paste0(plate_map$Plate_Name[1], "_Per_Image")
    cat("Loading image data from table '", table_name, "' from schema '", schema, "'\n", sep="")
    features_tbl <- con %>%
        dplyr::tbl(table_name) %>%
        dplyr::select(
           Image_Metadata_Field,
           Image_Metadata_Frame,
           Image_Metadata_PlateID,
           Image_Metadata_WellID,
           ImageNumber,
           Image_Classify_Negative_NumObjectsPerBin,
           Image_Classify_Negative_PctObjectsPerBin,
           Image_Classify_Positive_NumObjectsPerBin,
           Image_Classify_Positive_PctObjectsPerBin,
           Image_Count_Cells,
           Image_Count_Cytoplasm,
           Image_Count_Droplets,
           Image_Count_Nuclei,
           Image_Count_Nucleoli,
           Image_Count_syn_nucs,
           Image_Count_syncytia) %>%
        dplyr::collect(n = Inf)
    features_tbl %>%
        dplyr::left_join(
            plate_map,
            by = c("Image_Metadata_WellID" = "Well_ID"))
}

image_scores_CX5_20XX <- plate_map_20XX %>%
    dplyr::filter(schema=="covid19primary") %>%
    plyr::ddply(c("plate_id"), function(plate_map){
        load_image_features_CX5(con, plate_map$schema[1], plate_map)
    })

image_scores_CX5_20XX %>% save(file="intermediate_data/image_scores_CX5_20XX.Rdata")
image_scores_CX5_20XX %>% arrow::write_parquet("product/image_scores_CX5_20XX.parquet")


##########
image_scores_CX5_2021A <- load_image_features_CX5(con, 'covid19primary', plate_map_2021A)

image_scores_CX5_2021A %>% save(file="intermediate_data/image_scores_CX5_2021A.Rdata")
image_scores_CX5_2021A %>% arrow::write_parquet("product/image_scores_CX5_2021A.parquet")
