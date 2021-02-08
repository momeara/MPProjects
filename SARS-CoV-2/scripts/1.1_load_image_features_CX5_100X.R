library(plyr)
library(tidyverse)
library(readxl)
library(RMySQL)
library(magrittr)
library(tictoc)
library(arrow)

source("scripts/database.R")
con <- get_primary_database_connection()


plate_map <- readxl::read_excel("raw_data/FDA_Platemap_FULL_DrugAnnotation_Version3.xlsx")

load_image_features <- function(con, schema, master_plate_id, plate_id, dose_nM){
    con %>% set_schema(schema) 
    table_name <- paste0("SARS_", plate_id, "_Per_Image")
    cat("Loading image data from table '", table_name, "'\n", sep="")
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
        dplyr::collect(n=Inf)    
    features_tbl <- features_tbl %>%
        dplyr::left_join(
            plate_map %>%
            dplyr::filter(FDA_384_Barcode == master_plate_id) %>%
            dplyr::select(
                Compound,
                `Therapeutic Class`,
                `Smiles String`,
                CDR,
                FDA_384_Well_ID,
                row,
                column),
            by=c("Image_Metadata_WellID" = "FDA_384_Well_ID")) %>%
        dplyr::mutate(
            Compound = dplyr::case_when(
                column %in% c(1, 2) ~ "Negative Control",
                column %in% c(23, 24) ~ "Positive Control",
                TRUE ~ Compound),
            is_control = column %in% c(1, 2, 23, 24),
            plate_id = plate_id,
            master_plate_id = master_plate_id,
            dose_nM=dose_nM)
}

image_scores_CX5_100X <- readr::read_tsv("raw_data/plate_ids.tsv") %>%
    plyr::adply(1, function(df){
        cat("reading plate '", df$plate_id[1], "'\n", sep="")
        load_image_features(
            con=con,
            schema=df$schema[1],
            master_plate_id=df$master_plate_id[1],
            plate_id=df$plate_id[1],
            dose_nM=df$dose_nM[1])
    })

save(image_scores_CX5_100X, file="intermediate_data/image_scores_CX5_100X.Rdata")
image_scores_CX5_100X %>% arrow::write_parquet("product/image_scores_CX5_100X.parquet")








