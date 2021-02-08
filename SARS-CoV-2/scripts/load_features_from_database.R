
library(plyr)
library(tidyverse)
library(magrittr)
#library(argparser)

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

collect_cell_features <- function(con, schema, image_scores, plate_id) {
    con %>% set_schema(schema)
    cat(
        "Gathering cell/nuclei/cytoplasm features for plate '",
        plate_id, "' from schema '", schema, "'\n", sep = "")
    output_path <- paste0("product/", schema, "_SARS_", plate_id, "_Cell_MasterDataTable.parquet")
    cat("Writing output to '", output_path, "'\n", sep = "")
    table_name <- paste0("SARS_", plate_id, "_Per_Cells")
    cell_features <- con %T>%
        {tictoc::tic(paste0("collect ", table_name))} %>%
        dplyr::tbl(table_name) %>%
        dplyr::collect(n=Inf) %T>%
        {tictoc::toc()}
    cat("   nrow: ", cell_features %>% nrow(), ",  ncol: ", cell_features %>% ncol(), "\n", sep="")
    table_name <- paste0("SARS_", plate_id, "_Per_Nuclei")
    nuclei_features <- con %T>%
        {tictoc::tic(paste0("collect ", table_name))} %>%
        dplyr::tbl(table_name) %>%
        dplyr::collect(n=Inf) %T>%
        {tictoc::toc()}
    cat("   nrow: ", nuclei_features %>% nrow(), ",  ncol: ", nuclei_features %>% ncol(), "\n", sep="")
    table_name <- paste0("SARS_", plate_id, "_Per_Cytoplasm")
    cytoplasm_features <- con %T>%
        {tictoc::tic(paste0("collect ", table_name))} %>%        
        dplyr::tbl(table_name) %>%
        dplyr::collect(n=Inf) %T>%
        {tictoc::toc()}
    cat("   nrow: ", cytoplasm_features %>% nrow(), ",  ncol: ", cytoplasm_features %>% ncol(), "\n", sep="")
    cell_features <- image_scores %>%
        dplyr::filter(plate_id == !!plate_id) %>%
        dplyr::inner_join(
            dplyr::bind_cols(
                cell_features,
                nuclei_features %>% dplyr::select(-ImageNumber),
                cytoplasm_features %>% dplyr::select(-ImageNumber)),
            by="ImageNumber")
    cell_features %>%
        arrow::write_parquet(output_path)
}


collect_InfectedCells_features <- function(con, schema, plate_id, image_scores) {
    con %>% set_schema(schema)
    cat("Gathering InfectedCells features for plate '", plate_id, "' in schema '", schema, "'\n", sep="")
    table_name <- paste0("SARS_", plate_id, "_Per_InfectedCells")
    InfectedCells_features <- con %T>%
        {tictoc::tic(paste0("collect ", table_name))} %>%
        dplyr::tbl(table_name) %>%
        dplyr::collect(n = Inf) %T>%
        {tictoc::toc()}
    cat("   nrow: ", InfectedCells_features %>% nrow(), ",  ncol: ", InfectedCells_features %>% ncol(), "\n", sep = "")
    InfectedCells_features <- image_scores %>%
        dplyr::filter(plate_id == !!plate_id) %>%
        dplyr::inner_join(
            InfectedCells_features,
            by = "ImageNumber")
    InfectedCells_features %>%
        arrow::write_parquet(paste0("product/", schema, "_SARS_", plate_id, "_InfectedCells_MasterDataTable.parquet"))
}


collect_syncytia_features <- function(con, schema, plate_id, image_scores){
    con %>% set_schema(schema)
    cat("Gathering syncytia features for plate '", plate_id, "' in schema '", schema, "'\n", sep="")
    table_name <- paste0("SARS_", plate_id, "_Per_syncytia")
    syncytia_features <- con %T>%
        {tictoc::tic(paste0("collect ", table_name))} %>%        
        dplyr::tbl(table_name) %>%
        dplyr::collect(n=Inf) %T>%
        {tictoc::toc()}
    cat("   nrow: ", syncytia_features %>% nrow(), ",  ncol: ", syncytia_features %>% ncol(), "\n", sep="")
    syncytia_features <- image_scores %>%
        dplyr::filter(plate_id == !!plate_id) %>%
        dplyr::inner_join(
            syncytia_features,
            by="ImageNumber")
    syncytia_features %>% arrow::write_parquet(paste0("product/SARS_", plate_id, "_Syncytia_MasterDataTable.parquet"))
}

collect_nucleoli_features <- function(con, schema, plate_id, image_scores){
    con %>% set_schema(schema)
    cat("Gathering nucleoli features for plate '", plate_id, "' in schema '", schema, "'\n", sep="")
    table_name <- paste0("SARS_", plate_id, "_Per_Nucleoli")
    nucleoli_features <- con %T>%
        {tictoc::tic(paste0("collect ", table_name))} %>%        
        dplyr::tbl(table_name) %>%
        dplyr::collect(n=Inf) %T>%
        {tictoc::toc()}
    cat("   nrow: ", nucleoli_features %>% nrow(), ",  ncol: ", nucleoli_features %>% ncol(), "\n", sep="")
    nucleoli_features <- image_scores %>%
        dplyr::filter(plate_id == !!plate_id) %>%
        dplyr::inner_join(
            nucleoli_features,
            by="ImageNumber")
    nucleoli_features %>% arrow::write_parquet(paste0("product/SARS_", plate_id, "_Nucleoli_MasterDataTable.parquet"))
}
