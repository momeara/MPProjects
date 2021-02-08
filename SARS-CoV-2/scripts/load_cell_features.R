
library(plyr)
library(tidyverses)
library(magrittr)
library(tictoc)
library(arrow)

source("scripts/database.R")


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


collect_InfectedCells_features <- function(con, schema, plate_id, image_scores){
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
            by="ImageNumber")
    InfectedCells_features %>% arrow::write_parquet(paste0("product/SARS_", plate_id, "_InfectedCells_MasterDataTable.parquet"))
}
