
library(plyr)
library(tidyverse)
library(magrittr)
library(arrow)
library(MPStats)

source("parameters.R")

source("scripts/load_features_from_database.R")

load("intermediate_data/image_scores_CQ1_TS.Rdata")
load("intermediate_data/image_scores_CQ1_TS_202008.Rdata")

con <- MPStats::get_database_connection(
    database_parameters =  parameters$databases$covid19cq1)

##########################
# Time Series 1 (202006) #
##########################

plate_ids <- c("TS3h", "TS6h", "TS12h", "TS18h", "TS24h", "TS48h")

# collect Cell features
plate_ids %>% plyr::l_ply(function(plate_id){
    object_features <- 
        dplyr::bind_cols(
             con %>% collect_object_features_from_database(
                  table_name = paste0("SARS_", plate_id, "_Per_Cell")),
             con %>% collect_object_features_from_database(
                  table_name = paste0("SARS_", plate_id, "_Per_Nuclei")) %>%
                  dplyr::select(-ImageNumber),
             con %>% collect_object_features_from_database(
                  table_name = paste0("SARS_", plate_id, "_Per_Cytoplasm")) %>%
                  dplyr::select(-ImageNumber)) 
    
    image_scores_CQ1_TS %>%
        dplyr::filter(plate_id == !!plate_id) %>%
        dplyr::inner_join(object_features, by = "ImageNumber") %>%
        arrow::write_parquet(
            sink = paste0("product/", plate_id, "_Cell_MasterDataTable.parquet")
})

# collect InfectedCells features
plate_ids %>% plyr::l_ply(function(plate_id){
    object_features <- 
        con %>% collect_object_features_from_database(
          table_name = paste0("SARS_", plate_id, "_Per_InfectedCell")))

    image_scores_CQ1_TS %>%
        dplyr::filter(plate_id == !!plate_id) %>%
        dplyr::inner_join(object_features, by = "ImageNumber") %>%
        arrow::write_parquet(
            sink = paste0("product/", plate_id, "_InfectedCells_MasterDataTable.parquet")
})


##########################
# Time Series 2 (202008) #
##########################

plate_ids <- c("TS2PL1", "TS2PL2", "TS2PL3")

# collect Cell features
plate_ids %>% plyr::l_ply(function(plate_id){
    object_features <- 
        dplyr::bind_cols(
             con %>% collect_object_features_from_database(
                  table_name = paste0("SARS_", plate_id, "_Per_Cell")),
             con %>% collect_object_features_from_database(
                  table_name = paste0("SARS_", plate_id, "_Per_Nuclei")) %>%
                  dplyr::select(-ImageNumber),
             con %>% collect_object_features_from_database(
                  table_name = paste0("SARS_", plate_id, "_Per_Cytoplasm")) %>%
                  dplyr::select(-ImageNumber)) 
    
    image_scores_CQ1_TS_202008 %>%
        dplyr::filter(plate_id == !!plate_id) %>%
        dplyr::inner_join(object_features, by = "ImageNumber") %>%
        arrow::write_parquet(
            sink = paste0("product/", plate_id, "_Cell_MasterDataTable.parquet")
})

# collect InfectedCells features
plate_ids %>% plyr::l_ply(function(plate_id){
    object_features <- 
        con %>% collect_object_features_from_database(
          table_name = paste0("SARS_", plate_id, "_Per_InfectedCell")))

    image_scores_CQ1_TS_202008 %>%
        dplyr::filter(plate_id == !!plate_id) %>%
        dplyr::inner_join(object_features, by = "ImageNumber") %>%
        arrow::write_parquet(
            sink = paste0("product/", plate_id, "_InfectedCells_MasterDataTable.parquet")
})
