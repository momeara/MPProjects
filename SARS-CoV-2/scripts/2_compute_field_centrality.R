library(plyr)
library(tidyverse)
library(arrow)
library(monocle3)


cell_metadata_columns <- readr::read_tsv("raw_data/cell_metadata_columns_TS.tsv")

cell_feature_columns <- readr::read_tsv("raw_data/cell_feature_columns_TS.tsv")

plate_ids <- readr::read_tsv("raw_data/plate_ids.tsv") %>%
    dplyr::filter(plate_id %>% stringr::str_detect("^TS")) %>%
    magrittr::extract2("plate_id")


cell_metadata_dataframe <- plate_ids %>%
    plyr::ldply(function(plate_id){
        cat("Loading metadata for plate '", plate_id, "'\n", sep = "")
        arrow::read_parquet(
            paste0("product/covid19cq1_SARS_", plate_id, "_Cell_MasterDataTable.parquet"),
            col_select = c(
                cell_metadata_columns$feature,
                Nuclei_Number_Object_Number,
                Nuclei_Location_Center_X,
                Nuclei_Location_Center_Y))
    }) %>%
    dplyr::mutate(
        plate_id = factor(
            plate_id,
            levels=plate_ids,
            labels=plate_ids %>%
                stringr::str_replace("TS", "")))


cell_metadata_dataframe %>%
    dplyr::count(plate_id) %>%
    dplyr::arrange(plate_id)

# load cell features as a matrix with dimensions [feature, cell]
cell_features_dataframe <- plate_ids %>%
    plyr::ldply(function(plate_id){
        cat("Loading features for plate '", plate_id, "'\n", sep = "")
        arrow::read_parquet(
            paste0("product/covid19cq1_SARS_", plate_id, "_Cell_MasterDataTable.parquet"),
            col_select = cell_feature_columns$feature)
    })



infected_cells <- cell_data_set %>% monocle3::choose_cells()

library(BiocNeighbors)

cell_neighbors <- cell_metadata_dataframe %>%
    dplyr::filter(
        plate_id == "3h",
        Image_Metadata_WellID == "0001",
        Image_Metadata_FieldID == "0001") %>%
    dplyr::select(
        Nuclei_Location_Center_X,
        Nuclei_Location_Center_Y) %>%
    BiocNeighbors::findKNN(k = 8)

neighbor_graph <- cell_neighbors$distance %>% as("dgTMatrix")
neighbor_graph <- data.frame(
    from=neighbor_graph$i,
    to=neighbor_graph$x,
    value = neighbor_graph$j)
    
distance_graph <- cell_metadata_dataframe %>%
    dplyr::select(
        plate_id,
        Image_Metadata_WellID,
        Image_Metadata_FieldID,
        Nuclei_Location_Center_X,
        Nuclei_Location_Center_Y) %>%
    plyr::ddply(
        c("plate_id", "Image_Metadata_WellID", "Image_Metadata_FieldID"),
        function(field) {
            if (nrow(field) < 50) {
                return(data.frame())
            }
            cat("plate: ", as.character(field$plate_id[1]),
                " Well: ", field$Image_Metadata_WellID[1],
                " Field: ", field$Image_Metadata_FieldID[1], "\n",
                sep = "")
            knn_neighbors <- field %>%
                dplyr::select(
                    Nuclei_Location_Center_X,
                    Nuclei_Location_Center_Y) %>%            
                BiocNeighbors::findKNN(k = 8)
            index <- as(knn_neighbors$index, "dgTMatrix")
            distance <- as(knn_neighbors$distance, "dgTMatrix")
            data.frame(
                from = index@i,
                to = index@x,
                value = distance@x)
            })
              
   
distance_graph <- distance_graph %>%
    dplyr::left_join(
        
