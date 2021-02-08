
library(plyr)
library(tidyverse)
library(hdf5r)

source("scripts/mount_S3_bucket.R")
mount_S3_bucket()

source("scripts/load_dataset.R")
dataset_id <- "02f50f3f-4ec9-4192-a1f4-06ea93d50922"
if (!dir.exists(paste0("raw_data/", dataset_id))) {
    load_dataset(
        dataset_id = "02f50f3f-4ec9-4192-a1f4-06ea93d50922",
        verbose = TRUE)
}
dataset <- hdf5r::H5File$new(paste0("raw_data/", dataset_id, "/cpdata.h5"), "r")
#dataset$ls(recursive = TRUE)

plate_barcode <- dataset[['meta/experiment/assay_plate_barcode']][]
# "BR00110036"

object_columns <- tibble::tibble(feature=dataset[["object/cFeatureName"]][])
object_features <- dataset[["object/M"]][,]
colnames(object_features) <- object_columns$feature
object_features <- object_features %>% tibble::as_tibble()

cell_features <- object_features %>%
   dplyr::filter(
        Cell_Children_Cytoplasm_Count > 0,
        Nucleus_Children_Cells_Count > 0,
        Nucleus_Children_Cytoplasm_Count > 0)

# filter out cells without cytoplasms
# in this plate there are 14-24 per object that don't have children
cell_feature_columns <- object_columns %>%
    dplyr::filter(
        !(feature %in% c("Number_Object_Number")),
        !(feature %>% stringr::str_detect("_Parent_")),
        !(feature %>% stringr::str_detect("Center_[XY]")),
        !(feature %>% stringr::str_detect("_Count")))

cell_feature_columns <- cell_feature_columns %>%
    plyr::adply(1, function(df){
        feature <- df$feature[1]
        feature_values <- cell_features[feature]
        n_distinct <- feature_values %>% distinct() %>% nrow()
        n_na <- feature_values %>% dplyr::filter(is.na(!!sym(feature))) %>% nrow()
        n_nan <- feature_values %>% dplyr::filter(is.nan(!!sym(feature))) %>% nrow()
        return(data.frame(
            n_distinct = n_distinct,
            n_na = n_na,
            n_nan = n_nan,
            max_value = feature_values[[feature]] %>% max,
            min_value = feature_values[[feature]] %>% min))
    })

cell_feature_columns <- cell_feature_columns %>%
    dplyr::filter(
        n_distinct >= 50,
        n_nan <= 100)

cell_features <- object_features %>%
    dplyr::select(
        tidyselect::one_of(cell_feature_columns$feature))

cell_features %>%
    dplyr::filter_all(dplyr::any_vars(is.nan(.))) %>%
    nrow()

cell_features <- cell_features %>%
    tidyr::drop_na()

# check that it is not NaN
cell_features %>% sum()


cell_feature_columns <- cell_feature_columns %>%
    dplyr::mutate(
        object = feature %>% stringr::str_extract("^[A-Za-z]+"),
        feature_type = feature %>%
            stringr::str_match("^[A-Z][a-z]+_([A-Za-z]+)") %>%
            magrittr::extract(, 2))

cell_feature_columns <- cell_feature_columns %>%
    dplyr::mutate(transform = "identity")

cell_feature_columns %>%
    readr::write_tsv("raw_data/cell_feature_columns.tsv")


# remove cell index features
cell_feature_columns <- cell_features_columns %>%
    dplyr::filter(
        !(feature %>% stringr::str_detect("_Location_")),
        !(feature %>% stringr::str_detect("ClosestObjectNumber")))

cell_feature_columns %>%
    readr::write_tsv("raw_data/cell_feature_columns_20200814.tsv")
