

library(plyr)
library(tidyverse)
library(arrow)
library(hdf5r)
library(DBI)
library(RSQLite)

source("scripts/load_features_h5.R")

load("intermediate_data/plate_map_999A.Rdata")
load("intermediate_data/plate_map_1999B.Rdata")
load("intermediate_data/plate_map_2020A.Rdata")
load("intermediate_data/plate_map_2021A.Rdata")

system("
cd raw_data
mkdir infected_patches
cd infected_patches
aws s3 cp s3://sextoncov19/CPOutput_DR_Covid19.tar .
tar xvf CPOutput_DR_Covid19.tar
")

# input data
h5_dataset <- hdf5r::H5File$new("raw_data/infected_patches/CPOutput/DefaultOUT.h5", "r")
# h5_dataset$ls(recursive = TRUE)


# output path
output_path <- "intermediate_data/infected_patch_1999B_2020A_2021A_20201017"
if (!dir.exists(paths = output_path)) {
    dir.create(
        path = output_path,
        recursive = TRUE)
}



###################################
# Load features from HDF5 dataset #
###################################

# h5_dataset[["Measurements/2020-10-15-15-10-50/Image"]][['URL_NP']][["data"]][]
# 10368 rows like "file:/home/ubuntu/Desktop/Jonny/20200512T002847_2021A/Projection/W0384F0007T0001Z000C4.tif"
# but it looks like maybe two of the images failed to load?
# anyway, all we need is the ImageNumber to link the objects
# and the URL_NP has the well_id in it

image_features <- h5_dataset %>%
    load_features_h5(
        table_name = "Measurements/2020-10-15-15-10-50/Image",
        include_features = c("ImageNumber", "URL_NP")) %>%
    dplyr::mutate(
        plate_id = URL_NP %>%
            stringr::str_extract("[_].+Projection") %>%
            stringr::str_replace("_", "") %>%
            stringr::str_replace("/Projection", ""),
        well_id = URL_NP %>% stringr::str_extract("W[0-9][0-9][0-9][0-9]") %>%
            stringr::str_replace("W", ""),
        field_id = URL_NP %>% stringr::str_extract("F[0-9][0-9][0-9][0-9]") %>%
            stringr::str_replace("F", ""),
        row = floor((as.numeric(well_id) - 1) / 24) + 1,
        column = ((as.numeric(well_id) - 1) %% 24) + 1) %>%
    dplyr::select(-URL_NP)

viral_features <- h5_dataset %>% load_features_h5(
    table_name = "Measurements/2020-10-15-15-10-50/ViralObj",
    exclude_features = c("ObjectNumber"))
            
nuclei_features <- h5_dataset %>% load_features_h5(
    table_name = "Measurements/2020-10-15-15-10-50/Nuclei",
    exclude_features = c("ObjectNumber"))

syn_nuc_features <- h5_dataset %>% load_features_h5(
    table_name = "Measurements/2020-10-15-15-10-50/syn_nucs",
    exclude_features = c("ObjectNumber"))

puncta_features <- h5_dataset %>% load_features_h5(
    table_name = "Measurements/2020-10-15-15-10-50/puncta",
    exclude_features = c("ObjectNumber"))


###################
# join plate-maps #
###################
image_features <- image_features %>%
    dplyr::left_join(
        dplyr::bind_rows(
            plate_map_1999B %>%
                dplyr::transmute(
                    plate_id, row, column,
                    condition = Condition,
                    drug_1 = "Remdesivir",
                    drug_1_units = "µM",
                    drug_1_concentration = Remdesivir_Concentration,
                    drug_1_label = rem_label,
                    drug_2 = "Lactoferrin",
                    drug_2_units = "µg/mL",
                    drug_2_concentration = Lactoferrin_Concentration,
                    drug_2_label = lf_label),
            plate_map_2020A %>%
                dplyr::rename(condition = tidyselect::any_of("Fluid name")) %>%
                dplyr::transmute(
                    plate_id, row, column,
                    condition,
                    drug_1 = "Remdesivir",
                    drug_1_units = "µM",
                    drug_1_concentration = Remdesivir_Concentration,
                    drug_1_label = rem_label,
                    drug_2 = "Lactoferrin",
                    drug_2_units = "µg/mL",
                    drug_2_concentration = Lactoferrin_Concentration,
                    drug_2_label = lf_label),
            plate_map_2021A %>%
                dplyr::transmute(
                    plate_id, row, column,
                    condition = Condition,
                    drug_1 = "Hydroxychloroquine",
                    drug_1_units = "µM",
                    drug_1_concentration = Hydroxychloroquine_Concentration,
                    drug_1_label = hcq_label,
                    drug_2 = "Lactoferrin",
                    drug_2_units = "µg/mL",
                    drug_2_concentration = Lactoferrin_Concentration,
                    drug_2_label = lf_label)),
            by = c("plate_id", "row", "column"))

viral_features <- viral_features %>%
    dplyr::left_join(
        image_features,
        by = c("ImageNumber"))
            
nuclei_features <- nuclei_features %>%
    dplyr::left_join(
        image_features,
        by = c("ImageNumber"))
            
syn_nuc_features <- syn_nuc_features %>%
    dplyr::left_join(
        image_features,
        by = c("ImageNumber"))

puncta_features <- puncta_features %>%
    dplyr::left_join(
        image_features,
        by = c("ImageNumber"))


#
# Save features
#
viral_features %>%
    arrow::write_parquet(
        sink = paste0(output_path, "/viral_features.parquet"))
nuclei_features %>%
    arrow::write_parquet(
        sink = paste0(output_path, "/nuclei_features.parquet"))
syn_nuc_features %>%
    arrow::write_parquet(
        sink = paste0(output_path, "/syn_nuc_features.parquet"))
puncta_features %>%
    arrow::write_parquet(
        sink = paste0(output_path, "/puncta_features.parquet"))



#########################################
# identify feature and metadata columns #
#########################################
viral_feature_columns <- viral_features %>%
    dplyr::select(
        tidyselect::starts_with("AreaShape"),
        tidyselect::starts_with("Intensity"),
        Children_Nuclei_Count,
        Mean_Nuclei_Distance_Centroid_ViralObj,
        tidyselect::starts_with("RadialDistribution"),
        tidyselect::starts_with("Texture"),
        # remove plate location features
        -tidyselect::matches("_Center_"),
        # remove low-variance features
        -tidyselect::matches("RadialDistribution_ZernikePhase_.*_0"),
        -Mean_Nuclei_Distance_Centroid_ViralObj,
        -AreaShape_EulerNumber,
        -RadialDistribution_ZernikePhase_NP_0_0) %>%
    names() %>%
    tibble::tibble(feature = .) %>%
    tidyr::separate(
        col = feature,
        into = c("measure_type", "measure", "channel"),
        sep = "_",
        remove = FALSE,
        extra = "merge") %>%
    dplyr::mutate(
        transform = dplyr::case_when(
           # these should be log-transformed, but they have values <= 0 so use log1p
           measure_type == "Intensity" ~ "log1p",
           measure_type == "Texture" & measure %>% stringr::str_detect("Variance") ~ "log1p",
           measure_type == "Texture" & measure == "SumAverage" ~ "log",
           measure_type == "Texture" & measure == "Contrast" ~ "log1p",
           measure_type == "RadialDistribution" & measure_type == "ZernikeMagnitude" ~ "log",
           measure_type == "AreaShape" & measure == "Area" ~ "log",
           TRUE ~ "identity"))
viral_metadata_columns <- viral_features %>%
    dplyr::select(
        -tidyselect::any_of(viral_feature_columns$feature)) %>%
    names() %>%
    tibble::tibble(feature = .)

nuclei_feature_columns <- nuclei_features %>%
    dplyr::select(
        tidyselect::starts_with("AreaShape"),
        tidyselect::starts_with("Intensity"),
        Children_syn_nucs_Count,
        Distance_Centroid_ViralObj,
        -tidyselect::matches("_Center_")) %>%
    names() %>%
    tibble::tibble(feature = ., transformation = "identity")
nuclei_metadata_columns <- nuclei_features %>%
    dplyr::select(
        -tidyselect::any_of(nuclei_feature_columns$feature)) %>%
    names() %>%
    tibble::tibble(feature = .)


syn_nuc_feature_columns <- syn_nuc_features %>%
    dplyr::select(
        tidyselect::starts_with("AreaShape"),
        tidyselect::starts_with("Intensity"),
        -tidyselect::matches("_Center_")) %>%
    names() %>%
    tibble::tibble(feature = ., transformation = "identity")
syn_nuc_metadata_columns <- syn_nuc_features %>%
    dplyr::select(
        -tidyselect::any_of(syn_nuc_feature_columns$feature)) %>%
    names() %>%
    tibble::tibble(feature = .)

# puncta don't have any annotated morphological features 
puncta_metadata_columns <- puncta_features %>%
    names() %>%
    tibble::tibble(feature = .)

#
# save feature and metadata columns
#
viral_feature_columns %>% readr::write_tsv(
    path = paste0(output_path, "/viral_feature_columns.tsv"))
viral_feature_columns %>%
    dplyr::mutate(transform = 'identity') %>%
    readr::write_tsv(path = paste0(output_path, "/viral_feature_no_transform_columns.tsv"))
viral_metadata_columns %>% readr::write_tsv(
    path = paste0(output_path, "/viral_metadata_columns.tsv"))

nuclei_feature_columns %>% readr::write_tsv(
    path = paste0(output_path, "/nuclei_feature_columns.tsv"))
nuclei_metadata_columns %>% readr::write_tsv(
    path = paste0(output_path, "/nuclei_metadata_columns.tsv"))

syn_nuc_feature_columns %>% readr::write_tsv(
    path = paste0(output_path, "/syn_nuc_feature_columns.tsv"))
syn_nuc_metadata_columns %>% readr::write_tsv(
    path = paste0(output_path, "/syn_nuc_metadata_columns.tsv"))

puncta_metadata_columns %>% readr::write_tsv(
    path = paste0(output_path, "/puncta_metadata_columns.tsv"))


###################################
# viral patch data for plate 999A #
###################################

input_path <- "raw_data/infected_patches"
if (!dir.exists(paths = input_path)) {
    dir.create(
        path = input_path,
        recursive = TRUE)
}

# output path
output_path <- "intermediate_data/infected_patch_999A_20201112"
if (!dir.exists(paths = output_path)) {
    dir.create(
        path = output_path,
        recursive = TRUE)
}


# for some reason the SQL folder isn't working through the S3-fuse
system("
python $(which aws) s3 cp s3://sextoncov19/SQL/SARS_ViralPheno.sql raw_data/infected_patches/
")

# covert mysql dump to sqlite3 to read in with dbplyr
# https://stackoverflow.com/questions/5164033/export-a-mysql-database-to-sqlite-database
system("
pushd ~/opt
git clone https://github.com/dumblob/mysql2sqlite
popd
")

# this takes ~5-10 minutes
system("
time ~/opt/mysql2sqlite/mysql2sqlite raw_data/infected_patches/SARS_ViralPheno.sql | sqlite3 raw_data/infected_patches/SARS_ViralPheno.db3
")

con <- DBI::dbConnect(
    drv = RSQLite::SQLite(),
    dbname = "raw_data/infected_patches/SARS_ViralPheno.db3")

cat("Tables in 'raw_data/infected_patches/SARS_ViralPheno.db3':")
con %>% DBI::dbListTables() %>% print()

image_features <- con %>%
    dplyr::tbl("SARS_ViralPhenoPer_Image") %>%
    dplyr::select(
        ImageNumber,
        URL_NP = Image_URL_NP) %>%
    dplyr::collect(n=Inf) %>%
    dplyr::mutate(
        plate_id = URL_NP %>%
            stringr::str_extract("[_].+Projection") %>%
            stringr::str_replace("_", "") %>%
            stringr::str_replace("/Projection", ""),
        well_id = URL_NP %>% stringr::str_extract("W[0-9][0-9][0-9][0-9]") %>%
            stringr::str_replace("W", ""),
        field_id = URL_NP %>% stringr::str_extract("F[0-9][0-9][0-9][0-9]") %>%
            stringr::str_replace("F", ""),
        row = floor((as.numeric(well_id) - 1) / 24) + 1,
        column = ((as.numeric(well_id) - 1) %% 24) + 1) %>%
    dplyr::select(-URL_NP)
               
viral_features <- con %>%
    dplyr::tbl("SARS_ViralPhenoPer_ViralObj") %>%
    dplyr::collect(n=Inf) %>%
    dplyr::rename_with(~stringr::str_replace(., "^ViralObj_", ""))

nuclei_features <- con %>%
    dplyr::tbl("SARS_ViralPhenoPer_Nuclei") %>%
    dplyr::collect(n=Inf) %>%
    dplyr::rename_with(~stringr::str_replace(., "^Nuclei_", ""))

puncta_features <- con %>%
    dplyr::tbl("SARS_ViralPhenoPer_puncta") %>%
    dplyr::collect(n=Inf) %>%
    dplyr::rename_with(~stringr::str_replace(., "^puncta_", ""))

syn_nuc_features <- con %>%
    dplyr::tbl("SARS_ViralPhenoPer_syn_nucs") %>%
    dplyr::collect(n=Inf) %>%
    dplyr::rename_with(~stringr::str_replace(., "^syn_nucs_", ""))


#
# join plate-maps
#
image_features <- image_features %>%
    dplyr::left_join(
        plate_map_999A %>%
            dplyr::transmute(
                plate_id = Plate_ID,
                row, column,
                condition = COND,
                Compound, Concentration,
                drug = Compound,
                drug_units = Units,
                drug_concentration = Concentration,
                drug_label = Compound),
        by = c("plate_id", "row", "column"))

viral_features <- viral_features %>%
    dplyr::left_join(
        image_features,
        by = c("ImageNumber"))
            
nuclei_features <- nuclei_features %>%
    dplyr::left_join(
        image_features,
        by = c("ImageNumber"))
            
syn_nuc_features <- syn_nuc_features %>%
    dplyr::left_join(
        image_features,
        by = c("ImageNumber"))

puncta_features <- puncta_features %>%
    dplyr::left_join(
        image_features,
        by = c("ImageNumber"))


#
# Save features
#
viral_features %>%
    arrow::write_parquet(
        sink = paste0(output_path, "/viral_features.parquet"))
nuclei_features %>%
    arrow::write_parquet(
        sink = paste0(output_path, "/nuclei_features.parquet"))
syn_nuc_features %>%
    arrow::write_parquet(
        sink = paste0(output_path, "/syn_nuc_features.parquet"))
puncta_features %>%
    arrow::write_parquet(
        sink = paste0(output_path, "/puncta_features.parquet"))


# clean up raw input data
system("
rm raw_data/infected_patches/SARS_ViralPheno.db3
rm raw_data/infected_patches/SARS_ViralPheno.sql
")
