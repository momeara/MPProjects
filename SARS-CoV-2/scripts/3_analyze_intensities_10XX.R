
library(plyr)
library(tidyverse)
library(RMySQL)
library(magrittr)
library(tictoc)
library(arrow)

load("intermediate_data/drug_plate1_image_scores.Rdata")

con <- DBI::dbConnect(
    RMySQL::MySQL(),
    host="covid19cp.cgymeokijgns.us-east-1.rds.amazonaws.com",
    port=3306,
    user = "covid19cp",
    password = "Genes-brett-Flip-9Bottling")
con %>% DBI::dbSendQuery("USE covid19cp")

get_cells <- function(plate_prefix, parent_plate_barcode){
    cat("Getting Cell features for ", plate_prefix, "\n", sep="")
    features_tbl <- con %T>%
        {tictoc::tic("collect")} %>%
        dplyr::tbl(paste0(plate_prefix, "_Per_Cells")) %>%
        dplyr::collect(n=Inf) %T>%
        {tictoc::toc()} %>%
        dplyr::mutate(
             parent_plate_barcode = parent_plate_barcode,
             plate_prefix = plate_prefix)      

    features_tbl <- features_tbl %>%
        dplyr::left_join(
           drug_plate1_image_scores,
           by=c("parent_plate_barcode", "plate_prefix", "ImageNumber"))
}


get_features <- function(object_type, plate_prefix, parent_plate_barcode){
    n_row <- con %>% dplyr::tbl(paste0(plate_prefix, "_Per_", object_type)) %>% dplyr::tally()
    cat("Getting ", n_row$n[1], " ", object_type, " features for ", plate_prefix, "\n", sep="")
    features_tbl <- con %T>%
        {tictoc::tic("collect")} %>%
        dplyr::tbl(paste0(plate_prefix, "_Per_", object_type)) %>%
        dplyr::collect(n=Inf) %T>%
        {tictoc::toc()} %>%
        dplyr::mutate(
             parent_plate_barcode = parent_plate_barcode,
             plate_prefix = plate_prefix)      

    cat("  merging image meta data...\n")
    features_tbl <- features_tbl %>%
        dplyr::left_join(
           drug_plate1_image_scores,
           by=c("parent_plate_barcode", "plate_prefix", "ImageNumber"))
}


features_tbl <- arrow::read_parquet("intermediate_data/COVID_Plate2_Per_Cells.parquet")

##    
p <- ggplot2::ggplot(data=features_tbl) +
    ggplot2::theme_bw() +
    ggplot2::geom_histogram(
        mapping=ggplot2::aes(
            x=Cells_Intensity_IntegratedIntensity_Virus,
            y=..count..),                   
        bins=80) +
    facet_grid(row ~ column) +
    scale_y_log10("Well Cell Count") +
    scale_x_log10("Integrated Virus Intensity", limits=c(1,1000))

ggplot2::ggsave(
    "product/figures/drug_plate2_virus_intensity_by_well_histogram_200420.pdf",
    width=30, height=30)


p <- ggplot2::ggplot(data=features_tbl) +
    ggplot2::theme_bw() +
    ggplot2::geom_histogram(
        mapping=ggplot2::aes(
            x=Cells_Intensity_IntegratedIntensity_Virus,
            y=..count..),                   
        bins=130) +
    scale_y_log10("Well Cell Count") +
    scale_x_log10("Integrated Virus Intensity", limits=c(1,1000))

ggplot2::ggsave(
    "product/figures/drug_plate2_virus_intensity_histogram_200420.pdf",
    width=8, height=4)


###

##    
p <- ggplot2::ggplot(data=features_tbl) +
    ggplot2::theme_bw() +
    ggplot2::geom_histogram(
        mapping=ggplot2::aes(
            x=Cells_Intensity_MeanIntensity_Virus,
            y=..count..),                   
        bins=80) +
    facet_grid(row ~ column) +
    scale_y_log10("Well Cell Count") +
    scale_x_log10("Mean Virus Intensity")

ggplot2::ggsave(
    "product/figures/drug_plate2_virus_mean_intensity_by_well_histogram_200420.pdf",
    width=30, height=30)


########################################3
drug_plate1_cytoplasm_features <- dplyr::bind_rows(
    get_features(object_type="Cytoplasm", plate_prefix="COVID_Plate1_2", parent_plate_barcode=1001),
    get_features(object_type="Cytoplasm", plate_prefix="COVID_Plate2", parent_plate_barcode=1001),
    get_features(object_type="Cytoplasm", plate_prefix="COVID_Plate3", parent_plate_barcode=1001),
    get_features(object_type="Cytoplasm", plate_prefix="COVID_Plate4", parent_plate_barcode=1001),
    get_features(object_type="Cytoplasm", plate_prefix="COVID_Plate5", parent_plate_barcode=1001))

drug_plate1_cytoplasm_features <- drug_plate1_cytoplasm_features %>%
    dplyr::select(-tidyselect::contains("_Location_"))

drug_plate1_cytoplasm_features %>% arrow::write_parquet("intermediate_data/drug_plate1_cytoplasm_features.parquet")


p <- ggplot2::ggplot(data=drug_plate1_cytoplasm_features) +
    ggplot2::theme_bw() +
    ggplot2::geom_histogram(
        mapping=ggplot2::aes(
            x=Cytoplasm_Intensity_MeanIntensity_Virus,
            y=..count..),                   
        bins=80) +
    ggplot2::facet_wrap(~plate_prefix, ncol=1) +
    ggplot2::ggtitle("Cytoplasm_Intensity_MeanIntensity_Virus feature distribution") +
    ggplot2::scale_y_log10("Well Cell Count") +
    ggplot2::scale_x_log10("Mean Virus Intensity")

ggplot2::ggsave(
    "product/figures/drug_plate2_cytoplasm_virus_mean_intensity_by_dose_histogram_200420.pdf",
    width=7, height=12)

##
p <- ggplot2::ggplot(data=drug_plate1_cytoplasm_features %>% dplyr::filter(plate_prefix == "COVID_Plate2")) +
    ggplot2::theme_bw() +
    ggplot2::geom_histogram(
        mapping=ggplot2::aes(
            x=Cytoplasm_Intensity_MeanIntensity_Virus,
            y=..count..),                   
        bins=80) +
    ggplot2::facet_wrap(row~column) +
    facet_grid(row ~ column) +    
    ggplot2::ggtitle("Cytoplasm_Intensity_MeanIntensity_Virus feature distribution") +
    ggplot2::scale_y_log10("Well Cell Count") +
    ggplot2::scale_x_log10("Mean Virus Intensity")

ggplot2::ggsave(
    "product/figures/drug_plate2_cytoplasm_virus_mean_intensity_COVID_Plate2_by_well_histogram_200420.pdf",
    width=30, height=30)
