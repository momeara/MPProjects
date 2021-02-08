
library(plyr)
library(tidyverse)
library(fuzzyjoin)
library(ggplot2)
library(readxl)
library(MPStats)
library(arrow)

source("scripts/database.R")

plate_ids <- readr::read_tsv("raw_data/plate_ids.tsv")
cell_feature_columns <- readr::read_tsv("raw_data/cell_feature_columns.tsv")

con <- get_primary_database_connection("covid19cq1")
plate_id <- "2019A"
image_table_name <- paste0("SARS_", plate_id, "_Per_Image")
image_features <- con %>%
    dplyr::tbl(image_table_name) %>%
    dplyr::collect(n = Inf)
image_width <- image_features$Image_Width_Virus[1]
image_height <- image_features$Image_Height_Virus[1]


cell_features <- c(
    '2006A', '2007A', '2008A', '2009A',
    '2010A', '2010A',          '2012A',
    '2013A', '2014A', '2015A', '2016A',
    '2017A',          '2019A') %>%
    plyr::ldply(function(plate_id) {
        cat("well location batch effects for plate ", plate_id, " ...\n", sep = "")
        cell_features <- arrow::read_parquet(
            paste0("product/covid19cq1_SARS_", plate_id, "_Cell_MasterDataTable.parquet")) %>%
            dplyr::mutate_at(cell_feature_columns$feature, ~ scale(.) %>% as.vector) %>%
            dplyr::mutate(
                field_id = as.numeric(Image_Metadata_FieldID),
                field_x = field_id %% 3,
                field_y = floor(field_id / 3),
                location_x = Nuclei_Location_Center_X / image_width + (field_x - 3 / 2),
                location_y = Nuclei_Location_Center_Y / image_height + (field_y - 3 / 2),
                radial_length = sqrt(location_x**2 + location_y**2),
                angle = atan2(location_x, location_y))
    })


fits <- plyr::ddply("plate_id", function(df){
    cell_feature_columns$feature %>%
        plyr::ldply(function(column){
            cat("fitting model for feature ", column, "...\n")
            lm(
                formula = paste0(column, " ~ 0 + location_x + location_y") %>% as.formula,
                data = df) %>%
                broom::tidy() %>%
                dplyr::mutate(feature = !!column)
        })
})



most_skewed <- fits %>%
    dplyr::arrange(desc(estimate)) %>%
    head(6)

well_features <- cell_features %>%
    dplyr::select(
        plate_id,
        row,
        column,
        tidyselect::one_of(most_skewed$feature),
        location_x,
        location_y) %>%
    dplyr::group_by(plate_id, row, column) %>%
    dplyr::summarize_at(
        c(most_skewed$feature, "location_x", "location_y"),
        mean) %>%
    dplyr::ungroup()

well_features_long <- well_features %>%
    tidyr::pivot_longer(
        cols = tidyselect::one_of(most_skewed$feature),
        names_to = "feature",
        values_to = "feature_value") %>%
    tidyr::pivot_longer(
        cols = c(location_x, location_y),
        names_to = "term",
        values_to = "term_value") %>%
    dplyr::semi_join(
        most_skewed,
        by = c("plate_id", "term", "feature"))

ggplot2::ggplot(
    data = well_features_long) +
    ggplot2::theme_bw() +
    ggplot2::geom_point(
        mapping = ggplot2::aes(
            x = term_value,
            y = feature_value)) +
    ggplot2::facet_wrap(
        ~ plate_id + feature + term)

ggsave(
    "product/figures/batch_effects/field_batch_effects/most_skewed_20200627.pdf",
    width = 8,
    height = 8)
