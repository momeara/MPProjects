



library(plyr)
library(tidyverse)
library(fuzzyjoin)
library(ggplot2)
library(readxl)
library(MPStats)
library(arrow)


infected_ROI_membership <- arrow::read_parquet(
    "~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/top_hits_plate_scaled_200522a_umap2_2M_15_0.0/infected_ROI_membreship.parquet")

cell_features_20XX <- c(
    "2006A", "2007A", "2008A", "2009A",
    "2010A", "2010A",          "2012A",
    "2013A", "2014A", "2015A", "2016A",
    "2017A",          "2019A") %>%
    plyr::ldply(function(plate_id) {
        cat("Reading features from plate ", plate_id, " ...\n", sep = "")
        arrow::read_parquet(
            paste0("~/opt/MPLearn/vignettes/SARS-CoV-2/S25/raw_data/covid19cq1_SARS_", plate_id, "_plate_scaled_Cell_MasterDataTable.parquet"))
    })

cell_features_rem_lf <- arrow::read_parquet(
    "~/opt/MPLearn/vignettes/SARS-CoV-2/S25/raw_data/lf_rem_plate_scaled_1999B_2020A_Cell_MasterDataTable.parquet")

cell_features <- dplyr::bind_rows(
    cell_features_20XX,
    cell_features_rem_lf) %>%
    dplyr::bind_cols(
        infected_ROI_membership)

infected_cell_features <- cell_features %>%
    dplyr::filter(roi_0)
# nrow 2473367

# save/restore infected cell features for later
###############################################################33
infected_cell_features %>%
    arrow::write_parquet(
        "product/top_hits_infected_plate_scaled_200730_Cell_MasterDataTable.parquet")
infected_cell_features <-
    arrow::read_parquet(
        "product/top_hits_infected_plate_scaled_200730_Cell_MasterDataTable.parquet")
#################################################################

infected_cell_features <- infected_cell_features %>%
    dplyr::bind_cols(
        arrow::read_parquet(
            "~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/UMAP_embedding_top_hits_infected_plate_scaled_epochs=2000_re_embed_epochs=2000_200730/umap_embedding.parquet"),
        arrow::read_parquet(
            "~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/UMAP_embedding_top_hits_infected_plate_scaled_epochs=2000_re_embed_epochs=2000_200730/clusters.parquet"))



infected_cell_features %>%
    arrow::write_parquet(
        "product/top_hits_infected_plate_scaled_200730_Cell_MasterDataTable.parquet")


