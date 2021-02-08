
library(plyr)
library(tidyverse)
library(fuzzyjoin)
library(ggplot2)
library(readxl)
library(MPStats)
library(arrow)

plate_ids <- readr::read_tsv("raw_data/plate_ids.tsv")

cell_feature_columns <- readr::read_tsv("raw_data/cell_feature_columns.tsv")

top_hits <- readr::read_tsv("raw_data/top_hits_200522a.tsv")


infectivity_score_CQ1_features <- c(
    "Cells_Intensity_IntegratedIntensityEdge_Virus",
    "Cells_Intensity_MeanIntensityEdge_Virus",
    "Cells_Intensity_MaxIntensityEdge_Virus",
    "Cytoplasm_Intensity_StdIntensity_Virus")

add_infectivity_score_CQ1 <- function(cell_features) {
    cell_features %>%
        dplyr::mutate(
            infectivity_score = -0.59661 * Cells_Intensity_IntegratedIntensityEdge_Virus +
                3985.99907 * Cells_Intensity_MeanIntensityEdge_Virus +
                -54.14521 * Cells_Intensity_MaxIntensityEdge_Virus +
                -96.29369 * Cytoplasm_Intensity_StdIntensity_Virus)
}


c(
    '2006A', '2007A', '2008A', '2009A',
    '2010A', '2010A',          '2012A',
    '2013A', '2014A', '2015A', '2016A',
    '2017A',          '2019A') %>%
    plyr::l_ply(function(plate_id) {
        cat("Scaling plate ", plate_id, " ...\n", sep = "")
        arrow::read_parquet(
            paste0("product/covid19cq1_SARS_", plate_id, "_Cell_MasterDataTable.parquet")) %>%
            dplyr::mutate_at(cell_feature_columns$feature, ~ scale(.)[, 1]) %>%
            arrow::write_parquet(
                paste0("product/covid19cq1_SARS_", plate_id, "_plate_scaled_Cell_MasterDataTable.parquet"))
    })



cell_features <- top_hits %>%
    dplyr::distinct(plate_id) %>%
    dplyr::filter(plate_id != "1999B") %>%
    plyr::adply(1, function(df) {
        plate_id <- df$plate_id[1]
        cat("Reading cells from plate ", plate_id, "\n", sep = "")
        cell_features <- arrow::read_parquet(
            file = paste0("product/covid19cq1_SARS_", plate_id, "_Cell_MasterDataTable.parquet")) %>%
            add_infectivity_score_CQ1() %>%
            dplyr::mutate_at(cell_feature_columns$feature, scale) %>%
            dplyr::semi_join(
                dplyr::bind_rows(
                    top_hits,
                    data.frame(Compound = c("PC", "NC"))),
                by = "Compound")
    }) %>%
    dplyr::select(
        -COND,
        -Cell_Type,
        -Metadata_PlateID,
        -Concentration) %>%
    dplyr::mutate(
        master_plate_id = as.numeric(master_plate_id))

lf_rem_cell_features <- c("1999B", "2020A") %>%
    plyr::ldply(function(plate_id) {
        cat("Gathering cell features from plate ", plate_id, "\n", sep = "")
        lf_rem_cell_features <- arrow::read_parquet(
            file = paste0("product/covid19cq1_SARS_", plate_id, "_Cell_MasterDataTable.parquet")) %>%
            dplyr::select(
                -rem_label,
                -lf_label) %>%
            dplyr::mutate(
                master_plate_id = plate_id %>%
                    stringr::str_extract("^....") %>%
                    as.numeric(),
                Plate_Name = paste0("SARS_", plate_id),
                Plate_ID = NA) %>%
            add_infectivity_score_CQ1() %>%           
            dplyr::mutate_at(cell_feature_columns$feature, scale)
        if (plate_id == "1999B") {
            lf_rem_cell_features <- lf_rem_cell_features %>%
                dplyr::select(
                    -parent_plate_barcode)
        }
        if (plate_id == "2020A") {
            lf_rem_cell_features <- lf_rem_cell_features %>%
                dplyr::mutate(
                    Compound = `Fluid name`,
                    Condition = ifelse(Compound %in% c("PC", "NC"), Compound, "Treatment")) %>%
                dplyr::select(
#                    -Well_ID,
                    -`Fluid name`)
        }
        dplyr::bind_rows(
            lf_rem_cell_features %>%
                dplyr::filter(Remdesivir_Concentration == 0, Condition == "Treatment") %>%
                dplyr::select(-Remdesivir_Concentration) %>%
                dplyr::rename(
                    dose_nM=Lactoferrin_Concentration,
                    Units=Lactoferrin_Units) %>%
                dplyr::mutate(Compound = "Lactoferrin"),
            lf_rem_cell_features %>%
                dplyr::filter(Lactoferrin_Concentration == 0, Condition == "Treatment") %>%
                dplyr::select(-Lactoferrin_Concentration) %>%
                dplyr::rename(
                    dose_nM=Remdesivir_Concentration,
                    Units=`Remdesivir Units`) %>%
                dplyr::mutate(Compound = "Remdesivir"),
            lf_rem_cell_features %>%
                dplyr::filter(Compound %in% c("PC", "NC")) %>%
                dplyr::mutate(dose_nM = 0) %>%
                dplyr::select(
                    -Lactoferrin_Concentration,
                    -Remdesivir_Concentration,
                    -`Remdesivir Units`,
                    -Lactoferrin_Units))
    }) %>%
    dplyr::select(
        -Well_ID,
        -`Remdesivir Units`,
        -Lactoferrin_Units)


## just remove some PC cells?
#cell_features_pc <- cell_features %>% dplyr::filter(Compound == "PC")
#cell_features <- dplyr::bind_rows(
#    cell_features %>% dplyr::filter(Compound != "PC"),
#    cell_features %>%
#        dplyr::filter(Compound == "PC") %>%
#        dplyr::sample_n(421516))


## I messed up and forgot to map the dose_nM column for plate 999A
## so fix it here
#cell_features <- cell_features %>%
#    dplyr::left_join(
#        plate_map_999A %>%
#        dplyr::transmute(
#            plate_id = "0999A",
#            row,
#            column,
#            dose_nM_999A = dose_nM),
#        by = c("plate_id", "row", "column")) %>%
#    dplyr::mutate(
#        dose_nM = dplyr::case_when(
#            plate_id == "0999A" ~ dose_nM_999A,
#            TRUE ~ dose_nM)) %>%
#    dplyr::select(-dose_nM_999A)
#
#
#cell_features %>%
#    arrow::write_parquet("product/top_hit_cells_plate_scaled_200519.parquet")


############################################3

# check cell features from these plates can be joined

z_20XX <- data.frame(feature=cell_features %>% names())
z_lf <- data.frame(feature=lf_rem_cell_features %>% names())
# columns in z_20XX not in z_lf:
z_20XX %>% dplyr::anti_join(z_lf)
# columns in z_lf not in z_20XX
z_lf %>% dplyr::anti_join(z_20XX)

cell_features <- dplyr::bind_rows(
    cell_features,
    lf_rem_cell_features)

# check that cells for all the compounds were identified
cell_features %>%
    dplyr::filter(!(Compound %in% c("PC", "NC"))) %>%
    dplyr::anti_join(top_hits, by="Compound") %>%
    dplyr::count(Compound)

top_hits %>%
    dplyr::anti_join(
        cell_features %>% dplyr::count(Compound),
        by = "Compound")


# add overall infectivity score bin column
cell_features <- cell_features %>%
    dplyr::bind_cols(
        cell_features %>%
        dplyr::mutate(cell_index = dplyr::row_number()) %>%
        dplyr::arrange(infectivity_score) %>%
        dplyr::mutate(
            infectivity_score_bin =
                cut(x = c(1:dplyr::n()), breaks = 10, labels = c(1:10))) %>%
        dplyr::arrange(cell_index) %>%
        dplyr::select(infectivity_score_bin))

## just remove some PC cells?
#cell_features_pc <- cell_features %>% dplyr::filter(Compound == "PC")
#n_pc_to_keep <- nrow(cell_features_pc) - (nrow(cell_features) - 2000000)
#cell_features <- dplyr::bind_rows(
#    cell_features %>% dplyr::filter(Compound != "PC"),
#    cell_features %>%
#        dplyr::filter(Compound == "PC") %>%
#        dplyr::sample_n(n_pc_to_keep))



treatment_cells <- cell_features %>%
    dplyr::filter(!(Compound %in% c("NC", "PC"))) %>%
    dplyr::mutate(
        rand = runif(dplyr::n()),
        keep = dplyr::case_when(
            infectivity_score_bin == 1 ~ rand < .5,
            infectivity_score_bin == 2 ~ rand < .6,
            infectivity_score_bin == 3 ~ rand < .7,
#            infectivity_score_bin == 4 ~ rand < .4,
#            infectivity_score_bin == 5 ~ rand < .5,
#            infectivity_score_bin == 6 ~ rand < .6,
#            infectivity_score_bin == 7 ~ rand < .7,
            TRUE ~ TRUE)) %>%
    dplyr::filter(keep) %>%
    dplyr::select(-keep)

nc_cells <- cell_features %>%
    dplyr::filter(Compound == "NC") %>%
    dplyr::sample_n(50000)

pc_cells <- cell_features %>%
    dplyr::filter(Compound == "PC") %>%
    dplyr::sample_n(50000)

z <- dplyr::bind_rows(
        treatment_cells,
        nc_cells,
        pc_cells)

z %>% dplyr::count(Compound)


z %>% arrow::write_parquet("product/top_hit_cells_plate_scaled_200522a.parquet")


## compute embedding...




clusters <- arrow::read_parquet(
    "~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/top_hits_umap2_15_0.0/hdbscan_clustering_min100.parquet")

cell_features <- dplyr::bind_cols(
   cell_features,
   clusters)


# gather full lf_rem cell features to re-embed
lf_rem_full_cell_features <- dplyr::bind_rows(
    arrow::read_parquet(
        file = paste0("product/covid19cq1_SARS_1999B_Cell_MasterDataTable.parquet")),
    arrow::read_parquet(
        file = paste0("product/covid19cq1_SARS_2020A_Cell_MasterDataTable.parquet")) %>%
        dplyr::mutate(
            Condition = dplyr::case_when(
                Compound == "Positive Control" ~ "PC",
                Compound == "Negative Control" ~ "NC",
                TRUE ~ "Treatment"))) %>%
#    add_infectivity_score_CQ1() %>%
    dplyr::group_by(plate_id) %>%
    dplyr::mutate_at(cell_feature_columns$feature, ~ scale(.)[,1]) %>%
    dplyr::ungroup() %>%
    dplyr::filter(Condition != "BLANK")

lf_rem_full_cell_features %>%
    arrow::write_parquet("product/lf_rem_plate_scaled_1999B_2020A_Cell_MasterDataTable.parquet")










#######################################################################

cell_features <- plate_ids %>%
    dplyr::filter(
        schema == "covid19cq1",
        plate_id %in% c("1999B", "2020A")) %>%
    plyr::adply(1, function(df) {
        plate_id <- df$plate_id[1]
        cat("Reading cells from plate ", plate_id, "\n", sep="")
        cell_features <- arrow::read_parquet(
            file=paste0("product/covid19cq1_SARS_", plate_id, "_Cell_MasterDataTable.parquet"))
    })

cell_features <- cell_features %>%
    dplyr::mutate(infectivity_score = -6.758855e-01 +
        Cells_Intensity_IntegratedIntensityEdge_Virus*1.487025e-01 +
        Cells_Intensity_MeanIntensityEdge_Virus*-3.840196e+01 +
        Cells_Intensity_MaxIntensityEdge_Virus*4.270269e+01 +
        Cells_Intensity_MaxIntensity_Virus*4.254849e+01) %>%
    dplyr::arrange(infectivity_score) %>%
    dplyr::mutate(
        infectivity_score_bin = cut(x=c(1:dplyr::n()), breaks=10, labels=c(1:10)))

cell_features <- cell_features %>%
    dplyr::mutate(
        rand = runif(dplyr::n()),
        keep = dplyr::case_when(
            infectivity_score_bin == 1 ~ rand < .1,
            infectivity_score_bin == 2 ~ rand < .2,
            infectivity_score_bin == 3 ~ rand < .3,
            infectivity_score_bin == 4 ~ rand < .4,
            infectivity_score_bin == 5 ~ rand < .5,
            infectivity_score_bin == 6 ~ rand < .6,
            TRUE ~ TRUE))

cell_features <- cell_features %>%
    dplyr::filter(keep) %>%
    dplyr::select(-rand, -keep)

cell_features %>%
    arrow::write_parquet("product/covid19cq1_SARS_1999B_SARS_2020A_2M_Cell_MasterDataTable.parquet")
