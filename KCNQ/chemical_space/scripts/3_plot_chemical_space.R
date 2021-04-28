

library(plyr)
library(tidyverse)
library(arrow)
library(ggrepel)

source("parameters.R")
source("scripts/plot_embedding.R")

#tag <- "project_substances_decoys_chembl25_275k_20210323"
tag <- "project_substances_decoys_20210323"
dir.create(paste0("product/figures/", tag))

project_substance_ids <- arrow::read_parquet(
    file = "intermediate_data/project_substances_APDP_20210323/fingerprints.parquet",
    col_select = "substance_id") %>%
    dplyr::mutate(database = "project_substances")
    
project_substances <- readr::read_tsv(
    file = "raw_data/substances_20210323.tsv") %>%
    dplyr::select(
        substance_source,
        substance_name,
        substance_dock_id) %>%
    dplyr::mutate(
        substance_source = ifelse(substance_name == "SCR-2682", "", substance_source)) %>%
    dplyr::mutate(substance_dock_id = substance_dock_id %>% stringr::str_sub(-14, -1))
    

project_substances_docked_features <- readr::read_tsv(
    file = "intermediate_data/KCNQ2_7CR2_retigabine_AB_20201028,project_20210323,,20210402/docking_features.tsv") %>%
    dplyr::mutate(substance_dock_id = name %>% stringr::str_replace("[.][0-9]+$", "")) %>%
    dplyr::mutate(score_tranche = dplyr::case_when(
        total_energy < -28 ~ "Top",
        total_energy < -10 ~ "Mid",
        TRUE ~ "Bad"))

project_substances <- project_substances %>%
    dplyr::left_join(
        project_substances_docked_features,
        by = "substance_dock_id")




project_decoy_ids <- arrow::read_parquet(
    "intermediate_data/project_decoys_APDP_20210323/fingerprints.parquet",
    col_select = "substance_id") %>%
    dplyr::mutate(database = "project_decoys")
project_decoys <- readr::read_tsv(
    "intermediate_data/project_decoys_20210323.tsv") %>%
    dplyr::transmute(
        substance_source = "Decoy",
        substance_name = substance_zinc_id)


chembl25_substance_ids <- arrow::read_parquet(
    "intermediate_data/chembl25_substances_20210323/fingerprints.parquet",
    col_select = "substance_id") %>%
    dplyr::mutate(database = "project_decoys")
chembl25_substances <- arrow::read_parquet(
    "intermediate_data/chembl25_substances_20210323/fingerprints.parquet",
    col_select = "substance_id") %>%
    dplyr::transmute(
        substance_source = "Chembl25-275k",
        substance_name = substance_id)


substance_ids <- dplyr::bind_rows(
    project_substance_ids,
    project_decoy_ids)

    chembl25_substance_ids)

substances <- dplyr::bind_rows(
    project_substances,
    project_decoys)

    chembl25_substances)


source("scripts/plot_embedding.R")
for(dir in Sys.glob("intermediate_data/project_substances_decoys_APDP_a=1,b=1.8,metric=cosine_20210323")) {    
#for(dir in Sys.glob("intermediate_data/project_substances_decoys_chembl25_275k_a=1,b=1.8,n_neighbors=100,metric=cosine_20210323")) {
#for (dir in c("intermediate_data/project_substances_decoys_chembl25-50k_20210323")) {
    tag <- dir %>% stringr::str_replace("intermediate_data/", "")
    cat("plotting embedding for tag ", tag, " ...\n", sep = "")
    tryCatch({
        browser()
        substance_umap <- arrow::read_parquet(
            file = paste0("intermediate_data/", tag, "/umap_embedding.parquet"))
        substance_data <- project_substances %>%
            dplyr::inner_join(
                dplyr::bind_cols(
                    substance_ids,
                    substance_umap),
                by = c("substance_name" = "substance_id"))
        background_data <- dplyr::bind_rows(project_decoys) %>%        
#        background_data <- dplyr::bind_rows(project_decoys, chembl25_substances) %>%
            dplyr::inner_join(
                dplyr::bind_cols(
                    substance_ids,
                    substance_umap),
                by = c("substance_name" = "substance_id"))
        plot_embedding(
            substance_data = substance_data,
            label = substance_source,
            background_data = background_data,
            dataset_tag = tag,
            plot_tag = "source")
        plot_embedding(
            substance_data = substance_data,
            label = score_tranche,
            background_data = background_data,
            dataset_tag = tag,
            plot_tag = "score_tranche")
        substance_data <- substance_data %>%
            dplyr::mutate(
                to_classify = ifelse(
                    is.na(substance_source),
                    substance_name,
                    "classified"))
        plot_embedding(
            substance_data = substance_data,
            label = to_classify,
            background_data = substance_data,
            dataset_tag = tag,
            plot_tag = "to_classify")
    }, error = function(e) {
        cat("Error: ", e$message, "\n", sep = "")
    })
    cat("\n")
}

# check that sources are compact
source("scripts/plot_embedding.R")
tag <- "project_substances_decoys_chembl25-50k_20210323"
substance_data %>%
    dplyr::semi_join(
        substance_data %>%
            dplyr::count(substance_source) %>%
            dplyr::filter(n > 1),
        by = "substance_source") %>%
    dplyr::group_by(substance_source) %>%
    dplyr::do({
        substance_group <- .
        cat("Plotting substances for source ", substance_group$substance_source[1], "\n", sep = "")
        plot_embedding(
            substance_data = substance_group,
            label = substance_source,
            background_data = substance_data,
            dataset_tag = tag,
            plot_tag = paste0(substance_group$substance_source[1], "_with_background"))
        tryCatch({
            plot_embedding(
                substance_data = substance_group,
                label = substance_name,
                dataset_tag = tag,
                plot_tag = paste0(substance_group$substance_source[1], "_by_substance"))
        }, error = function(msg) {
            cat("Skipping: ", msg$message, "\n", sep = "")
        })
        data.frame()
    })

