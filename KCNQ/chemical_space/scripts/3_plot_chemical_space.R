

library(plyr)
library(tidyverse)
library(arrow)
library(ggrepel)

source("parameters.R")
source("scripts/plot_embedding.R")

tag <- "project_substances_decoys_APDP_20210323"
dir.create(paste0("product/figures/", tag))

project_substance_ids <- arrow::read_parquet(
    "intermediate_data/project_substances_APDP_20210323/fingerprints.parquet",
    col_select = "substance_id") %>%
    dplyr::mutate(database = "project_substances")
project_substances <- readr::read_tsv(
    "raw_data/substances_20210323.tsv") %>%
    dplyr::select(
        substance_source,
        substance_name)

project_decoy_ids <- arrow::read_parquet(
    "intermediate_data/project_decoys_APDP_20210323/fingerprints.parquet",
    col_select = "substance_id") %>%
    dplyr::mutate(database = "project_decoys")
project_decoys <- readr::read_tsv(
    "intermediate_data/project_decoys_20210323.tsv") %>%
    dplyr::transmute(
        substance_source = "Decoy",
        substance_name = substance_zinc_id)

substance_ids <- dplyr::bind_rows(
    project_substance_ids,
    project_decoy_ids)

substances <- dplyr::bind_rows(
    project_substances,
    project_decoys)


source("scripts/plot_embedding.R")
#for(dir in Sys.glob("intermediate_data/project_substances_*")) {
for (dir in c("intermediate_data/project_substances_decoys_APDP_a=1,b=1_20210323")) {
    tag <- dir %>% stringr::str_replace("intermediate_data/", "")
    cat("plotting embedding for tag ", tag, " ...\n", sep = "")
    substance_umap <- arrow::read_parquet(
        file = paste0("intermediate_data/", tag, "/umap_embedding.parquet"))
    substance_data <- project_substances %>%
        dplyr::inner_join(
            dplyr::bind_cols(
                substance_ids,
                substance_umap),
            by = c("substance_name" = "substance_id"))
    background_data <- project_decoys %>%
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
    cat("\n")
}

# check that sources are compact
source("scripts/plot_embedding.R")
tag <- "project_substances_a=1,b=1_20210323"
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
    })

