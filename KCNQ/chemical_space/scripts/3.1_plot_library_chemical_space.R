
library(plyr)
library(tidyverse)
library(arrow)
library(ggrepel)

source("parameters.R")
source("scripts/plot_embedding.R")


###############

dataset_tag <- "library_substances_APDP"
dir.create(paste0("product/figures/", dataset_tag))


substance_ids <- dplyr::bind_rows(
    readr::read_tsv("intermediate_data/AID2239_APDP_20210501/substance_info.tsv") %>%
        dplyr::transmute(
            substance_id = PUBCHEM_SUBSTANCE_ID,
            library = "Johns Hopkins Ion Channel Center"),
    readr::read_tsv("intermediate_data/UPCMLD_1_APDP_20210501/substance_info.tsv") %>%
        dplyr::transmute(substance_id = CdId, library = "UPCMDL"),
    readr::read_tsv("intermediate_data/UPCMLD_2_APDP_20210501/substance_info.tsv") %>%
        dplyr::transmute(substance_id = CdId, library = "UPCMDL"))

substance_umap <- arrow::read_parquet(
    file = paste0("intermediate_data/", dataset_tag, "/umap_embedding.parquet"))


substance_data <- dplyr::bind_cols(
    substance_ids,
    substance_umap)

plot_tag <- "full"
plot <- ggplot2::ggplot(data = substance_data) +
    ggplot2::geom_point(
        data = substance_data,
        mapping = ggplot2::aes(
            x = UMAP_1,
            y = UMAP_2,
            color = library),
        size = .05,
        alpha = .5,
        shape = 16) +
    MPStats::geom_indicator(
        data = substance_data %>%
            dplyr::count(library) %>%
            dplyr::mutate(
                indicator = paste0("Library size: ", n)),
        mapping = ggplot2::aes(
            indicator = indicator),
        group = 1,
        xpos = "right",
        ypos = "top") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position="none") +
    ggplot2::scale_x_continuous("UMAP 1") +
    ggplot2::scale_y_continuous("UMAP 2") +
    ggplot2::coord_fixed() +
    ggplot2::facet_wrap(facets = dplyr::vars(library))
    

fname <- paste0(
    "product/figures/",
    dataset_tag, "/",
    "embedding_", plot_tag, ".png")
cat("saving plot '", fname, "' ...\n", sep = "")
ggplot2::ggsave(
    filename = fname,
    width = 12,
    height = 6,
    dpi = 1200)
    
