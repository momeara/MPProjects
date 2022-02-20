

library(plyr)
library(dplyr)
library(monocle3)
library(MPStats)

source("parameters.R")
marker_genes <- readr::read_tsv(parameters$marker_genes_fname)

# off chip control
load("intermediate_data/cds_CZ.Rdata")

plot <- cds_CZ %>% monocle3::plot_cells()
ggplot2::ggsave(filename = "product/figures/CZ/umap_20210208.pdf")

hepatocyte_iPSC <- marker_genes %>%
    dplyr::filter(gene_set == "hepatocyte", tissue_type == "iPSC")

hepatocyte_iPSC %>%
    dplyr::group_by(gene_set) %>%
    dplyr::do({
        data <- .
        cat("Plotting gene set '", data$gene_set[1], "' \n", sep = "")
        plot <- cds_CZ %>%
            monocle3::plot_cells(
                show_trajectory_graph = FALSE,
                genes = unique(data$gene),
                cell_size = 0.5,
                alpha = .8,
                rasterize = TRUE) +
            ggplot2::coord_fixed() +
            ggplot2::theme_bw()  +
            ggplot2::theme(legend.position = "bottom") +
            ggplot2::scale_color_viridis_c(
                "Log10 Expression",
                option = "A")
        ggplot2::ggsave(
            filename = paste0(
                "product/figures/CZ/umap_", data$gene_set[1], "_", data$tissue_type[1], "_", MPStats::date_code(), ".pdf"),
            width = 15,
            height = 14,
            plot = plot)
        data.frame()
    })


###################
# on chip control #
###################
load("intermediate_data/cds_CZ_2.Rdata")
plot <- cds_CZ_2 %>% monocle3::plot_cells()
ggplot2::ggsave(filename = "product/figures/CZ_2/umap_20210208.pdf")


marker_genes %>%
    dplyr::group_by(gene_set) %>%
    dplyr::do({
        data <- .
        cat("Plotting gene set ", data$gene_set[1], "\n", sep = "")
        plot <- cds_CZ_2 %>%
            monocle3::plot_cells(genes = unique(data$gene)) +
            ggplot2::coord_fixed() +
            ggplot2::theme_bw() +
            ggplot2::theme(legend.position = "bottom")
        ggplot2::ggsave(
            filename = paste0(
                "product/figures/CZ_2/umap_", data$gene_set[1], "_", MPStats::date_code(), ".png"),
            width = 15,
            height = 14,
            plot = plot)
        data.frame()
    })



# compare off and on chip
load("intermediate_data/cds_CZ_CZ_2.Rdata")
plot <- cds_CZ_CZ_2 %>% monocle3::plot_cells()
ggplot2::ggsave(
    filename = "product/figures/CZ_CZ_2/umap_aligned_k=200_clusters_5e-4_20210208.pdf")

plot <- cds_CZ_CZ_2 %>% monocle3::plot_cells(color_cells_by = "sample")
ggplot2::ggsave(
    filename = "product/figures/CZ_CZ_2/umap_aligned_k=200_20210208.pdf")


# all on-chip samples
load("intermediate_data/cds_CZ_x.Rdata")

plot <- cds_CZ_x %>%
    monocle3::plot_cells(color_cells_by = "sample")
ggplot2::ggsave(
    filename = "product/figures/CZ_x/umap_aligned_k=200_1e-5_20210208.pdf")

plot <- cds_CZ_x %>%
    monocle3::plot_cells()
ggplot2::ggsave(
    filename = "product/figures/CZ_x/umap_aligned_k=200_clusters_1e-5_20210208.pdf")


hepatocyte_iPSC <- marker_genes %>%
    dplyr::filter(gene_set == "hepatocyte", tissue_type == "iPSC")
        


hepatocyte_iPSC %>%
    dplyr::group_by(gene_set) %>%
    dplyr::do({
        data <- .
        cat("Plotting gene set '", data$gene_set[1], "' \n", sep = "")
        plot <- cds_CZ_x %>%
            monocle3::plot_cells(genes = unique(data$gene)) +
            ggplot2::coord_fixed() +
            ggplot2::theme_bw()  +
            ggplot2::theme(legend.position = "bottom")
        ggplot2::ggsave(
            filename = paste0(
                "product/figures/CZ_x/umap_", data$gene_set[1], "_", data$tissue_type[1], "_", MPStats::date_code(), ".png"),
            width = 15,
            height = 14,
            plot = plot)
        data.frame()
    })


# stellate markers
cds <- cds_CZ_x
marker_genes %>%
    dplyr::filter(gene_set == "stellate", tissue_type == "iPSC") %>%
    dplyr::filter(gene == "VIM") %>%
    dplyr::group_by(gene) %>%
    dplyr::do({
        data <- .
        cat("Plotting gene '", data$gene[1], "'  \n", sep = "")
        plot <- cds_CZ_x %>%
            monocle3::plot_cells(
                genes = unique(data$gene),
                show_trajectory_graph = FALSE,
                scale_to_range = TRUE) +
            ggplot2::coord_fixed() +
            ggplot2::theme_bw()  +
            ggplot2::scale_color_viridis_c(
                "Normalized Expression",
                option = "A") +# trans = scales::exp_trans()) +
            ggplot2::theme(legend.position = "bottom")
        ggplot2::ggsave(
            filename = paste0(
                "product/figures/CZ_x/umap_", data$gene_set[1], "_", data$tissue_type[1], "_", data$gene[1], "_", MPStats::date_code(), ".png"),
            width = 10,
            height = 10,
            plot = plot)
        data.frame()
    })


# extract sample and cluster ids for cells
cds_CZ_x_clusters <- data.frame(
    cluster_id = cds_CZ_x %>% monocle3::clusters(),
    barcode = cds_CZ_x %>% monocle3::clusters() %>% names()) %>%
    dplyr::mutate(
        sample_id = barcode %>%
            as.character() %>%
            stringr::str_replace("^[^_]+_", ""))
cds_CZ_x_clusters %>%
    dplyr::count(cluster_id, sample_id) %>%
    tidyr::pivot_wider(
        id_cols = "sample_id",
        names_from = "cluster_id",
        values_from = "n")
cds_CZ_x_clusters %>%
    readr::write_tsv("product/figures/CZ_x/cluster_labels_20210323.tsv")

hapatocyte_genes <- marker_genes %>%
    dplyr::filter(gene_set == "hepatocyte")


# Ouchi2019
load("intermediate_data/cds_Ouchi2019.Rdata")
plot <- cds_Ouchi2019 %>%
    monocle3::plot_cells(color_cells_by = "cluster")
ggplot2::ggsave(
    filename = "product/figures/Ouchi2019/umap_a=10,b=1_1e-4_20210208.pdf")

# Ouchi2019 + on_chip control
load("intermediate_data/cds_Ouchi2019_CZ_2.Rdata")
plot <- cds_Ouchi2019_CZ_2 %>%
    monocle3::plot_cells(
        color_cells_by = "sample")
ggplot2::ggsave(
    filename = "product/figures/Ouchi2019_CZ_2/umap_aligned_k=200_1e-5_20210208.pdf")

plot <- cds_Ouchi2019_CZ_2 %>%
    monocle3::plot_cells()
ggplot2::ggsave(
    filename = "product/figures/Ouchi2019_CZ_2/umap_aligned_k=200_clusters_1e-5_20210208.pdf")

# Ouchi2019 + all on chip samples
load("intermediate_data/cds_Ouchi2019_CZ_x.Rdata")

plot <- cds_Ouchi2019_CZ_x %>%
    monocle3::plot_cells(
        color_cells_by = "sample")
ggplot2::ggsave(
    filename = "product/figures/Ouchi2019_CZ_x/umap_aligned_k=200_1e-5_20210208.pdf")

plot <- cds_Ouchi2019_CZ_x %>%
    monocle3::plot_cells()
ggplot2::ggsave(
    filename = "product/figures/Ouchi2019_CZ_x/umap_aligned_k=200_clusters_1e-5_20210208.pdf")



z <- cds_Ouchi2019_CZ_x[, cds_Ouchi2019_CZ_x$sample  == "Ouchi2019"]
plot <- z %>%
    monocle3::plot_cells()
ggplot2::ggsave(
    filename = "product/figures/Ouchi2019_CZ_x/umap_Ouchi2019_only_aligned_k=200_clusters_1e-5_20210208.pdf")



####################################
# Plot re-embedding of hepatocytes #
####################################

# Off Chip
load("intermediate_data/cds_CZ_hepatocyte.Rdata")

plot <- cds_CZ_hepatocyte %>%
    monocle3::plot_cells(
        cell_size = 0.5,
        alpha = .8,
        rasterize = TRUE,
        show_trajectory_graph = FALSE) +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw()
ggplot2::ggsave(
    filename = paste0("product/figures/CZ_hepatocyte/umap_", MPStats::date_code(), ".pdf"),
    width = 6,
    height = 4)
    


OXPHOG_markers <- marker_genes %>%
    dplyr::filter(gene_set == "OXPHOG")

plot <- cds_CZ_hepatocyte %>%
    monocle3::plot_cells(
        genes = OXPHOG_markers$gene,
        show_trajectory_graph = FALSE,
        cell_size = .07,
        raster = TRUE) +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw() +
    ggplot2::scale_color_viridis_c(
        "Log10 Expression",
        option = "A") +
    ggplot2::theme(legend.position = "bottom")

ggplot2::ggsave(
    filename = paste0("product/figures/CZ_hepatocyte/umap_OXPHOG_all_", MPStats::date_code(), ".pdf"),
    width = 15,
    height = 15)


plot <- cds_CZ_hepatocyte %>%
    monocle3::plot_cells(
        genes = factor(
            x = c("GAPDH", "NDUFA4", "HIF1A", "MYC", "CHCHD10", "APOC3"),
            levels = c("GAPDH", "NDUFA4", "HIF1A", "MYC", "CHCHD10", "APOC3"),
            labels = c("GAPDH", "NDUFA4", "HIF1A", "MYC", "CHCHD10", "APOC3")),
        cell_size = 0.5,
        alpha = .8,
        raster = TRUE,
        show_trajectory_graph = FALSE) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle("Distribution of selected OXPHOG expression across hepatocyte cells") +
    ggplot2::coord_fixed() +
    ggplot2::scale_color_viridis_c(
        "Log10 Expression",
        option = "A") +
    ggplot2::theme(legend.position = "bottom")

ggplot2::ggsave(
    filename = paste0("product/figures/CZ_hepatocyte/umap_OXPHOG_some_", MPStats::date_code(), ".pdf"),
    width = 8,
    height = 6)


# On-chip

load("intermediate_data/cds_CZ_x_hepatocyte.Rdata")
plot <- cds_CZ_x_hepatocyte[,
    SummarizedExperiment::colData(cds_CZ_x_hepatocyte)[["sample"]] %in% c(
        "Control", "Fialuridine", "Tenofovir", "Tenofovir_Inarigivir")] %>%
    monocle3::plot_cells(
        color_cells_by = "sample",
        show_trajectory_graph = FALSE) +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw() +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(
    filename = "product/figures/CZ_x_hepatocyte/umap_treat_20210812.pdf",
    width = 6,
    height = 4)
    


OXPHOG_markers <- marker_genes %>%
    dplyr::filter(gene_set == "OXPHOG")

plot <- cds_CZ_x_hepatocyte %>%
    monocle3::plot_cells(
        genes = OXPHOG_markers$gene,
        show_trajectory_graph = FALSE,
        cell_size = .07,
        raster = TRUE) +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")

ggplot2::ggsave(
    filename = "product/figures/CZ_x_hepatocyte/umap_OXPHOG_20210611.pdf",
    width = 15,
    height = 15)


plot <- cds_CZ_x_hepatocyte %>%
    monocle3::plot_cells(
        genes = factor(
            x = c("GAPDH", "NDUFA4", "HIF1A", "MYC", "CHCHD10", "APOC3"),
            levels = c("GAPDH", "NDUFA4", "HIF1A", "MYC", "CHCHD10", "APOC3"),
            labels = c("GAPDH", "NDUFA4", "HIF1A", "MYC", "CHCHD10", "APOC3")),
        cell_size = 0.5,
        alpha = .8,
        show_trajectory_graph = FALSE) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle("Distribution of selected OXPHOG expression across hepatocyte cells") +
    ggplot2::coord_fixed() +
    ggplot2::scale_color_viridis_c(
        "Log10 Expression",
        option = "A") +
    ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(
    filename = "product/figures/CZ_x_hepatocyte/umap_OXPHOG_20210611.png",
    width = 6,
    height = 4)
