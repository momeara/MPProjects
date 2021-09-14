


library(plyr)
library(dplyr)
library(monocle3)
library(MPStats)

source("parameters.R")
marker_genes <- readr::read_tsv(parameters$marker_genes_fname)

# off chip control
load("intermediate_data/cds_CZ.Rdata")

cds_CZ_filtered <-  cds_CZ %>%
    MPStats::compute_qc_covariates()

cds_CZ_filtered <- cds_CZ[,
    SummarizedExperiment::colData(cds_CZ)[["count_depth"]] >= 10000,
    SummarizedExperiment::colData(cds_CZ)[["mt_fraction"]] < .3]



# on chip control
load("intermediate_data/cds_CZ_2.Rdata")

cds_CZ_2_filtered <-  cds_CZ_2 %>%
    MPStats::compute_qc_covariates()

cds_CZ_2_filtered <- cds_CZ_2[,
    SummarizedExperiment::colData(cds_CZ_2)[["count_depth"]] >= 10000,
    SummarizedExperiment::colData(cds_CZ_2)[["mt_fraction"]] < .3]


load("intermediate_data/cds_CZ_x.Rdata")

load("intermediate_data/cds_CZ_x_hepatocyte.Rdata")


bishr_genes <- marker_genes %>%
    dplyr::filter(gene_set == "Bishr Recommendations")


plot <- cds_CZ %>%
    monocle3::plot_cells(genes = unique(bishr_genes$gene)) +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw()  +
    ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(
    filename = paste0(
        "product/figures/bishr_genes/CZ_umap_", MPStats::date_code(), ".png"),
    width = 15,
    height = 8,
    plot = plot)
data.frame()

plot <- cds_CZ_filtered %>%
    monocle3::plot_cells(genes = unique(bishr_genes$gene)) +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw()  +
    ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(
    filename = paste0(
        "product/figures/bishr_genes/CZ_filtered_umap_", MPStats::date_code(), ".png"),
    width = 15,
    height = 8,
    plot = plot)
data.frame()


plot <- cds_CZ_2 %>%
    monocle3::plot_cells(genes = unique(bishr_genes$gene)) +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw()  +
    ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(
    filename = paste0(
        "product/figures/bishr_genes/CZ_2_umap_", MPStats::date_code(), ".png"),
    width = 15,
    height = 8,
    plot = plot)
data.frame()


plot <- cds_CZ_2_filtered %>%
    monocle3::plot_cells(genes = unique(bishr_genes$gene)) +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw()  +
    ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(
    filename = paste0(
        "product/figures/bishr_genes/CZ_2_filtered_umap_", MPStats::date_code(), ".png"),
    width = 15,
    height = 8,
    plot = plot)
data.frame()

plot <- cds_CZ_x %>%
    monocle3::plot_cells(genes = unique(bishr_genes$gene)) +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw()  +
    ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(
    filename = paste0(
        "product/figures/bishr_genes/CZ_x_umap_", MPStats::date_code(), ".pdf"),
    width = 15,
    height = 8,
    useDingbats = FALSE,
    plot = plot)
data.frame()


plot <- cds_CZ_x_hepatocyte %>%
    monocle3::plot_cells(genes = unique(bishr_genes$gene)) +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw()  +
    ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(
    filename = paste0(
        "product/figures/bishr_genes/CZ_x_hepatocyte_umap_", MPStats::date_code(), ".png"),
    width = 15,
    height = 8,
    plot = plot)
data.frame()


plot <- cds_CZ_x_hepatocyte_filtered %>%
    monocle3::plot_cells(genes = unique(bishr_genes$gene)) +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw()  +
    ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(
    filename = paste0(
        "product/figures/bishr_genes/CZ_x_hepatocyte_filtered_umap_", MPStats::date_code(), ".png"),
    width = 15,
    height = 8,
    plot = plot)
data.frame()
