

library(plyr)
library(dplyr)
library(monocle3)
library(MPStats)
library(Seurat)
library(SummarizedExperiment)

load("intermediate_data/cds_CZ_x.Rdata")
cds <- cds_CZ_x

cds <- cds %>%
    MPStats::compute_qc_covariates()


ribo_genes <- SummarizedExperiment::rowData(cds) %>% data.frame %>%
    dplyr::mutate(index = dplyr::row_number()) %>%
    dplyr::filter(gene_short_name %>% stringr::str_detect("^RP[LS][P]?[0-9]+[AXY]?$"))
ribo_count_depth <- cds[ribo_genes$index, ] %>%
    SingleCellExperiment::counts() %>%
    Matrix::colSums()
SummarizedExperiment::colData(cds)[["ribo_fraction"]] <-
    ribo_count_depth / SummarizedExperiment::colData(cds)[["count_depth"]]


cell_metadata <- SummarizedExperiment::colData(cds) %>%
    as.data.frame() %>%
    dplyr::mutate(is_hepatocyte =
        ifelse(
            monocle3::clusters(cds, reduction_method = "UMAP") %in% c(1, 4, 7, 9, 11),
            "Hepatocyte",
            "Non-Hepatocyte"))
                
plot <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_vline(
        xintercept = log10(10000),
        color = "grey70") +
    ggplot2::geom_hline(
        yintercept = log10(.3),
        color = "grey70") +
    ggplot2::geom_point(
        data = cell_metadata %>%
            dplyr::filter(
                count_depth < 10000 |
                mt_fraction >= .3),
        mapping = ggplot2::aes(
            x = log10(count_depth),
            y = log10(mt_fraction)),
        size = .1,
        alpha = .2) +
    ggplot2::geom_point(
        data = cell_metadata %>%
            dplyr::filter(
                count_depth >= 10000,
                mt_fraction < .3),
        mapping = ggplot2::aes(
            x = log10(count_depth),
            y = log10(mt_fraction)),
        color = "blue",
        size = .1,
        alpha = .2) +
    ggplot2::facet_grid(
        rows = dplyr::vars(sample),
        cols = dplyr::vars(is_hepatocyte)) +
    ggplot2::scale_x_continuous(
        "Reads per barcode",
        breaks = log10(c(300, 1000, 3000, 10000, 30000, 100000, 300000)),
        labels = c("300", "1k", "3k", "10k", "30k", "100k", "300k")) +
    ggplot2::scale_y_continuous(
        "Fraction mitochondrial genes",
        breaks = log10(c(.001, .01, 0.1, 1)),
        labels = c("0.1%", "1%", "10%", "100%"))

ggplot2::ggsave(
    filename = "product/figures/CZ_x/hepatocyte_qc_mt_vs_count_depth_20210611.png",
    plot = plot,
    width = 6,
    height = 12)


###############
plot <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_vline(
        xintercept = log10(10000),
        color = "grey70") +
#    ggplot2::geom_hline(
#        yintercept = log10(.3),
#        color = "grey70") +
    ggplot2::geom_point(
        data = cell_metadata %>%
            dplyr::filter(
                count_depth < 10000 |
                mt_fraction >= .3),
        mapping = ggplot2::aes(
            x = log10(count_depth),
            y = log10(ribo_fraction)),
        size = .1,
        alpha = .2) +
    ggplot2::geom_point(
        data = cell_metadata %>%
            dplyr::filter(
                count_depth >= 10000,
                mt_fraction < .3),
        mapping = ggplot2::aes(
            x = log10(count_depth),
            y = log10(ribo_fraction)),
        color = "blue",
        size = .1,
        alpha = .2) +
    ggplot2::geom_hline(
        yintercept = cell_metadata %>%
            dplyr::filter(sample == "Control") %>%
            dplyr::summarize(mean_ribo_fraction = mean(ribo_fraction)) %>%
            purrr::pluck("mean_ribo_fraction") %>%
            log10(),
        color = "lightgray",
        size = 1) +
    ggplot2::geom_hline(
        data = cell_metadata %>%
            dplyr::group_by(sample) %>%
            dplyr::summarize(mean_ribo_fraction = mean(ribo_fraction)),
        mapping = ggplot2::aes(
            yintercept = log10(mean_ribo_fraction)),
        color = "green",
        size = 1.5) +
    ggplot2::facet_grid(
        rows = dplyr::vars(sample),
        cols = dplyr::vars(is_hepatocyte)) +
    ggplot2::scale_x_continuous(
        "Reads per barcode",
        breaks = log10(c(300, 1000, 3000, 10000, 30000, 100000, 300000)),
        labels = c("300", "1k", "3k", "10k", "30k", "100k", "300k")) +
    ggplot2::scale_y_continuous(
        "Fraction ribosomal protein genes",
        breaks = log10(c(.001, .01, 0.1, 1)),
        labels = c("0.1%", "1%", "10%", "100%"))

ggplot2::ggsave(
    filename = "product/figures/CZ_x/hepatocyte_qc_ribo_vs_count_depth_20210611.png",
    plot = plot,
    width = 6,
    height = 12)



