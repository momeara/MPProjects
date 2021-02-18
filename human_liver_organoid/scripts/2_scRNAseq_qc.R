

library(plyr)
library(tidyverse)
library(MPStats)
library(googlesheets4)
library(monocle3)
library(GGally)

source("parameters.R")
source("scripts/get_dataset.R")

datasets <- readr::read_tsv(
    "raw_data/datasets_20210208.tsv",
    col_types = readr::cols())

#TODO
##what are size factors?
#
#   # run when creating a cell_data_se
#   cds <- estimate_size_factors(cds)
#             method = c("mean-geometric-mean-total", "mean-geometric-mean-log-total")
#
#   estimate_sf_sparse
#      cell_total <- Matrix::colSums(counts)
#      sfs <- cell_total / exp(mean(log(cell_total)))


#   cds <- estimateDispersions(cds)    
#  when combining datasets, what is going on with 'conf' conflicts?


#' Compute important quality control covariates suggested by (Lueken and Theis 2019)
#'
#' Add count_depth, n_genes, mt_fraction values to
#' the column data for the input cell_data_set
#'
#' @param cds cell_data_set
#'@export
compute_qc_covariates <- function(cds) {
    count_depth <- cds %>%
        SingleCellExperiment::counts() %>%
        Matrix::colSums()
    SummarizedExperiment::colData(cds)[["count_depth"]] <- count_depth
    
    # compute number of non-zero entries per-column for a dgCMatrix
    # https://stackoverflow.com/a/51560622/198401
    SummarizedExperiment::colData(cds)[["n_genes"]] <- SingleCellExperiment::counts(cds)@p %>% diff()
    mt_genes <- SummarizedExperiment::rowData(cds) %>%
        data.frame %>%
        dplyr::mutate(index = dplyr::row_number()) %>%
        dplyr::filter(gene_short_name %>% stringr::str_starts("MT-"))
         
    mt_count_depth <- cds[mt_genes$index, ] %>%
        SingleCellExperiment::counts() %>%
        Matrix::colSums()

    SummarizedExperiment::colData(cds)[["mt_fraction"]] <- mt_count_depth / count_depth
    cds
}

qc_covariates <- datasets %>%
    dplyr::filter(type == "scRNAseq") %>%
    dplyr::rowwise() %>%
    dplyr::do({
        sample_name <- .$sample_name
        sample_id <- .$sample_id
        get_dataset(sample_id) %>%
            compute_qc_covariates() %>%
            SummarizedExperiment::colData() %>%
            data.frame() %>%
            dplyr::mutate(
                sample_name = sample_name,
                log_count_depth = log10(count_depth),
                log_mt_fraction = log10(mt_fraction)) %>%
            dplyr::select(
                sample_name,
                log_count_depth,
                n_genes,
                log_mt_fraction)
    }) %>%
    dplyr::ungroup()

qc_covariates <- qc_covariates %>%
    dplyr::mutate(
        sample_name = sample_name %>%
            dplyr::recode_factor(
                `Ouchi2019` = "Ouchi2019",
                `Off-chip control` = "Off-chip control",
                `On-chip control` = "On-chip control",
                `Tenofovir` = "Tenofovir",
                `Inarigivir` = "Inarigivir",
                `Tenofovir/Inarigivir` = "Tenofovir/Inarigivir",
                `Acetaminophen` = "Acetaminophen",
                `Fialuridine` = "Fialuridine",
                .ordered = TRUE))
        
ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_vline(
        xintercept = log10(10000),
        color = "grey70") +
    ggplot2::geom_hline(
        yintercept = log10(.3),
        color = "grey70") +
    ggplot2::geom_point(
        data = qc_covariates %>%
            dplyr::filter(
                log_count_depth < log10(10000) |
                log_mt_fraction >= log10(.3)),
        mapping = ggplot2::aes(
            x = log_count_depth,
            y = log_mt_fraction),
        size = .1,
        alpha = .2) +
    ggplot2::geom_point(
        data = qc_covariates %>%
            dplyr::filter(
                log_count_depth >= log10(10000),
                log_mt_fraction < log10(.3)),
        mapping = ggplot2::aes(
            x = log_count_depth,
            y = log_mt_fraction),
        color = "blue",
        size = .1,
        alpha = .2) +
    ggplot2::facet_wrap(~sample_name) +
    ggplot2::scale_x_continuous(
        "Reads per barcode",
        breaks = log10(c(300, 1000, 3000, 10000, 30000, 100000, 300000)),
        labels = c("300", "1k", "3k", "10k", "30k", "100k", "300k")) +
    ggplot2::scale_y_continuous(
        "Fraction mitochondrial genes",
        breaks = log10(c(.001, .01, 0.1, 1)),
        labels = c("0.1%", "1%", "10%", "100%"))

ggplot2::ggsave(
    "product/figures/all/qc_covariates_pairs_gated_20210214.pdf")




#' Compute pearson residuals for each gene
#' Jan Lause, Philipp Berens, and Dmitry Kobak
# Analytic Pearson residuals for normalization of single-cell RNA-seq UMI data


#
