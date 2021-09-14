
library(tidyverse)
library(SummarizedExperiment)
library(monocle3)
library(MPStats)
library(GGally)

#######################
# On-chip vs Off-chip #
#######################

# off chip
load("intermediate_data/cds_CZ.Rdata")
# on chip control
load("intermediate_data/cds_CZ_2.Rdata")

cds <- monocle3::combine_cds(
    cds_list = list(
        off_chip = cds_CZ,
        on_chip = cds_CZ_2))

cds <- cds %>%
    MPStats::compute_qc_covariates()

cds <- cds[
    # genes that are observed in at least 5 cells
    SingleCellExperiment::counts(cds) %>% Matrix::rowSums() >= 5,
    # cells that have at least 1000 reads which are at most 30% mitochondrial
    (SummarizedExperiment::colData(cds)[["count_depth"]] >= 10000) &
    (SummarizedExperiment::colData(cds)[["mt_fraction"]] < .3)]

gene_summary <- SummarizedExperiment::colData(cds) %>%
    data.frame() %>%
    dplyr::distinct(sample) %>%
    dplyr::rowwise() %>%
    dplyr::do({
        sample <- .$sample
        z <- cds[, SummarizedExperiment::colData(cds)[["sample"]] == sample] %>%
            SingleCellExperiment::counts()
        result <- SingleCellExperiment::rowData(cds) %>%
            data.frame() %>%
            dplyr::mutate(
                sample = sample,
                expression_mean = rowMeans2(z),
                expression_variance = rowVars(z))
    })


gene_summary_wide <- gene_summary %>%
    dplyr::transmute(
        id, sample,
        log_expression_mean = log(expression_mean),
        log_expression_variance = log(expression_variance)) %>%
    tidyr::pivot_wider(
        id_cols = c("id"),
        names_from = "sample",
        values_from = c("log_expression_mean", "log_expression_variance"))

GGally::ggpairs(
    data = gene_summary_wide %>% dplyr::select(tidyselect::starts_with("log_expression")),
    lower = list(continuous = GGally::wrap("points", alpha = 0.3, size = 0.1)),
    upper = list(continuous = GGally::wrap("points", alpha = 0.3, size = 0.1))) +
    ggplot2::theme_bw()

ggplot2::ggsave(
    filename = "product/figures/CZ_CZ_2/differential_expression/mean_expression_by_condition_20210612.pdf",
    width = 10,
    height = 9)



###########
# On-Chip #
###########


load("intermediate_data/cds_CZ_x.Rdata")
cds <- cds_CZ_x

cds <- cds %>%
    MPStats::compute_qc_covariates()

cds <- cds[
    # genes that are observed in at least 5 cells
    SingleCellExperiment::counts(cds) %>% Matrix::rowSums() >= 5,
    # hepatocyte cluster
    (monocle3::clusters(cds, reduction_method = "UMAP") %in% c(1, 4, 7, 9, 11)) &
    # cells that have at least 1000 reads which are at most 30% mitochondrial
    (SummarizedExperiment::colData(cds)[["count_depth"]] >= 10000) &
    (SummarizedExperiment::colData(cds)[["mt_fraction"]] < .3)]



umi_matrix <- cds %>%
    SingleCellExperiment::counts() %>%
    as.matrix()


gene_summary <- SummarizedExperiment::colData(cds) %>%
    data.frame() %>%
    dplyr::distinct(sample) %>%
    dplyr::rowwise() %>%
    dplyr::do({
        sample <- .$sample
        z <- cds[, SummarizedExperiment::colData(cds)[["sample"]] == sample] %>%
            SingleCellExperiment::counts()
        result <- SingleCellExperiment::rowData(cds) %>%
            data.frame() %>%
            dplyr::mutate(
                sample = sample,
                expression_mean = rowMeans2(z),
                expression_variance = rowVars(z))
    })


plot <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(
        label = "Mean gene expression by condition") +
    ggplot2::geom_histogram(
        data = gene_summary,
        mapping = ggplot2::aes(
            x = log10(expression_mean)),
        bins = 120,
        fill = "darkblue") +
    ggplot2::scale_x_continuous(
        "Mean expression across cells") +
    ggplot2::facet_wrap(
        facets = dplyr::vars(sample)) +
    ggplot2::theme(legend.position = "bottom")

ggplot2::ggsave(
    filename = "product/figures/CZ_x/differential_expression_hepatocyte/mean_expression_by_condition_20210609.pdf",
    width = 10,
    height = 6)

load("intermediate_data/de_data.Rdata")

gene_summary <- gene_summary %>%
    dplyr::left_join(
        de_data %>% dplyr::filter(ident_1 == "Control"),
        by = c(
            "sample" = "ident_2",
            "id" = "gene_id"))


ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_abline(
        intercept = 0,
        slope = 1,
        color = "darkgray",
        size = 1) +
    ggplot2::geom_point(
        data = gene_summary %>% dplyr::filter((sample == "Control") | (test == "DESeq2")),
        mapping = ggplot2::aes(
            x = log(expression_mean),
            y = log(expression_variance) - log(expression_mean),
            color = n_log_p_val_adj),
        alpha = .6,
        size = .8) +
    ggplot2::ggtitle(
        label = "Per-gene mean vs. variance",
        subtitle = "Hepatocytes: DESeq2 differential expression vs. control") +
    ggplot2::scale_x_continuous(
        "log(Expression mean)") +
    ggplot2::scale_y_continuous(
        "log(Expression variance)") +
    ggplot2::scale_color_continuous(
        "-Log10(adjusted P-vlaue)") +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::facet_wrap(
        facets = dplyr::vars(sample))

ggplot2::ggsave(
    filename = "product/figures/CZ_x/differential_expression_hepatocyte/mean_vs_variance_DESeq2_20210609.pdf",
    width = 10,
    heigh = 7,
    useDingbats = FALSE)
