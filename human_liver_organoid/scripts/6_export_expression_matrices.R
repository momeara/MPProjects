
library(tidyverse)
library(MPStats)
library(monocle3)
library(Seurat)
library(biomaRt)

hsapiens_genes <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol", "description"),
    mart = biomaRt::useEnsembl(
        biomart = "genes",
        dataset = "hsapiens_gene_ensembl"))



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


reducedDims(cds) <- NULL
# get error about missing logcounts, resolve by setting data = NULL
# https://github.com/satijalab/seurat/issues/3746#issuecomment-731419868
cds <- Seurat::as.Seurat(cds, data = NULL) %>%
    Seurat::SCTransform(
        assay = "originalexp",
        vst.flavor = "v2")

Seurat::Idents(cds) <- cds@meta.data$sample

counts <- dplyr::bind_cols(
    data.frame(
        gene_id = cds@assays$SCT@counts@Dimnames[[1]]) %>%
        dplyr::left_join(
            hsapiens_genes %>%
            dplyr::distinct(ensembl_gene_id, .keep_all = TRUE),
            by = c("gene_id" = "ensembl_gene_id")),
    cds@assays$SCT@counts %>%
        as.matrix() %>%
        as.data.frame())


# dim y: 21725  2313


counts %>% as.data.frame() %>% readr::write_tsv("product/figures/CZ_x_hepatocyte/normalized_counts.tsv")


