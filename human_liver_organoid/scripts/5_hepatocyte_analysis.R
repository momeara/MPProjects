

library(plyr)
library(dplyr)
library(monocle3)
library(MPStats)
library(tidymodels)
suppressPackageStartupMessages(library(DESeq2))

############
# Off-Chip #
############
load("intermediate_data/cds_CZ.Rdata")

cds_CZ <- cds_CZ %>%
    MPStats::compute_qc_covariates()

cds_CZ_hepatocyte <- cds_CZ[,
    monocle3::clusters(cds_CZ, reduction_method = "UMAP") %in% c(1, 3, 4, 6, 7)]

cds_CZ_hepatocyte <- cds_CZ_hepatocyte %>%
    monocle3::preprocess_cds(num_dim = 100)

cds_CZ_hepatocyte <- cds_CZ_hepatocyte %>%
    monocle3::reduce_dimension(
        preprocess_method = "PCA",
        verbose = TRUE)

save(cds_CZ_hepatocyte, file = "intermediate_data/cds_CZ_hepatocyte.Rdata")


###########
# On Chip #
###########

load("intermediate_data/cds_CZ_x.Rdata")

# by inspecting the distribution of marker genes across the UMAP clusters
# decide that the clusters 1, 4, 7, 9 and 11 are likely enriched for hepatocyte cells
cds_CZ_x_hepatocyte <- cds_CZ_x[,
    which(
        monocle3::clusters(cds_CZ_x, reduction_method = "UMAP") %in% c(1, 4, 7, 9, 11),
        arr.ind = TRUE)]

cds_CZ_x_hepatocyte <- cds_CZ_x_hepatocyte %>%
    MPStats::compute_qc_covariates()

cds_CZ_x_hepatocyte_filtered <- cds_CZ_x_hepatocyte[,
    SummarizedExperiment::colData(cds_CZ_x_hepatocyte)[["count_depth"]] >= 10000,
    SummarizedExperiment::colData(cds_CZ_x_hepatocyte)[["mt_fraction"]] < .3]


save(cds_CZ_x_hepatocyte, file = "intermediate_data/cds_CZ_x_hepatocyte.Rdata")

col_data <- cds_CZ_x_hepatocyte_filtered %>%
    SummarizedExperiment::colData() %>%
    data.frame()

col_data %>% dplyr::count(sample)

# randomly split cells into 6 technical replicas
col_data <- col_data %>%
    dplyr::mutate(
        group_id = rep(1:6, length.out = nrow(col_data)) %>% sample())

dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts(cds_CZ_x_hepatocyte_filtered),
    colData = col_data,
    design = ~group_id + sample)


# this takes maybe 20 minutes
dds_default <- dds %>% DESeq2::DESeq()

# these are the possible contrasts
dds_default %>% DESeq2::resultsNames()
dds_default %>% DESeq2::results()


dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts(cds_CZ_x_hepatocyte_filtered),
    colData = col_data,
    design = ~sample)


# keep genes with at least 10 reads
# reduce from 36348 to 19807 genes
dds <- dds[rowSums(counts(dds)) >= 10, ]

dds <- scran::computeSumFactors(dds)

# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# Recommendations for single-cell analysis
#   fitType:
#     - defatult 'parametric' type => warning that dispersion trend was not well captured
#       and that 'local' regression was substituted
#     - recommendation: use 'glmGamPoi' for faster dispersion and parameter estimation routines
#       for single-cell data (Ahlmann-Eltze and Huber 2020)
#   test:
#     - default: test="Wald"
#     - recommentation: test="LRT" aka Likelihood ratio test to analyze all levels of a factor at once
#     - likelihood ratio test requires a 'reduced' design   
dds <- dds %>%
    DESeq2::DESeq(
        fitType = "glmGamPoi",
        test = "LRT",
        useT = TRUE,
        minmu = 1e-6,
        minReplicatesForReplace = Inf,
        reduced = ~ 1)




# https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html

# use the apeglm to shrink estimates
# Zhu, A., Ibrahim, J.G., Love, M.I. (2018)
# Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895

# https://support.bioconductor.org/p/98833/
# So my current recommendation would be to use the p-values from
# un-shrunken LFC and then use the shrunken LFC for visualization or
# ranking of genes. This is the table you get with default DESeq =>
# results => lfcShrink.

res <- dds %>% DESeq2::results()

res <- dds %>%
    DESeq2::results() %>%
    data.frame() %>%
    tibble::rownames_to_column(var = "gene") %>%
    as_tibble() %>%
    filter(padj < .05)

