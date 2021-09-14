

library(plyr)
library(dplyr)
library(monocle3)
library(MPStats)
library(Seurat)
library(biomaRt)


hsapiens_genes <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol", "description"),
    mart = biomaRt::useEnsembl(
        biomart = "genes",
        dataset = "hsapiens_gene_ensembl"))

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


# get warning and errors about the names of the reducedDim objects (e.g. PCA, UMAP) names
# so remote them for now
reducedDims(cds) <- NULL
# get error about missing logcounts, resolve by setting data = NULL
# https://github.com/satijalab/seurat/issues/3746#issuecomment-731419868
cds <- Seurat::as.Seurat(cds, data = NULL)

Seurat::Idents(cds) <- cds@meta.data$sample

tests <- c("DESeq2", "negbinom", "wilcox")

tests <- c("wilcox")

for (test.use in tests) {
    cat("Computing analysis for ", test.use, "\n", sep = "")

    output_path <- paste0("product/figures/CZ_CZ_2/differential_expression_sctransform_v2/reads10000_", test.use)
    if (!dir.exists(output_path)) {
        dir.create(output_path)
    }

    de_on_chip_vs_off_chip <- cds %>%
        Seurat::SCTransform(
            assay = "originalexp",
            vst.flavor = "v2") %>%
        Seurat::FindMarkers(
            ident.1 = "on_chip",
            ident.2 = "off_chip",
            test.use = test.use) %>%
        tibble::rownames_to_column(var = "gene_id")
    de_on_chip_vs_off_chip %>% readr::write_tsv(
        paste0(output_path, "/de_on_chip_vs_off_chip.tsv"))
}



output_path <- "product/figures/CZ_CZ_2/differential_expression_sctransform_v2/reads10000_DESeq2"
de_data_DESeq2 <- readr::read_tsv(
        paste0(output_path, "/de_on_chip_vs_off_chip.tsv")) %>%
    dplyr::mutate(ident_1 = "On Chip", ident_2 = "Off Chip") %>%
    dplyr::mutate(test = "DESeq2")


de_data <- dplyr::bind_rows(
    de_data_DESeq2) %>%
    dplyr::mutate(ident = paste0(ident_1, " vs ", ident_2)) %>%
    dplyr::group_by(ident) %>%
    dplyr::mutate(avg_log2FC_normed = avg_log2FC - mean(avg_log2FC)) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(
        hsapiens_genes,
        by = c("gene_id" = "ensembl_gene_id")) %>%
    dplyr::mutate(
        n_log_p_val_adj = ifelse(
            test = p_val_adj == 0,
            yes = 310 + runif(dplyr::n(), -2, 2),
            no = -log10(p_val_adj)))

save(de_data, file = "intermediate_data/differential_expression_sctransform_v2_on_chip_off_chip/de_data.Rdata")


de_data %>%
    dplyr::filter(!(hgnc_symbol %>% stringr::str_detect("^MT-"))) %>%
    dplyr::filter(p_val_adj <= 1e-20) %>%
    dplyr::arrange(test, ident_1, ident_2, p_val_adj, avg_log2FC) %>%
    readr::write_tsv(
        file = "product/figures/CZ_CZ_2/differential_expression_sctransform_v2/DESeq2_top_1e-20_20210711.tsv")

de_data %>%
    dplyr::filter(!(hgnc_symbol %>% stringr::str_detect("^MT-"))) %>%
    dplyr::filter(p_val_adj <= 1e-60) %>%
    dplyr::arrange(test, ident_1, ident_2, p_val_adj, avg_log2FC) %>%
    readr::write_tsv(
        file = "product/figures/CZ_CZ_2/differential_expression_sctransform_v2/DESeq2_top_1e-60_20210711.tsv")




###########
# On chip #
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


# get warning and errors about the names of the reducedDim objects (e.g. PCA, UMAP) names
# so remote them for now
reducedDims(cds) <- NULL
# get error about missing logcounts, resolve by setting data = NULL
# https://github.com/satijalab/seurat/issues/3746#issuecomment-731419868
cds <- Seurat::as.Seurat(cds, data = NULL) %>%
    Seurat::SCTransform(
        assay = "originalexp",
        vst.flavor = "v2")
# total step1 genes: 21343
# Total Step 1 genes: 21343
# Total overdispersed genes: 18066
# Excluding 3277 genes from Step 1 because they are not overdispersed.
# Variance stabilizing transformation of count matrix of size 21725 by 2313
# Model formula is y ~ log_umi
# Get Negative Binomial regression parameters per gene
# Using 2000 genes, 2000 cells
#   |======================================================================| 100%
# Setting estimate of  35 genes to inf as theta_mm/theta_mle < 1e-3
# # of step1 poisson genes (variance < mean): 0
# # of low mean genes (mean < 0.001): 0
# Total # of Step1 poisson genes (theta=Inf; variance < mean): 35
# Total # of poisson genes (theta=Inf; variance < mean): 3615
# Calling offset model for all 3615 poisson genes
# Found 41 outliers - those will be ignored in fitting/regularization step
# 
# Ignoring theta inf genes
# Replacing fit params for 3615 poisson genes by theta=Inf
# Second step: Get residuals using fitted parameters for 21725 genes
#   |======================================================================| 100%
# Computing corrected count matrix for 21725 genes
#   |======================================================================| 100%
# Calculating gene attributes
# Wall clock passed: Time difference of 24.29755 secs
# Determine variable features
# Place corrected count matrix in counts slot
# Centering data matrix
#   |======================================================================| 100%
# Set default assay to SCT
# 
Seurat::Idents(cds) <- cds@meta.data$sample

tests <- c("DESeq2") #, "negbinom", "wilcox")


for (test.use in tests) {
    cat("Computing analysis for ", test.use, "\n", sep = "")

    output_path <- paste0("product/figures/CZ_x/differential_expression_sctransform_v2_hepatocyte/reads10000_", test.use)
    if (!dir.exists(output_path)) {
        dir.create(output_path)
    }


    de_CvF <- cds %>%
        Seurat::FindMarkers(
            ident.1 = "Fialuridine",
            ident.2 = "Control",
            test.use = test.use) %>%
        tibble::rownames_to_column(var = "gene_id")
    de_CvF %>% readr::write_tsv(
        paste0(output_path, "/de_Fialuridine_vs_Control.tsv"))
    
    de_CvA <- cds %>%
        Seurat::SCTransform(
            assay = "originalexp",
            vst.flavor = "v2") %>%        
        Seurat::FindMarkers(
            ident.1 = "Acetaminophen",
            ident.2 = "Control",
            test.use = test.use) %>%
        tibble::rownames_to_column(var = "gene_id")
    de_CvA %>% readr::write_tsv(
        paste0(output_path, "/de_Acetaminophen_vs_Control.tsv"))
    
    de_CvI <- cds %>%
        Seurat::SCTransform(
            assay = "originalexp",
            vst.flavor = "v2") %>%        
        Seurat::FindMarkers(
            ident.1 = "Inarigivir",
            ident.2 = "Control",
            test.use = test.use) %>%
        tibble::rownames_to_column(var = "gene_id")
    de_CvI %>% readr::write_tsv(
        paste0(output_path, "/de_Inarigivir_vs_Control.tsv"))
    
    de_CvT <- cds %>%
        Seurat::SCTransform(
            assay = "originalexp",
            vst.flavor = "v2") %>%        
        Seurat::FindMarkers(
            ident.1 = "Tenofovir",
            ident.2 = "Control",
            test.use = test.use) %>%
        tibble::rownames_to_column(var = "gene_id")
    de_CvT %>% readr::write_tsv(
        paste0(output_path, "/de_Tenofovir_vs_Control.tsv"))
    
    de_CvTI <- cds %>%
        Seurat::SCTransform(
            assay = "originalexp",
            vst.flavor = "v2") %>%        
        Seurat::FindMarkers(
            ident.1 = "Tenofovir_Inarigivir",
            ident.2 = "Control",
            test.use = test.use) %>%
        tibble::rownames_to_column(var = "gene_id")
    de_CvTI %>% readr::write_tsv(
        paste0(output_path, "/de_Tenofovir_Inarigivir_vs_Control.tsv"))
   
    de_FvTI <- cds %>%
        Seurat::SCTransform(
            assay = "originalexp",
            vst.flavor = "v2") %>%        
        Seurat::FindMarkers(
            ident.1 = "Tenofovir_Inarigivir",
            ident.2 = "Fialuridine",
            test.use = test.use) %>%
        tibble::rownames_to_column(var = "gene_id")
    de_FvTI %>% readr::write_tsv(
        paste0(output_path, "/de_Fialuridine_vs_Tenofovir_Inarigivir.tsv"))
}



output_path <- "product/figures/CZ_x/differential_expression_sctransform_v2_hepatocyte/reads10000_DESeq2"
de_data_DESeq2 <- dplyr::bind_rows(
    readr::read_tsv(
        paste0(output_path, "/de_Fialuridine_vs_Control.tsv")) %>%
        dplyr::mutate(ident_2 = "Control", ident_1 = "Fialuridine"),
    readr::read_tsv(
        paste0(output_path, "/de_Acetaminophen_vs_Control.tsv")) %>%
        dplyr::mutate(ident_2 = "Control", ident_1 = "Acetaminophen"),
    readr::read_tsv(
        paste0(output_path, "/de_Inarigivir_vs_Control.tsv")) %>%
        dplyr::mutate(ident_2 = "Control", ident_1 = "Inarigivir"),
    readr::read_tsv(
        paste0(output_path, "/de_Tenofovir_vs_Control.tsv")) %>%
        dplyr::mutate(ident_2 = "Control", ident_1 = "Tenofovir"),
    readr::read_tsv(
        paste0(output_path, "/de_Tenofovir_Inarigivir_vs_Control.tsv")) %>%
        dplyr::mutate(ident_2 = "Control", ident_1 = "Tenofovir_Inarigivir"),
    readr::read_tsv(
        paste0(output_path, "/de_Tenofovir_Inarigivir_vs_Fialuridine.tsv")) %>%
        dplyr::mutate(ident_2 = "Fialuridine", ident_1 = "Tenofovir_Inarigivir")) %>%
    dplyr::mutate(test = "DESeq2")

output_path <- "product/figures/CZ_x/differential_expression/reads10000_wilcox"
de_data_wilcox <- dplyr::bind_rows(
    readr::read_tsv(
        paste0(output_path, "/de_Fialuridine_vs_Control.tsv")) %>%
        dplyr::mutate(ident_1 = "Control", ident_2 = "Fialuridine"),
    readr::read_tsv(
        paste0(output_path, "/de_Acetaminophen_vs_Control.tsv")) %>%
        dplyr::mutate(ident_1 = "Control", ident_2 = "Acetaminophen"),
    readr::read_tsv(
        paste0(output_path, "/de_Inarigivir_vs_Control.tsv")) %>%
        dplyr::mutate(ident_1 = "Control", ident_2 = "Inarigivir"),
    readr::read_tsv(
        paste0(output_path, "/de_Tenofovir_vs_Control.tsv")) %>%
        dplyr::mutate(ident_1 = "Control", ident_2 = "Tenofovir"),
    readr::read_tsv(
        paste0(output_path, "/de_Tenofovir_Inarigivir_vs_Control.tsv")) %>%
        dplyr::mutate(ident_1 = "Control", ident_2 = "Tenofovir_Inarigivir"),
    readr::read_tsv(
        paste0(output_path, "/de_Tenofovir_Inarigivir_vs_Fialuridine.tsv")) %>%
        dplyr::mutate(ident_1 = "Fialuridine", ident_2 = "Tenofovir_Inarigivir")) %>%
    dplyr::mutate(test = "wilcox")

output_path <- "product/figures/CZ_x/differential_expression/reads10000_negbinom"
de_data_negbinom <- dplyr::bind_rows(
    readr::read_tsv(
        paste0(output_path, "/de_Fialuridine_vs_Control.tsv")) %>%
        dplyr::mutate(ident_1 = "Control", ident_2 = "Fialuridine"),
    readr::read_tsv(
        paste0(output_path, "/de_Acetaminophen_vs_Control.tsv")) %>%
        dplyr::mutate(ident_1 = "Control", ident_2 = "Acetaminophen"),
    readr::read_tsv(
        paste0(output_path, "/de_Inarigivir_vs_Control.tsv")) %>%
        dplyr::mutate(ident_1 = "Control", ident_2 = "Inarigivir"),
    readr::read_tsv(
        paste0(output_path, "/de_Tenofovir_vs_Control.tsv")) %>%
        dplyr::mutate(ident_1 = "Control", ident_2 = "Tenofovir"),
    readr::read_tsv(
        paste0(output_path, "/de_Tenofovir_Inarigivir_vs_Control.tsv")) %>%
        dplyr::mutate(ident_1 = "Control", ident_2 = "Tenofovir_Inarigivir"),
    readr::read_tsv(
        paste0(output_path, "/de_Tenofovir_Inarigivir_vs_Fialuridine.tsv")) %>%
        dplyr::mutate(ident_1 = "Fialuridine", ident_2 = "Tenofovir_Inarigivir")) %>%
    dplyr::mutate(test = "negbinom")

de_data <- dplyr::bind_rows(
#    de_data_wilcox,
    de_data_DESeq2) %>%
#    de_data_negbinom) %>%
    dplyr::mutate(ident = paste0(ident_1, " vs ", ident_2)) %>%
    dplyr::group_by(ident) %>%
    dplyr::mutate(avg_log2FC_normed = avg_log2FC - mean(avg_log2FC)) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(
        hsapiens_genes,
        by = c("gene_id" = "ensembl_gene_id")) %>%
    dplyr::mutate(
        n_log_p_val_adj = ifelse(
            test = p_val_adj == 0,
            yes = 310 + runif(dplyr::n(), -2, 2),
            no = -log10(p_val_adj)))

save(de_data, file = "intermediate_data/differential_expression_sctransform_v2_hepatocyte/de_data.Rdata")



de_data %>%
    dplyr::filter(!(hgnc_symbol %>% stringr::str_detect("^MT-"))) %>%
    dplyr::filter(p_val_adj <= 1e-20) %>%
    dplyr::arrange(test, ident_1, ident_2, p_val_adj, avg_log2FC) %>%
    readr::write_tsv(
        file = "product/figures/CZ_x/differential_expression_sctransform_v2_hepatocyte/DESeq2_top_1e-20_20210716.tsv")

de_data %>%
    dplyr::filter(!(hgnc_symbol %>% stringr::str_detect("^MT-"))) %>%
    dplyr::filter(p_val_adj <= 1e-60) %>%
    dplyr::arrange(test, ident_1, ident_2, p_val_adj, avg_log2FC) %>%
    readr::write_tsv(
        file = "product/figures/CZ_x/differential_expression_sctransform_v2_hepatocyte/DESeq2_top_1e-60_20210716.tsv")

