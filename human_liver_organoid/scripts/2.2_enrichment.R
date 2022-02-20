
library(tidyverse)
library(WebGestaltR)




########################
# On Chip vs. Off Chip #
########################
load("intermediate_data/differential_expression_sctransform_v2_on_chip_off_chip/de_data.Rdata")
de_data %>%
    dplyr::filter(p_val_adj == 0, avg_log2FC > 1) %>%
    #    dplyr::arrange(p_val_adj, desc(avg_log2FC)) %>%
#    dplyr::filter(avg_log2FC > 2) %>%
#    head() %>%
    dplyr::group_by(test, ident) %>%
    dplyr::do({
        data <- .
        cat("Analyzing ", data$ident[1], " with ", nrow(data), " genes\n", sep = "")
        z <- WebGestaltR::WebGestaltR(
            enrichMethod = "ORA",
            referenceGene = data$hgnc_symbol,
            referenceGeneType = "genesymbol",
            enrichDatabase = "geneontology_Biological_Process",
            fdrThr = 1,
            interestGene = data %>%
                purrr::pluck("hgnc_symbol"),
            interestGeneType = "genesymbol",
            projectName = paste0("ORA_on_chip_vs_off_chip_top20_", data$test[1], "_", data$ident[1]))
        })




######################
# On Chip Hepatocyte #
######################

load("intermediate_data/differential_expression_hepatocyte/de_data.Rdata")

de_data %>%
    dplyr::filter(p_val_adj <= 1e-30) %>%
    dplyr::group_by(test, ident) %>%
    dplyr::do({
        data <- .
        cat("Analyzing ", data$ident[1], "\n", sep = "")
        z <- WebGestaltR::WebGestaltR(
            enrichMethod = "ORA",
            referenceGene = data$hgnc_symbol,
            referenceGeneType = "genesymbol",
            enrichDatabase = "geneontology_Biological_Process",
            interestGene = data %>%
                purrr::pluck("hgnc_symbol"),
            interestGeneType = "genesymbol",
            projectName = paste0("ORA_1e-30_", data$test[1], "_", data$ident[1]))
        })

de_data %>%
    dplyr::group_by(test, ident) %>%
    dplyr::do({
        data <- .
        cat("Analyzing ",data$test[1], " ", data$ident[1], "\n", sep = "")
        z <- WebGestaltR::WebGestaltR(
            enrichMethod = "GSEA",
            enrichDatabase = "geneontology_Biological_Process",
            interestGene = data %>%
                dplyr::transmute(
                    gene = hgnc_symbol,
                    score = -log10(p_val_adj)),
            interestGeneType = "genesymbol",
            projectName = paste0("GSEA_", data$test[1], "_", data$ident[1]))
        })





umi_matrix <- cds %>%
    SingleCellExperiment::counts() %>%
    as.matrix()

vst_out <- umi_matrix %>%
    sctransform::vst(
        return_gene_attr = TRUE,
        return_cell_attr = TRUE,
        method = "glmGamPoi",
        verbosity = 2)
cds_CZ_x_filtered_vst_poisson <- vst_out
save(cds_CZ_x_filtered_vst, file = "intermediate_data/cds_CZ_x_filtered_vst_poisson.Rdata")
n_genes <- dim(cds)[1]
n_cells <- dim(cds)[2]


ln_p_g <- matrix(
    log(rowMeans2(counts(cds)))[row(cds)] -
    log(mean(colSums2(counts(cds)))),
    nrow = n_genes,
    ncol = n_cells)
ln_n_c <- matrix(
    log(rowSums2(counts(cds)))[row(cds)],
    nrow = n_genes,
    ncol = n_cells)
mu_cg <- ln_p_g + ln_n_c

theta <- 10
pearson_residual_variance <- mu_cg + mu_cg^2 / theta
pearson_residual <- (cds_CZ_x - mu_cg) / sqrt(mu_cg + mu_cg^2 / theta)


