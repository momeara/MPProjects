
library(tidyverse)

marker_genes <- readr::read_tsv("raw_data/marker_genes_20210306.tsv")

load("intermediate_data/differential_expression_hepatocyte/de_data.Rdata")


de_data %>%
    dplyr::filter(p_val_adj < 1e-10) %>%
    dplyr::semi_join(
        marker_genes %>%
        dplyr::filter(gene_set == "OXPHOG"),
        by = c("hgnc_symbol" = "gene")) %>%
    dplyr::select(
        test, ident,
        hgnc_symbol,
        description, p_val_adj, avg_log2FC) %>%
    dplyr::arrange(
        test,
        ident,
        p_val_adj) %>%
    readr::write_tsv(
        file = "product/figures/CZ_x/differential_expression_hepatocyte/reads10000_DESeq2/OXPHOG_enrichment.tsv")
