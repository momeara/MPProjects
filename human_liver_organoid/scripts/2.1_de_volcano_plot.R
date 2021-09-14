library(tidyverse)

#######################
# On-Chip vs Off-Chip #
#######################

load("intermediate_data/differential_expression_sctransform_v2_on_chip_off_chip/de_data.Rdata")

de_significant <- de_data %>%
    dplyr::group_by(ident) %>%
    dplyr::filter(
        ((rank(avg_log2FC) / dplyr::n() <= .05) |
        (rank(avg_log2FC) / dplyr::n() >= .95)) &
        (-log10(p_val_adj) > 60)) %>%
    dplyr::ungroup()

de_not_significant <- de_data %>%
    dplyr::anti_join(de_significant, by = "gene_id")


marker_genes <- readr::read_tsv("raw_data/marker_genes_20210306.tsv")

besties <- de_data %>%
    dplyr::filter((hgnc_symbol == "NDUFA4") | !stringr::str_detect(hgnc_symbol, "^NDUF")) %>%
    dplyr::filter((hgnc_symbol == "COX8A") | !stringr::str_detect(hgnc_symbol, "^COX")) %>%
    dplyr::filter((hgnc_symbol == "UQCR10") | !stringr::str_detect(hgnc_symbol, "^UQCR")) %>%
    dplyr::filter(p_val_adj < 1e-100) %>%
    dplyr::filter() %>%    
    dplyr::semi_join(
        marker_genes %>%
        dplyr::filter(gene_set == "CZ-ONCHIP-DEPICKS"),
        by = c("hgnc_symbol" = "gene"))

ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_point(
        data = de_not_significant,
        mapping = ggplot2::aes(
            x = avg_log2FC,
            y = n_log_p_val_adj),
        size = .4,
        color = "grey50",
        alpha = .3) +
    ggplot2::geom_point(
        data = de_significant,
        mapping = ggplot2::aes(
            x = avg_log2FC,
            y = n_log_p_val_adj),
        color = "darkgreen",
        size = .4,
        alpha = .5) +
    ggplot2::geom_point(
        data = besties,
        mapping = ggplot2::aes(
            x = avg_log2FC,
            y = n_log_p_val_adj),
        size = 1,
        alpha = .8) +
    ggrepel::geom_text_repel(
        data = besties,
        mapping = ggplot2::aes(
            x = avg_log2FC,
            y = n_log_p_val_adj,
            label = hgnc_symbol),
        size = 1.8,
        force = 10,
        min.segment.length = 0) +
    ggplot2::scale_x_continuous(expression("Average log[2] FC")) +
    ggplot2::scale_y_continuous(
        expression("-log[10] Corrected P-value"))


besties %>% readr::write_tsv("product/figures/CZ_CZ_2/differential_expression_sctransform_v2/volcano_CZ-ONCHIP-DEPICKS_20210711.tsv")

ggplot2::ggsave(
    filename = paste0(
        "product/figures/CZ_CZ_2/differential_expression_sctransform_v2/volcano_CZ-ONCHIP-DEPICKS_20210711.pdf"),
        width = 5,
        heigh = 4,
        useDingbats = FALSE)




###########
# On-Chip #
###########

tests <- c("DESeq2")

for (test in tests) {
    load("intermediate_data/differential_expression_sctransform_v2_hepatocyte/de_data.Rdata")
    de_data <- de_data %>%
        dplyr::filter(test == {{test}})

    de_significant <- de_data %>%
        dplyr::group_by(test, ident) %>%
        dplyr::filter(
            ((rank(avg_log2FC) / dplyr::n() <= .05) |
            (rank(avg_log2FC) / dplyr::n() >= .95)) &
            (-log10(p_val_adj) > 60)) %>%
        dplyr::ungroup()

    de_not_significant <- de_data %>%
        dplyr::anti_join(de_significant, by = "gene_id")

    besties <- de_data %>%
        dplyr::filter(!(hgnc_symbol %>% stringr::str_detect("^MT-"))) %>%
        dplyr::group_by(test, ident_1, ident_2) %>%
        dplyr::arrange(p_val_adj, avg_log2FC) %>%
        dplyr::slice(1:10) %>%
        dplyr::ungroup()

    ggplot2::ggplot() +
        ggplot2::theme_bw() +
        ggplot2::geom_point(
            data = de_not_significant,
            mapping = ggplot2::aes(
                x = avg_log2FC,
                y = n_log_p_val_adj),
            size = .4,
            color = "grey50",
            alpha = .3) +
        ggplot2::geom_point(
            data = de_significant,
            mapping = ggplot2::aes(
                x = avg_log2FC,
                y = n_log_p_val_adj),
            color = "darkgreen",
            size = .4,
            alpha = .5) +
        ggplot2::geom_point(
            data = besties,
            mapping = ggplot2::aes(
                x = avg_log2FC,
                y = n_log_p_val_adj),
            size = 1,
            alpha = .8) +
        ggrepel::geom_text_repel(
            data = besties,
            mapping = ggplot2::aes(
                x = avg_log2FC,
                y = n_log_p_val_adj,
                label = hgnc_symbol),
            size = 1.8,
            force = 10,
            min.segment.length = 0) +
        ggplot2::facet_wrap(facets = ggplot2::vars(ident)) +
        ggplot2::scale_x_continuous(expression("Average log[2] FC")) +
        ggplot2::scale_y_continuous(
            expression("-log[10] Corrected P-value"))


    ggplot2::ggsave(
        filename = paste0(
            "product/figures/CZ_x/differential_expression_sctransform_v2_hepatocyte/", test, "_volcano_20210716.pdf"),
        width = 10,
        heigh = 5,
        useDingbats = FALSE)
}
