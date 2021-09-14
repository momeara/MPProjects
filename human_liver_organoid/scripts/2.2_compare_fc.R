

library(tidyverse)
library(GGally)

load("intermediate_data/differential_expression_hepatocyte/de_data.Rdata")
compare_fc <- de_data %>%
    dplyr::select(
        test,
        ident,
        gene_id,
        hgnc_symbol,
        description,
        p_val_adj,
        avg_log2FC) %>%
    tidyr::pivot_wider(
        id_cols = c(gene_id, hgnc_symbol, description, ident),
        names_from = test,
        values_from = c(p_val_adj, avg_log2FC))


plot <- GGally::ggpairs(
    data = compare_fc %>%
        dplyr::select(
            tidyselect::starts_with("avg_log2FC")),
    lower = list(continuous = wrap("points", alpha = 0.3, size = 0.1)),
    upper = list(continuous = wrap("points", alpha = 0.3, size = 0.1)))

ggplot2::ggsave(
    filename = "product/figures/CZ_x/differential_expression_hepatocyte/compare_fc_20210609.png",
    plot = plot,
    width = 12,
    height = 12)
