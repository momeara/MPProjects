
library(plyr)
library(tidyverse)
library(ggplot2)
library(readxl)
library(MPStats)

nucleoli_features <- arrow::read_parquet(
    paste0("product/SARS_1999B_Nucleoli_MasterDataTable.parquet"))


# count per cell
data <- nucleoli_features %>%
    dplyr::count(
       Condition,
       Remdesivir_Concentration,
       Lactoferrin_Concentration,
       rem_label,
       lf_label,
       Nucleoli_Parent_Nuclei)

p <- ggplot2::ggplot(data=data) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position ="bottom") + 
    ggplot2::geom_density(
        data = data %>% dplyr::filter(Condition == "NC"),
        mapping=ggplot2::aes(
            x=log(n)),
            color="red", size=1.5) +
    ggplot2::geom_density(
        data = data %>% dplyr::filter(Condition == "PC"),
        mapping=ggplot2::aes(
            x=log10(n)),
            color="blue", size=1.5) +
    ggplot2::geom_density(
        data = data %>% dplyr::filter(Remdesivir_Concentration > 0 | Lactoferrin_Concentration > 0),
        mapping=ggplot2::aes(
            x=log(n),
            group=lf_label,
            color=Lactoferrin_Concentration)) +
    ggplot2::facet_wrap(~rem_label) +
    ggplot2::scale_x_continuous(
        "number of nucleoli per cell",
        breaks = log10(c(1, 10, 100, 1000, 10000, 100000)),
        labels = c("1", "10", "100", "1k", "10k", "100k")) +
    ggplot2::scale_color_continuous("Lactoferrin Concentration (µg/mL)") +
    ggplot2::ggtitle(
        label="Nucleoli per cell",
        subtitle="Plate 1999B  (Red: NC, Blue: PC)")

ggplot2::ggsave(
    "product/figures/nucleoli_analysis/plate_1999B/nucleoli_per_Cell_200506.pdf",
    width=6, height=6)



