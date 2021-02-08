library(plyr)
library(tidyverse)
library(ggplot2)
library(readxl)
library(MPStats)


cell_features <- arrow::read_parquet(
    "product/SARS_1999B_Cell_MasterDataTable.parquet")

cell_features <- cell_features %>%
    dplyr::mutate(jonny_score =  -6.758855e-01 +
        Cells_Intensity_IntegratedIntensityEdge_Virus*1.487025e-01 +
        Cells_Intensity_MeanIntensityEdge_Virus*-3.840196e+01 +
        Cells_Intensity_MaxIntensityEdge_Virus*4.270269e+01 +
        Cells_Intensity_MaxIntensity_Virus*4.254849e+01)

score_densities <- cell_features %>%
    plyr::ddply(c("Condition", "Remdesivir_Concentration", "Lactoferrin_Concentration", "row", "column"), function(df){
        dens <- density(df$jonny_score)
        tibble::tibble(
            plate_id = plate_id,
            value = dens$x,
            density = dens$y)
    })               


p <- ggplot2::ggplot(
    data=score_densities) +
    ggplot2::geom_line(
        mapping=ggplot2::aes(x=value, y=log10(density+1), group=paste0(row, column, sep=" "))) +
    ggplot2::facet_grid(Remdesivir_Concentration ~ Lactoferrin_Concentration) +
    ggplot2::scale_x_log10("Jonny Cell-infectivity score") +
    ggplot2::scale_y_continuous("log(density + 1)", limits=c(0, .25))

ggplot2::ggsave(
    "product/figures/feature_densities/jonny_score_200504.pdf",
    width=15, height=15)
