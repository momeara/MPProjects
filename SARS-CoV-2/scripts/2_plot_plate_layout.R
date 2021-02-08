library(plyr)
library(tidyverse)
library(ggplot2)
library(MPStats)





load(file="intermediate_data/plate_map_20XX.Rdata")
plate_map_20XX %>%
    plyr::d_ply("Plate_Name", function(plate_map){
        cat("Making plate map plot for plate ", plate_map$Plate_Name[1], "\n", sep="")
        plot <- ggplot2::ggplot(data=plate_map) +
            ggplot2::theme_bw() +
            MPStats::geom_indicator(
                mapping=ggplot2::aes(
                    indicator=ifelse(
                        Condition == "Treatment",
                        paste0(Compound, "\n", signif(dose_nM, 3), " nM"),
                        Compound)),
                xpos="left", ypos=.9, group=1, size=2) +
            ggplot2::facet_grid(row ~ column) +
            ggplot2::ggtitle(label=paste0(plate_map$Plate_Name[1], " Plate Map"))
        ggplot2::ggsave(
            paste0("product/figures/plate_layout/plate_", plate_map$Plate_Name[1],"_layout_200510.pdf"),
            width=25, height=5)
})


#################
## Plate 1999B ##
#################


load("intermediate_data/plate_map_1999B.Rdata")

plot <- ggplot2::ggplot(data=plate_map_1999B) +
    ggplot2::theme_bw() +
    MPStats::geom_indicator(
        mapping=ggplot2::aes(
            indicator=ifelse(Condition == "Treatment", as.character(lf_label), Condition)),
        xpos="left", ypos=.9, group=1, size=2) +
    MPStats::geom_indicator(
        mapping=ggplot2::aes(
            indicator=ifelse(Condition == "Treatment", as.character(rem_label), "")),
        xpos="left", ypos=.1, group=1, size=2) +
    facet_grid(row ~ column)
ggplot2::ggsave(
    "product/figures/plate_1999B_layout_200503.pdf",
    width=25, height=5)
