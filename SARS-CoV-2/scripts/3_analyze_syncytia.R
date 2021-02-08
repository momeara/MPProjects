
library(plyr)
library(tidyverse)
library(ggplot2)
library(readxl)
library(MPStats)

syncytia_features <- arrow::read_parquet(
    paste0("product/SARS_1999B_Syncytia_MasterDataTable.parquet"))

# count nuclei children per syncytia
data <- syncytia_features %>%
    dplyr::count(
       Condition,
       Remdesivir_Concentration,
       Lactoferrin_Concentration,
       rem_label,
       lf_label,
       syncytia_Children_Nuclei_Count)

p <- ggplot2::ggplot(data=data) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position ="bottom") + 
    ggplot2::geom_density(
        data = data %>% dplyr::filter(Condition == "NC"),
        mapping=ggplot2::aes(
            x=log(n),
            y=log(..density..+1)),
            color="red", size=1.5) +
    ggplot2::geom_density(
        data = data %>% dplyr::filter(Condition == "PC"),
        mapping=ggplot2::aes(
            x=log10(n),
            y=log(..density..+1)),
            color="blue", size=1.5) +
    ggplot2::geom_density(
        data = data %>% dplyr::filter(Remdesivir_Concentration > 0 | Lactoferrin_Concentration > 0),
        mapping=ggplot2::aes(
            x=log(n),
            y=log(..density..+1),                             
            group=lf_label,
            color=Lactoferrin_Concentration)) +
    ggplot2::facet_wrap(~rem_label) +
    ggplot2::scale_y_continuous(
        "Density",
        breaks=log10(c(0,1,3, 10, 30)+1),
        labels=c("0", "1", "3", "10", "30")) +
    ggplot2::scale_x_continuous(
        "number of nuclei per syncytia",
        breaks = log10(c(1, 10, 100, 1000, 10000, 100000)),
        labels = c("1", "10", "100", "1k", "10k", "100k")) +
    ggplot2::scale_color_continuous("Lactoferrin Concentration (µg/mL)") +
    ggplot2::ggtitle(
        label="Nuclei per syncytia",
        subtitle="Plate 1999B  (Red: NC, Blue: PC)")

ggplot2::ggsave(
    "product/figures/syncytia_analysis/plate_1999B/nuclei_per_syncytia_200510.pdf",
    width=6, height=6)

############
## Zprime ##
############

syncytia_features <- readr::read_tsv("raw_data/plate_ids.tsv") %>%
    plyr::adply(1, function(df){
        cat("Gathering Syncytia features for plate ", df$plate_id[1], "\n", sep="")
        arrow::read_parquet(
            paste0("product/SARS_", df$plate_id[1], "_Syncytia_MasterDataTable.parquet"))
    })

syncytia_zprime_10XX <- syncytia_features %>%
    dplyr::filter(plate_id %>% stringr::str_detect("^10")) %>%
    dplyr::group_by(plate_id, Compound, row, column) %>%
    dplyr::summarize(mean_nuclei = sum(syncytia_Children_Nuclei_Count)) %>%
    dplyr::ungroup() %>%
    plyr::ddply(c("plate_id"), function(features){
        Zprime(
            positives = features %>%
                dplyr::filter(Compound == "Positive Control") %>%
                magrittr::extract2("mean_nuclei"),
            negatives = features %>%
                dplyr::filter(Compound == "Negative Control") %>%
                magrittr::extract2("mean_nuclei")) %>%
            data.frame(Zprime=.)
    })
syncytia_zprime_10XX %>% dplyr::summarize(mean_zprime=mean(Zprime, na.rm=TRUE))    


syncytia_zprime_20XX <- syncytia_features %>%
    dplyr::filter(plate_id %>% stringr::str_detect("^20")) %>%
    dplyr::group_by(plate_id, Compound, row, column) %>%
    dplyr::summarize(mean_nuclei = sum(syncytia_Children_Nuclei_Count)) %>%
    dplyr::ungroup() %>%
    plyr::ddply(c("plate_id"), function(features){
        Zprime(
            positives = features %>%
                dplyr::filter(Compound == "Positive Control") %>%
                magrittr::extract2("mean_nuclei"),
            negatives = features %>%
                dplyr::filter(Compound == "Negative Control") %>%
                magrittr::extract2("mean_nuclei")) %>%
            data.frame(Zprime=.)
    })
syncytia_zprime_20XX %>% dplyr::summarize(mean_zprime=mean(Zprime, na.rm=TRUE))    



