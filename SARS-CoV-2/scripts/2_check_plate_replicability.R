library(plyr)
library(tidyverse)
library(RMySQL)
library(magrittr)
library(tictoc)
library(arrow)
library(MPStats)

#####################
## Viral Intensity ##
#####################

load("intermediate_data/image_scores_CX5_100X.Rdata")
viral_intensity_well_scores <- image_scores_CX5_100X %>%
    dplyr::select(
        plate_id,
        Compound,
        dose_nM,
        row,
        column,
        Image_Count_Cells,
        Image_Classify_Positive_PctObjectsPerBin) %>%
    dplyr::group_by(plate_id, dose_nM) %>%
        dplyr::mutate(
            cell_count_baseline = mean(Image_Count_Cells),
            prob_pose_baseline =
                mean(Image_Classify_Positive_PctObjectsPerBin[which(Compound != "Positive Control")])) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
        normed_prob_pos =
            Image_Classify_Positive_PctObjectsPerBin / prob_pose_baseline) %>%
    dplyr::group_by(plate_id, dose_nM, Compound, row, column) %>%
        dplyr::summarize(
            mean_normed_prob_pos = mean(normed_prob_pos)) %>%
    dplyr::ungroup()


viral_intensity_scores <- viral_intensity_well_scores %>%
    dplyr::filter(plate_id %>% stringr::str_detect("^1005")) %>%
    dplyr::mutate(master_plate_id = 1005) %>%
    dplyr::mutate(replica = plate_id %>% stringr::str_extract(".$")) %>%
    dplyr::select(-plate_id) %>%
    tidyr::pivot_wider(
        names_from = replica,
        names_prefix = "replica_",
        values_from = mean_normed_prob_pos)

spearman_cor <- cor(
    x = viral_intensity_scores$replica_B,
    y = viral_intensity_scores$replica_C,
    method = "spearman")


ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom") +
    MPStats::geom_indicator(
        data = data.frame(indicator = paste0("Spearman: ", signif(spearman_cor, 2))),
        mapping = ggplot2::aes(indicator = indicator),
        group = 1) +
    ggplot2::geom_point(
        data = viral_intensity_scores %>%
            dplyr::filter(!(Compound %in% c("Positive Control", "Negative Control"))),
        mapping = ggplot2::aes(
            x = replica_B,
            y = replica_C,
            color = dose_nM)) +
    ggplot2::geom_point(
        data = viral_intensity_scores %>%
            dplyr::filter(Compound == "Negative Control"),
        mapping = ggplot2::aes(
            x = replica_B,
            y = replica_C),
        color = "purple") +
    ggplot2::geom_point(
        data = viral_intensity_scores %>%
            dplyr::filter(Compound == "Positive Control"),
        mapping = ggplot2::aes(
            x = replica_B,
            y = replica_C),
        color = "green") +
    ggplot2::geom_smooth(
        data = viral_intensity_scores %>%
            dplyr::filter(!(Compound %in% c("Positive Control", "Negative Control"))),
        mapping = ggplot2::aes(
            x = replica_B,
            y = replica_C),
        method = "lm") +
    ggplot2::coord_cartesian(
        xlim = c(0, 3),
        ylim = c(0, 3)) +
    ggplot2::scale_color_continuous("Dose (nM)") +
    ggplot2::scale_x_continuous("Normalized viral intensity score replicate 1") +
    ggplot2::scale_y_continuous("Nomralized viral intensity score replicate 2") +
    ggplot2::ggtitle(
        label = "Viral intensity score replicability",
        subtitle = "qHTS plate 1005")

ggplot2::ggsave(
    file = glue::glue(
        "product/figures/batch_effects/",
        "plate_replicability_1005B_vs_1005C_viral_intensity_score_200525.pdf"))




#######################
## Infectivity Score ##
#######################


infectivity_scores <- arrow::read_parquet(
    "product/infectivity_score_well_10XX_200519.parquet") %>%
    dplyr::mutate(master_plate_id = 1005) %>%
    dplyr::filter(plate_id %>% stringr::str_detect("^1005")) %>%
    dplyr::mutate(replica = plate_id %>% stringr::str_extract(".$")) %>%
    dplyr::select(-plate_id, -infectivity_score_well_sem) %>%
    tidyr::pivot_wider(
        names_from = replica,
        names_prefix = "replica_",
        values_from = infectivity_score_well_mean)

spearman_cor <- cor(
    x = infectivity_scores$replica_B,
    y = infectivity_scores$replica_C,
    method = "spearman")


ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom") +
    MPStats::geom_indicator(
        data = data.frame(indicator = paste0("Spearman: ", signif(spearman_cor, 2))),
        mapping = ggplot2::aes(indicator = indicator),
        group = 1) +
    ggplot2::geom_point(
        data = infectivity_scores %>%
            dplyr::filter(!(Compound %in% c("Positive Control", "Negative Control"))),
        mapping = ggplot2::aes(
            x = replica_B,
            y = replica_C,
            color = dose_nM)) +
    ggplot2::geom_point(
        data = infectivity_scores %>%
            dplyr::filter(Compound == "Positive Control"),
        mapping = ggplot2::aes(
            x = replica_B,
            y = replica_C),
        color = "green") +
    ggplot2::geom_point(
        data = infectivity_scores %>%
            dplyr::filter(Compound == "Negative Control"),
        mapping = ggplot2::aes(
            x = replica_B,
            y = replica_C),
        color = "purple") +
    ggplot2::geom_smooth(
        data = infectivity_scores,
        mapping = ggplot2::aes(
            x = replica_B,
            y = replica_C),
        method = "lm") +
    ggplot2::scale_color_continuous("Dose (nM)") +
    ggplot2::scale_x_continuous("Normalized Infectivity Score Replicate 1") +
    ggplot2::scale_y_continuous("Nomralized Infectivity Score Replicate 2") +
    ggplot2::ggtitle(
        label = "Infectivity score replicability",
        subtitle = "qHTS plate 1005")

ggplot2::ggsave(
    file = glue::glue(
        "product/figures/batch_effects/",
        "plate_replicability_1005B_vs_1005C_infectivity_score_200525.pdf"))




############################
## RF Score Replicability ##
############################


load("intermediate_data/rf_scores_field_10XX.Rdata")
rf_scores <- rf_scores_field_10XX %>%
    dplyr::filter(master_plate_id == 1005) %>%
    dplyr::mutate(replica = Plate_Name %>% stringr::str_extract(".$")) %>%
    dplyr::select(
        replica,
        dose_nM,
        Compound,
        row, column,
        Image_Metadata_Field,
        infectivity_probpos_field) %>%
    dplyr::group_by(replica, dose_nM, Compound, row, column) %>%
    dplyr::summarize(infectivity_probpos_well = mean(1-infectivity_probpos_field)) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_wider(
        names_from = replica,
        names_prefix = "replica_",
        values_from = infectivity_probpos_well)

spearman_cor <- cor(
    x = rf_scores$replica_B,
    y = rf_scores$replica_C,
    method = "spearman")

ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom") +
    MPStats::geom_indicator(
        data = data.frame(indicator = paste0("Spearman: ", signif(spearman_cor, 2))),
        mapping = ggplot2::aes(indicator = indicator),
        ypos = 0.99,
        group = 1) +
    ggplot2::geom_point(
        data = rf_scores %>%
            dplyr::filter(!(Compound %in% c("PC", "NC"))),
        mapping = ggplot2::aes(
            x = replica_B,
            y = replica_C,
            color = dose_nM)) +
    ggplot2::geom_point(
        data = rf_scores %>% dplyr::filter(Compound == "PC"),
        mapping = ggplot2::aes(
            x = replica_B,
            y = replica_C),
        color = "green") +
    ggplot2::geom_point(
        data = rf_scores %>% dplyr::filter(Compound == "NC"),
        mapping = ggplot2::aes(
            x = replica_B,
            y = replica_C),
        color = "purple") +
    ggplot2::geom_smooth(
        data = rf_scores,
        mapping = ggplot2::aes(
            x = replica_B,
            y = replica_C),
        method = "lm") +
    ggplot2::scale_color_continuous("Dose (nM)") +
    ggplot2::scale_x_continuous("Normalized Infectivity Score Replicate 1") +
    ggplot2::scale_y_continuous("Nomralized Infectivity Score Replicate 2") +
    ggplot2::ggtitle(
        label = "Random forest infectivity score",
        subtitle = "qHTS plate 1005")

ggplot2::ggsave(
    file = glue::glue(
        "product/figures/batch_effects/",
        "plate_replicability_1005B_vs_1005C_rf_score_well_200508.pdf"))
