
library(plyr)
library(tidyverse)
library(magrittr)
library(ggplot2)
library(drc)

source("scripts/geom_indicator.R")


load("intermediate_data/image_scores.Rdata")
meta_well_scores <- image_scores %>%
    dplyr::rename(master_plate_id = parent_plate_barcode) %>%
    dplyr::distinct(
        master_plate_id,
        Compound,
        dose_nM,
        row,
        column)

viral_intensity_well_scores <- image_scores %>%
    dplyr::select(
        master_plate_id=parent_plate_barcode,
        Compound,
        dose_nM,       
        Image_Count_Cells,
        Image_Classify_Positive_PctObjectsPerBin) %>%
    dplyr::group_by(master_plate_id, dose_nM) %>%
        dplyr::mutate(
            cell_count_baseline = mean(Image_Count_Cells),
            prob_pose_baseline = mean(Image_Classify_Positive_PctObjectsPerBin[which(Compound == "Negative Control")])) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(master_plate_id, dose_nM, Compound) %>%
        dplyr::mutate(
            normed_cell_count = Image_Count_Cells / cell_count_baseline,
            normed_prob_pos = Image_Classify_Positive_PctObjectsPerBin / prob_pose_baseline) %>%
        dplyr::summarize(
            mean_normed_prob_pos = mean(normed_prob_pos),
            sem_normed_prob_pos = sd(normed_prob_pos)/sqrt(dplyr::n()),
            mean_normed_cell_count = mean(normed_cell_count),
            sem_normed_cell_count = sd(normed_cell_count)/sqrt(dplyr::n())) %>%
    dplyr::ungroup()


Compound <- "Niclosamide"
compound_scores <- viral_intensity_well_scores %>%
    dplyr::filter(Compound == !!Compound) %>%
    dplyr::select(
        master_plate_id,
        Compound,
        dose_nM,
        viral_intensity_prob_pos = mean_normed_prob_pos,
        viral_intensity_prob_pos_sem = sem_normed_prob_pos,
        cell_count = mean_normed_cell_count,
        cell_count_sem = sem_normed_cell_count)


p <- ggplot2::ggplot(data=compound_scores) +
    ggplot2::theme_bw() +
    ggplot2::geom_errorbar(
       mapping=ggplot2::aes(
         x=dose_nM,
         ymin=cell_count-cell_count_sem,
         ymax=cell_count+cell_count_sem),
       color="black",
       width=.1) +
    ggplot2::geom_line(
       mapping=ggplot2::aes(x=dose_nM, y=cell_count),
       color="black",
       size=1.4) +    
    ggplot2::geom_errorbar(
       mapping=ggplot2::aes(
         x=dose_nM,
         ymin=viral_intensity_prob_pos-viral_intensity_prob_pos_sem,
         ymax=viral_intensity_prob_pos+viral_intensity_prob_pos_sem),
       color="red",
       width=.1) +
    ggplot2::geom_line(
       mapping=ggplot2::aes(x=dose_nM, y=viral_intensity_prob_pos),
       color="red",
       size=1.3) +
    ggplot2::scale_x_log10(
       "Drug Dose",
       breaks=c(50, 250, 500, 1000, 2000),
       labels=c("50nM", "250nM", "0.5uM", "1uM", "2uM")) +
    ggplot2::scale_y_continuous("Normalized Intensity")

ggplot2::ggsave(
    "product/figures/dose_response/Niclosamide_dose_response_200501.pdf",
    width=3, height=2)

#######################

compound_scores <- compound_scores %>%
    dplyr::mutate(
        log_dose = log10(dose_nM) - 9)

fit <- drc::drm(
    formula=viral_intensity_prob_pos ~ log_dose,
    data=compound_scores,
    fct=drc::L.4(fixed=c(NA, 0, NA, NA)))
log_dose <- seq(log10(50)-9, log10(2000)-9, length.out=200)
pred_value <- predict(fit, expand.grid(log_dose, 1))
fit_data <- data.frame(log_dose, pred_value) %>%
  dplyr::mutate(
    slope=fit$coefficients[1],
    bottom=0,
    top=fit$coefficients[2],
    ic50=fit$coefficients[3],
    dose_nM = 10^(log_dose + 9))

p <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_errorbar(
       data=compound_scores,
       mapping=ggplot2::aes(
         x=dose_nM,
         ymin=cell_count-cell_count_sem,
         ymax=cell_count+cell_count_sem),
       color="black",
       width=.1) +
    ggplot2::geom_line(
       data=compound_scores,                 
       mapping=ggplot2::aes(x=dose_nM, y=cell_count),
       color="black",
       size=1.4) +    
    ggplot2::geom_errorbar(
       data=compound_scores,                 
       mapping=ggplot2::aes(
         x=dose_nM,
         ymin=viral_intensity_prob_pos-viral_intensity_prob_pos_sem,
         ymax=viral_intensity_prob_pos+viral_intensity_prob_pos_sem),
       color="red",
       width=.1) +
    ggplot2::geom_point(
       data=compound_scores,                 
       mapping=ggplot2::aes(
         x=dose_nM,
         y=viral_intensity_prob_pos),
       color="red") +
    ggplot2::geom_line(
       data=fit_data,
       mapping=ggplot2::aes(x=dose_nM, y=pred_value),
       color="red",
       size=1.3) +
    MPStats::geom_indicator(
       data=data.frame(ic50=paste0("IC50: ", signif(10^(fit_data$ic50[1] + 9), 3), " nM")),
       mapping=ggplot2::aes(indicator=ic50),
       group=1) +
    ggplot2::scale_x_log10(
       "Drug Dose",
       breaks=c(50, 250, 500, 1000, 2000),
       labels=c("50nM", "250nM", "0.5uM", "1uM", "2uM")) +
    ggplot2::scale_y_continuous("Normalized Intensity")

ggplot2::ggsave(
    "product/figures/dose_response/Niclosamide_dose_response_fit_200501.pdf",
    width=3, height=2)

