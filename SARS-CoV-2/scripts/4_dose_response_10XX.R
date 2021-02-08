
library(plyr)
library(tidyverse)
library(ggplot2)

load("intermediate_data/image_scores.Rdata")


well_scores <- image_scores %>%
    dplyr::select(
        parent_plate_barcode,
        Compound,
        dose_nM,       
        Image_Count_Cells,
        Image_Classify_Positive_PctObjectsPerBin) %>%
    dplyr::group_by(parent_plate_barcode, dose_nM) %>%
        dplyr::mutate(
            cell_count_baseline = mean(Image_Count_Cells),
            prob_pose_baseline = mean(Image_Classify_Positive_PctObjectsPerBin[which(Compound == "Negative Control")])) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(parent_plate_barcode, dose_nM, Compound) %>%
        dplyr::mutate(
            normed_cell_count = Image_Count_Cells / cell_count_baseline,
            normed_prob_pos = Image_Classify_Positive_PctObjectsPerBin / prob_pose_baseline) %>%
        dplyr::summarize(
            mean_normed_prob_pos = mean(normed_prob_pos),
            sem_normed_prob_pos = sd(normed_prob_pos)/sqrt(dplyr::n()),
            mean_normed_cell_count = mean(normed_cell_count),
            sem_normed_cell_count = sd(normed_cell_count)/sqrt(dplyr::n())) %>%
    dplyr::ungroup()

make_dose_response_plot <- function(well_scores, subtitle) {
    p <- ggplot2::ggplot(data=well_scores) +
        ggplot2::theme_bw() +
        ggplot2::geom_errorbar(
           mapping=ggplot2::aes(
             x=dose_nM,
             ymin=mean_normed_cell_count-sem_normed_cell_count,
             ymax=mean_normed_cell_count+sem_normed_cell_count),
           color="black",
           width=.1) +
        ggplot2::geom_line(
           mapping=ggplot2::aes(x=dose_nM, y=mean_normed_cell_count),
           color="black",
           size=1.5) +
        ggplot2::geom_errorbar(
           mapping=ggplot2::aes(
             x=dose_nM,
             ymin=mean_normed_prob_pos-sem_normed_prob_pos,
             ymax=mean_normed_prob_pos+sem_normed_prob_pos),
           color="red",
           width=.1) +
        ggplot2::geom_line(
           mapping=ggplot2::aes(x=dose_nM, y=mean_normed_prob_pos),
           color="red",
           size=1.5) +
        ggplot2::facet_wrap(~Compound) +
        ggplot2::ggtitle(
                     label="Dose Response",
                     subtitle=subtitle) +
        ggplot2::scale_x_log10(
           "Dose nM",
           breaks=c(50, 250, 500, 1000, 2000),
           labels=c("50nM", "250nM", "0.5uM", "1uM", "2uM")) +
        ggplot2::scale_y_continuous(
           "Fraction of per-plate negative control",
           limits=c(0, 5))          
}

p <- make_dose_response_plot(
    well_scores = well_scores %>% dplyr::filter(parent_plate_barcode==1002),
    subtitle="Plate 1002")
ggplot2::ggsave(
    "product/figures/plate_1002_scores_probpos_cell_count_200422.pdf",
    width=30, height=30)


p <- make_dose_response_plot(
    well_scores = well_scores %>% dplyr::filter(parent_plate_barcode==1003),
    subtitle="Plate 1003")
ggplot2::ggsave(
    "product/figures/plate_1003_scores_probpos_cell_count_200422.pdf",
    width=30, height=30)

    
    

