

library(plyr)
library(tidyverse)
library(ggplot2)
library(readxl)
library(MPStats)
library(arrow)

load("intermediate_data/plate_map_2021A.Rdata")
load("intermediate_data/image_scores_CX5_2021A.Rdata")

treatment_scores <- image_scores_CX5_2021A %>%
    dplyr::filter(
        !(Condition == "Treatment" &
        Hydroxychloroquine_Concentration == 0 &
        Lactoferrin_Concentration == 0)) %>%
    dplyr::mutate(
        Condition = ifelse(Condition == "NC", "Treatment", Condition)) %>%
    dplyr::group_by(
        Condition,
        Hydroxychloroquine_Concentration,
        Lactoferrin_Concentration) %>%
    dplyr::summarize(
        percent_infected_treatment_mean =
            sum(Image_Classify_Positive_NumObjectsPerBin) / sum(Image_Count_Cells)) %>%
    dplyr::ungroup()
        

# plate summary
treatment_scores %>%
    dplyr::filter(Condition != "PC") %>%
    dplyr::mutate(
        percent_infected_treatment_mean =
            signif(percent_infected_treatment_mean, 3)) %>%
    tidyr::pivot_wider(
        id_cols = Lactoferrin_Concentration,
        names_from = Hydroxychloroquine_Concentration,
        values_from = percent_infected_treatment_mean,
        names_prefix = "H") %>%
    data.frame


space_before_0_dose <- 0.0
hcq_area <- tibble::tibble(
    dose = treatment_scores %>%
        dplyr::distinct(Hydroxychloroquine_Concentration) %>%
        magrittr::extract2("Hydroxychloroquine_Concentration")) %>%
    dplyr::arrange(dose) %>%
    dplyr::mutate(
        log_dose = log10(dose),
        label = dplyr::case_when(
            dose == 0 ~ 0,
            dose > 0 ~ signif(dose, 3)) %>%
            factor() %>%
            reorder(dose),
        dose_index = dplyr::row_number(),
        # put the zero dose similar spacing to the rest
        # but add a little bit of space to sparate it
        log_dose = dplyr::case_when(
            dose_index == 1 ~ log_dose[2] - (log_dose[3] - log_dose[2]) - space_before_0_dose,
            TRUE ~ log_dose),
        low = dplyr::case_when(
            dose_index == 1 ~ log_dose[1] - (log_dose[3] - log_dose[2]) / 2,
            dose_index == 2 ~ log_dose[2] - (log_dose[3] - log_dose[2]) / 2,
            TRUE ~ log_dose - (log_dose - dplyr::lag(log_dose)) / 2),
        high = dplyr::case_when(
            dose_index == 1 ~ log_dose[1] + (log_dose[3] - log_dose[2]) / 2,
            dose_index == dplyr::n() ~
                log_dose[dplyr::n()] +
                (log_dose[dplyr::n()] - log_dose[dplyr::n() - 1]) / 2,
            TRUE ~ log_dose + (dplyr::lead(log_dose) - log_dose) / 2))

lf_area <- tibble::tibble(
    dose = treatment_scores %>%
        dplyr::distinct(Lactoferrin_Concentration) %>%
        magrittr::extract2("Lactoferrin_Concentration")) %>%
    dplyr::arrange(dose) %>%
    dplyr::mutate(
        log_dose = log10(dose),
        label = dplyr::case_when(
            dose == 0 ~ 0,
            dose > 0 ~ signif(dose*11.49, 3)) %>%
            factor() %>%
            reorder(dose),
        dose_index = dplyr::row_number(),
        # put the zero dose similar spacing to the rest
        # but add a little bit of space to sparate it
        log_dose = dplyr::case_when(
            dose_index == 1 ~ log_dose[2] - (log_dose[3] - log_dose[2]) - space_before_0_dose,
            TRUE ~ log_dose),
        low = dplyr::case_when(
            dose_index == 1 ~ log_dose[1] - (log_dose[3] - log_dose[2]) / 2,
            dose_index == 2 ~ log_dose[2] - (log_dose[3] - log_dose[2]) / 2,
            TRUE ~ log_dose - (log_dose - dplyr::lag(log_dose)) / 2),
        high = dplyr::case_when(
            dose_index == 1 ~ log_dose[1] + (log_dose[3] - log_dose[2]) / 2,
            dose_index == dplyr::n() ~
                log_dose[dplyr::n()] +
                (log_dose[dplyr::n()] - log_dose[dplyr::n() - 1]) / 2,
      TRUE ~ log_dose + (dplyr::lead(log_dose) - log_dose) / 2))

treatment_scores <- treatment_scores %>%
  dplyr::left_join(
    hcq_area %>%
      dplyr::select(
        Hydroxychloroquine_Concentration = dose,
        log_hcq_dose = log_dose,
        hcq_low = low,
        hcq_high = high),
    by = "Hydroxychloroquine_Concentration") %>%
  dplyr::left_join(
    lf_area %>%
      dplyr::select(
        Lactoferrin_Concentration = dose,
        log_lf_dose = log_dose,
        lf_low = low,
        lf_high = high),
    by = "Lactoferrin_Concentration")

make_dose_response_plot_checkerboard <- function(
  treatment_scores, subtitle) {
    plot <- ggplot2::ggplot(data = treatment_scores) +
        ggplot2::theme_bw() +
        theme(
            panel.background =
                ggplot2::element_rect(fill = "#555555", colour = "#555555"),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            legend.position = "bottom") +
        ggplot2::geom_rect(
            mapping = ggplot2::aes(
                xmin = lf_low, xmax = lf_high,
                ymin = hcq_low, ymax = hcq_high,
                fill = percent_infected_treatment_mean)) +
        ggplot2::scale_x_continuous(
            "Lactoferrin (nM)",
            breaks = lf_area$log_dose,
            labels = lf_area$label,
            expand = c(0, 0)) +
        ggplot2::scale_y_continuous(
            "Hydroxychloroquine (uM)",
            breaks = hcq_area$log_dose,
            labels = hcq_area$label,
            expand = c(0, 0)) +
        ggplot2::scale_fill_gradient(
            "% infected cells/treatment",
            low = "#00274C", # MPStats::umich_colors['blue'],
            high = "#FF00FF", # magenta MPStats::umich_colors['maize'],
      expand = c(0, 0))
}

p <- make_dose_response_plot_checkerboard(
  treatment_scores = treatment_scores,
  subtitle = "Plate 2021A")
p

ggplot2::ggsave(
    "product/figures/dose_response/plate_2021A/percent_infected_checkerboard_200523.pdf",
    width = 4, height = 4)

treatment_scores %>%
    readr::write_tsv(
        "product/figures/dose_response/plate_2021A/percent_infected_checkerboard_source_data_200523.tsv")
