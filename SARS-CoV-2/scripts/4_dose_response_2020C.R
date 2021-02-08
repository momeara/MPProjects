

library(plyr)
library(tidyverse)
library(ggplot2)
library(readxl)
library(MPStats)
library(arrow)


field_scores <- readr::read_csv(
    "raw_data/plate_2020C_Remdesivir_and_Lactoferrin_200523.csv") %>%
    dplyr::select(
        Well_ID,
        Condition,
        Compound,
        Remdesivir_Concentration = `Remdesivir Concentration`,
        Remdesivir_Units = `Remdesivir Units`,
        Lactoferrin_Concentration = `Lactoferrin Concentration`,
        Lactoferrin_Units,
        percent_infected = `Percent Infected`,
        Count_Nuclei,
        Count_syn_nucs)

well_scores <- field_scores %>%
    dplyr::group_by(Remdesivir_Concentration, Lactoferrin_Concentration) %>%
    dplyr::summarize(
        percent_infected_well_mean = mean(percent_infected) %>% signif(3)) %>%
    dplyr::ungroup()

# plate summary
well_scores %>%
    tidyr::pivot_wider(
        id_cols = Lactoferrin_Concentration,
        names_from = Remdesivir_Concentration,
        values_from = percent_infected_well_mean,
        names_prefix = "R") %>%
    data.frame


space_before_0_dose <- 0.0
rem_area <- tibble::tibble(
    dose = well_scores %>%
        dplyr::distinct(Remdesivir_Concentration) %>%
        magrittr::extract2("Remdesivir_Concentration")) %>%
    dplyr::arrange(dose) %>%
    dplyr::mutate(
        log_dose = log10(dose),
        label = dplyr::case_when(
            dose == 0 ~ 0,
            dose > 0 ~ signif(dose * 10^9, 3)) %>%
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
    dose = well_scores %>%
        dplyr::distinct(Lactoferrin_Concentration) %>%
        magrittr::extract2("Lactoferrin_Concentration")) %>%
    dplyr::arrange(dose) %>%
    dplyr::mutate(
        log_dose = log10(dose),
        label = dplyr::case_when(
            dose == 0 ~ 0,
            dose > 0 ~ signif(dose * 10^9, 3)) %>%
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

well_scores <- well_scores %>%
  dplyr::left_join(
    rem_area %>%
      dplyr::select(
        Remdesivir_Concentration = dose,
        log_rem_dose = log_dose,
        rem_low = low,
        rem_high = high),
    by = "Remdesivir_Concentration") %>%
  dplyr::left_join(
    lf_area %>%
      dplyr::select(
        Lactoferrin_Concentration = dose,
        log_lf_dose = log_dose,
        lf_low = low,
        lf_high = high),
    by = "Lactoferrin_Concentration")

make_dose_response_plot_checkerboard <- function(
  well_scores, subtitle) {
    plot <- ggplot2::ggplot(data = well_scores) +
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
                ymin = rem_low, ymax = rem_high,
                fill = percent_infected_well_mean)) +
        ggplot2::scale_x_continuous(
            "Lactoferrin (nM)",
            breaks = lf_area$log_dose,
            labels = lf_area$label,
            expand = c(0, 0)) +
        ggplot2::scale_y_continuous(
            "Remdesivir (nM)",
            breaks = rem_area$log_dose,
            labels = rem_area$label,
            expand = c(0, 0)) +
        ggplot2::scale_fill_gradient(
            "Percent infected cells/well",
            low = "#00274C", # MPStats::umich_colors['blue'],
            high = "#FF00FF", # magenta MPStats::umich_colors['maize'],
      expand = c(0, 0))
}

p <- make_dose_response_plot_checkerboard(
  well_scores = well_scores,
  subtitle = "Plate 2020C")
p

ggplot2::ggsave(
    "product/figures/dose_response/plate_2020C/percent_infected_checkerboard_200523.pdf",
    width = 4, height = 4)

well_scores %>%
    readr::write_tsv(
        "product/figures/dose_response/plate_2020C/percent_infected_checkerboard_source_data_200523.tsv")
