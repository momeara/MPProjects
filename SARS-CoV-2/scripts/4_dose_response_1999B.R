library(plyr)
library(tidyverse)
library(ggplot2)
library(readxl)
library(MPStats)

load("intermediate_data/plate_map_1999B.Rdata")
load("intermediate_data/image_scores_1999B.Rdata")
meta_scores <- plate_map_1999B %>%
    dplyr::distinct(
        Condition,
        Compound,       
        Remdesivir_Concentration,
        Lactoferrin_Concentration,
        rem_label,
        lf_label,
        row,
        column) 

###################
## Image probpos ##
###################

viral_intensity_scores <- image_scores_1999B %>%
    dplyr::group_by(
         Condition,
         Remdesivir_Concentration,
         Lactoferrin_Concentration,
         row,
         column) %>%
    dplyr::summarize(
         well_cell_count = sum(Image_Count_Cells),
         well_nucleoli_count = sum(Image_Count_Nucleoli),               
         well_infected_count = sum(Image_Classify_Positive_NumObjectsPerBin)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
          plate_cell_count_mean = mean(well_cell_count),
          plate_nucleoli_count_mean = mean(well_nucleoli_count),
          plate_infected_count_mean = mean(well_infected_count)) %>%
    dplyr::group_by(
         Condition,
         Remdesivir_Concentration,
         Lactoferrin_Concentration,
         plate_cell_count_mean,
         plate_nucleoli_count_mean,
         plate_infected_count_mean) %>%
    dplyr::summarize(
         cell_count_mean = mean(well_cell_count),
         cell_count_sem = sd(well_cell_count)/sqrt(dplyr::n()),
         nucleoli_count_mean = mean(well_nucleoli_count),
         nucleoli_count_sem = sd(well_nucleoli_count)/sqrt(dplyr::n()),
         infected_count_mean = mean(well_infected_count),
         infected_count_sem = sd(well_infected_count)/sqrt(dplyr::n())) %>%
    dplyr::ungroup()

#################
## Jonny score ##
#################
source("scripts/classifiers.R")
manual_score <- cell_features %>%
    dplyr::mutate(manual_score = manual_score(cell_features)) %>%
    dplyr::group_by(
         Condition,
         Remdesivir_Concentration,
         Lactoferrin_Concentration,
         row,
         column) %>%
    dplyr::summarize(
         manual_score_well_mean = mean(manual_score)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
          manual_score_plate_mean = mean(manual_score_well_mean)) %>%
    dplyr::group_by(
         Condition,
         Remdesivir_Concentration,
         Lactoferrin_Concentration,
         manual_score_plate_mean) %>%
    dplyr::summarize(
         manual_score_condition_mean = mean(manual_score_well_mean),
         manual_score_condition_mean = sd(manual_score_well_mean)/sqrt(dplyr::n())) %>%
    dplyr::ungroup()


####################
## Combine Scores ##
####################

well_scores <- meta_scores %>%
    dplyr::left_join(
        viral_intensity_scores %>%
        dplyr::select(
            Condition,
            Remdesivir_Concentration,
            Lactoferrin_Concentration,
            plate_cell_count_mean,
            plate_nucleoli_count_mean,
            plate_infected_count_mean,
            cell_count_mean,
            cell_count_sem,
            nucleoli_count_mean,
            nucleoli_count_sem,
            infected_count_mean,
            infected_count_sem),
        by=c("Remdesivir_Concentration", "Lactoferrin_Concentration"))

well_scores %>% readr::write_tsv("product/plate_1999B_Lf_vs_Rem_well_scores.tsv")



make_dose_response_plot <- function(
    scores, subtitle, score_types=c()) {
    plot <- ggplot2::ggplot(data=scores) +
        ggplot2::theme_bw() +
        theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust=.5)) +
        ggplot2::geom_hline(yintercept=1, size=.7, color="darkgray")
    if("cell_count" %in% names(score_types)){
        cat("Adding cell count layer with color ", score_types['cell_count'], "\n", sep="")
        plot <- plot +
            ggplot2::geom_errorbar(
               mapping=ggplot2::aes(
                 x=Lactoferrin_Concentration*10^-6,
                 ymin=(cell_count_mean-cell_count_sem) / plate_cell_count_mean,
                 ymax=(cell_count_mean+cell_count_sem) / plate_cell_count_mean,
                 group=Remdesivir_Concentration),
               color=score_types['cell_count'],
               width=.1) +
            ggplot2::geom_line(
               mapping=ggplot2::aes(
                   x=Lactoferrin_Concentration*10^-6,
                   y=cell_count_mean / plate_cell_count_mean,
                   group=Remdesivir_Concentration),
               color=score_types['cell_count'],
               size=1)
    }
    if("nucleoli_count" %in% names(score_types)){
        cat("Adding nucleoli count layer with color ", score_types['nucleoli_count'], "\n", sep="")
        plot <- plot +
            ggplot2::geom_errorbar(
               mapping=ggplot2::aes(
                 x=Lactoferrin_Concentration*10^-6,
                 ymin=(nucleoli_count_mean-nucleoli_count_sem) / plate_nucleoli_count_mean,
                 ymax=(nucleoli_count_mean+nucleoli_count_sem) / plate_nucleoli_count_mean,
                 group=Remdesivir_Concentration),
               color=score_types['nucleoli_count'],
               width=.1) +
            ggplot2::geom_line(
               mapping=ggplot2::aes(
                   x=Lactoferrin_Concentration*10^-6,
                   y=nucleoli_count_mean / plate_nucleoli_count_mean,
                   group=Remdesivir_Concentration),
               color=score_types['nucleoli_count'],
               size=1)
    }
    if("viral_intensity" %in% names(score_types)){
        cat("Adding viral intensity layer with color ", score_types['viral_intensity'], "\n", sep="")        
        plot <- plot +
            ggplot2::geom_errorbar(
               mapping=ggplot2::aes(
                 x=Lactoferrin_Concentration*10^-6,
                 ymin=(infected_count_mean-infected_count_sem)/plate_infected_count_mean,
                 ymax=(infected_count_mean+infected_count_sem)/plate_infected_count_mean),
               color="red",
               width=.1) +
            ggplot2::geom_line(
               mapping=ggplot2::aes(
                   x=Lactoferrin_Concentration*10^-6,
                   y=infected_count_mean/plate_infected_count_mean),
               color="red",
               size=1)
    }
    plot <- plot +
        facet_wrap(~rem_label)+
        ggplot2::ggtitle(
           label="Dose Response (Image_Classify_Positive_NumObjectsPerBin)",
           subtitle=subtitle) +
        ggplot2::scale_x_log10(
           "Lactoferrin Dose (µg/mL)",
           breaks=c(0, 0.39, 0.78, 1.56, 3.12, 6.25, 12.5, 25)*10^-6,
           labels=c("0", "0.39", "0.78", "1.56", "3.12", "6.25", "12.5", "25")) +
        ggplot2::scale_y_continuous(
           paste0(score_types, ": ", names(score_types), collapse="   "))
}


p <- make_dose_response_plot(
    scores = scores,
    subtitle="Plate 1999B",
    score_types=c(
        cell_count="black",
        nucleoli_count="green",
        viral_intensity="red"))
ggplot2::ggsave(
    "product/figures/dose_response/plate_1999B/plate_cell_count_nucleoli_count_viral_intensity_dose_response_lf_by_rem_200505.pdf",
    width=6, height=6)



##################
## Checkerboard ##
##################



rem_area <- tibble::tibble(
  dose = well_scores %>%
    dplyr::distinct(Remdesivir_Concentration) %>%
    magrittr::extract2("Remdesivir_Concentration")) %>%
  dplyr::arrange(dose) %>%
  dplyr::mutate(log_dose = log10(dose*10^-6)) %>%
  dplyr::mutate(
    label=c("0", "3.2", "6.2", "12.4", "25", "48", "100", "200"),
    dose_index = dplyr::row_number(),
    # put the zero dose similar spacing to the rest, but add a little bit of space to sparate it
    log_dose = dplyr::case_when(
      dose_index == 1 ~ log_dose[2] - (log_dose[3] - log_dose[2]) - 0.01,
      TRUE ~ log_dose)) %>%
  dplyr::mutate(
    low = dplyr::case_when(
      dose_index == 1 ~ log_dose[1] - (log_dose[3] - log_dose[2])/2,
      dose_index == 2 ~ log_dose[2] - (log_dose[3] - log_dose[2])/2,
      TRUE ~ log_dose - (log_dose - dplyr::lag(log_dose))/2),
    high = dplyr::case_when(
      dose_index == 1 ~ log_dose[1] + (log_dose[3] - log_dose[2])/2,
      dose_index == dplyr::n() ~ log_dose[dplyr::n()] + (log_dose[dplyr::n()] - log_dose[dplyr::n() - 1])/2,
      TRUE ~ log_dose + (dplyr::lead(log_dose) - log_dose)/2))

lf_area <- tibble::tibble(
  dose = well_scores %>%
    dplyr::distinct(Lactoferrin_Concentration) %>%
    magrittr::extract2("Lactoferrin_Concentration")) %>%
  dplyr::arrange(dose) %>%
  dplyr::mutate(log_dose = log10(dose*10^-6)) %>%
  dplyr::mutate(
    label = c("0", "0.39", "0.78", "1.56", "3.12", "6.25", "12.5", "25"),
    dose_index = dplyr::row_number(),
    # put the zero dose similar spacing to the rest, but add a little bit of space to sparate it
    log_dose = dplyr::case_when(
      dose_index == 1 ~ log_dose[2] - (log_dose[3] - log_dose[2]) - 0.01,
      TRUE ~ log_dose)) %>%
  dplyr::mutate(
    low = dplyr::case_when(
      dose_index == 1 ~ log_dose[1] - (log_dose[3] - log_dose[2])/2,
      dose_index == 2 ~ log_dose[2] - (log_dose[3] - log_dose[2])/2,
      TRUE ~ log_dose - (log_dose - dplyr::lag(log_dose))/2),
    high = dplyr::case_when(
      dose_index == 1 ~ log_dose[1] + (log_dose[3] - log_dose[2])/2,
      dose_index == dplyr::n() ~ log_dose[dplyr::n()] + (log_dose[dplyr::n()] - log_dose[dplyr::n() - 1])/2,
      TRUE ~ log_dose + (dplyr::lead(log_dose) - log_dose)/2))

well_scores <- well_scores %>%
  dplyr::left_join(
    rem_area %>%
      dplyr::select(
        Remdesivir_Concentration = dose,
        log_rem_dose = log_dose,
        rem_low = low,
        rem_high = high),
    by="Remdesivir_Concentration") %>%
  dplyr::left_join(
    lf_area %>%
      dplyr::select(
        Lactoferrin_Concentration = dose,
        log_lf_dose = log_dose,
        lf_low = low,
        lf_high = high),
    by="Lactoferrin_Concentration")

make_dose_response_plot_checkerboard <- function(
  well_scores, subtitle) {
  plot <- ggplot2::ggplot(data=well_scores) +
    ggplot2::theme_bw() +
    theme(legend.position = "bottom") +
    ggplot2::geom_rect(
      mapping=ggplot2::aes(
        xmin=lf_low, xmax=lf_high,
        ymin=rem_low, ymax=rem_high,
        fill=infected_count_mean)) +
    ggplot2::ggtitle(
      label="Dose Response (Image_Classify_Positive_NumObjectsPerBin)",
      subtitle=subtitle) +
    ggplot2::scale_x_continuous(
      "Lactoferrin Dose (µg/mL)",
      breaks=lf_area$log_dose,
      labels=lf_area$label,
      expand = c(0, 0)) +
    ggplot2::scale_y_continuous(
      "Remdesivir (nM)",
        breaks=rem_area$log_dose,
        labels=rem_area$label,
        expand = c(0, 0)) +
    ggplot2::scale_fill_gradient(
      "infected cells/well",
      low="#00274C", # MPStats::umich_colors['blue'],
      high="#FF00FF", # magenta MPStats::umich_colors['maize'],
      expand = c(0, 0))
}

p <- make_dose_response_plot_checkerboard(
  well_scores = well_scores,
  subtitle="Plate 1999B")
p

ggplot2::ggsave(
  "product/figures/dose_response/plate_1999B/plate_viral_intensity_checkerboard_200503.pdf",
  width=6, height=6)
