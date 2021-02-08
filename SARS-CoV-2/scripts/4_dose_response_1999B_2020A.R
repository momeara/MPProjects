library(plyr)
library(tidyverse)
library(ggplot2)
library(readxl)
library(MPStats)
library(arrow)


load("intermediate_data/plate_map_1999B.Rdata")
load("intermediate_data/plate_map_2020A.Rdata")
load("intermediate_data/image_scores_1999B.Rdata")
image_scores_2020A <- arrow::read_parquet("product/image_scores_2020A.parquet")

meta_scores <- dplyr::bind_rows(
    plate_map_1999B %>%
        dplyr::mutate(Plate_Name = "SARS_1999B") %>%
        dplyr::distinct(
            Condition,
            Remdesivir_Concentration,
            Lactoferrin_Concentration,
            rem_label,
            lf_label),
    plate_map_2020A %>%
        dplyr::mutate(
            Plate_Name = "SARS_2020A",
            Condition = `Fluid name`) %>%
        dplyr::distinct(
            Condition,
            Remdesivir_Concentration,
            Lactoferrin_Concentration,
            rem_label,
            lf_label)) %>%
    dplyr::filter(Condition == "Treatment") %>%
    dplyr::mutate(
        Remdesivir_Concentration = signif(Remdesivir_Concentration, 3),
        Lactoferin_Concentration = signif(Lactoferrin_Concentration, 3)) %>%
    dplyr::distinct(
        Condition,
        Remdesivir_Concentration,
        Lactoferrin_Concentration,
        rem_label,
        lf_label)

###################
## Image probpos ##
###################

# average plate normalized frame scores across conditions
viral_intensity_scores_treatment <- dplyr::bind_rows(
    image_scores_1999B,
    image_scores_2020A) %>%
    dplyr::group_by(Plate_Name) %>%
    dplyr::mutate(
        cell_count_plate_mean = mean(Image_Count_Cells),
        nucleoli_count_plate_mean = mean(Image_Count_Nucleoli),
        infected_count_plate_mean = mean(Image_Classify_Positive_NumObjectsPerBin)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(
        Condition,
        Remdesivir_Concentration,
        Lactoferrin_Concentration) %>%
    dplyr::summarize(
        cell_count_treatment_mean = mean(Image_Count_Cells/cell_count_plate_mean),
        cell_count_treatment_sem = sd(Image_Count_Cells/cell_count_plate_mean)/sqrt(dplyr::n()),
        nucleoli_count_treatment_mean = mean(Image_Count_Nucleoli/nucleoli_count_plate_mean),
        nucleoli_count_treatment_sem = sd(Image_Count_Nucleoli/nucleoli_count_plate_mean)/sqrt(dplyr::n()),
        infected_count_treatment_mean = mean(Image_Classify_Positive_NumObjectsPerBin/infected_count_plate_mean),
        infected_count_treatment_sem = sd(Image_Classify_Positive_NumObjectsPerBin/infected_count_plate_mean)/sqrt(dplyr::n())) %>%
    dplyr::ungroup() %>%
    dplyr::filter(Condition == "Treatment") %>%
    dplyr::mutate(
        Remdesivir_Concentration = signif(Remdesivir_Concentration, 3),
        Lactoferin_Concentration = signif(Lactoferrin_Concentration, 3))
    



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

###############################
### Stratominer Score 200515 ##
###############################

stratominer_scores_treatment <- arrow::read_parquet(
    "intermediate_data/stratominer_well_scores_1999B_2020A_200515.parquet") %>%
    dplyr::rename(Condition=COND) %>%
    dplyr::filter(Condition == "Treatment") %>%
    dplyr::mutate(
        Remdesivir_Concentration = signif(Remdesivir_Concentration, 3),
        Lactoferin_Concentration = signif(Lactoferrin_Concentration, 3)) %>%
    dplyr::group_by(
        Condition,
        Remdesivir_Concentration,
        Lactoferrin_Concentration) %>%
    dplyr::summarize(
        stratominer_score_treatment = mean(probPOSITIVE)) %>%
    dplyr::ungroup()
    

stratominer_well_scores_1999B_2020A %>% dplyr::filter(Plate_Name == "SARS_1999B") %>%
    dplyr::mutate(
        Remdesivir_Concentration = signif(Remdesivir_Concentration, 3),
        Lactoferin_Concentration = signif(Lactoferrin_Concentration, 3)) %>%
    dplyr::group_by(Remdesivir_Concentration, Lactoferrin_Concentration) %>%
    dplyr::summarize(mean=mean(probPOSITIVE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(mean=signif(mean, 3)) %>%
    tidyr::pivot_wider(
        id_cols=Lactoferrin_Concentration,
        names_from=Remdesivir_Concentration,
        values_from=mean,names_prefix="R") %>%
    data.frame

# Lf_conc    R0 R0.0032 R0.0062 R0.0124 R0.025 R0.048  R0.1  R0.2
#   0.000 0.616   0.681   0.754   0.742  0.749  0.720 0.746 0.762
#   0.390 0.577   0.649   0.661   0.733  0.681  0.677 0.717 0.738
#   0.780 0.564   0.560   0.618   0.671  0.693  0.694 0.712 0.701
#   1.562 0.608   0.634   0.678   0.712  0.601  0.751 0.661 0.734
#   3.125 0.605   0.612   0.615   0.633  0.640  0.727 0.732 0.725
#   6.250 0.608   0.747   0.732   0.721  0.703  0.738 0.687 0.692
#  12.500 0.675   0.709   0.751   0.735  0.705  0.762 0.747 0.769
#  25.000 0.709   0.803   0.770   0.803  0.730  0.784 0.786 0.857


stratominer_well_scores_1999B_2020A %>% dplyr::filter(Plate_Name == "SARS_2020A") %>%
    dplyr::mutate(
        Remdesivir_Concentration = signif(Remdesivir_Concentration, 3),
        Lactoferin_Concentration = signif(Lactoferrin_Concentration, 3)) %>%
    dplyr::group_by(Remdesivir_Concentration, Lactoferrin_Concentration) %>%
    dplyr::summarize(mean=mean(probPOSITIVE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(mean=signif(mean, 3)) %>%
    tidyr::pivot_wider(
        id_cols=Lactoferrin_Concentration,
        names_from=Remdesivir_Concentration,
        values_from=mean,names_prefix="R") %>%
    data.frame

#   Lactoferrin_Concentration    R0 R2e.04 R4e.04 R8e.04 R0.0016 R0.0032 R0.0062 R0.0124 R0.025
# 1                         0 0.994  0.636  0.186  0.397      NA   0.657   0.739   0.820  0.746
# 2                        25 0.747  0.662     NA  0.780      NA      NA      NA      NA     NA
# 3                        50 0.636     NA  0.701     NA      NA      NA      NA      NA     NA
# 4                        75 0.770     NA  0.793  0.704   0.783      NA      NA      NA     NA
# 5                       100 0.705     NA     NA     NA   0.627   0.789      NA      NA     NA
# 6                       125 0.771     NA     NA     NA      NA   0.722   0.588      NA     NA
# 7                       150 0.765     NA     NA     NA      NA      NA   0.731   0.705     NA
# 8                       175 0.610     NA     NA     NA      NA      NA      NA   0.759  0.791
# 9                       200 0.816     NA     NA     NA      NA      NA      NA      NA  0.768

stratominer_well_scores_1999B_2020A %>%
    dplyr::mutate(
        Remdesivir_Concentration = signif(Remdesivir_Concentration, 3),
        Lactoferin_Concentration = signif(Lactoferrin_Concentration, 3)) %>%
    dplyr::group_by(Remdesivir_Concentration, Lactoferrin_Concentration) %>%
    dplyr::summarize(mean=mean(probPOSITIVE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(mean=signif(mean, 3)) %>%
    tidyr::pivot_wider(
        id_cols=Lactoferrin_Concentration,
        names_from=Remdesivir_Concentration,
        values_from=mean,names_prefix="R") %>%
    data.frame

#       Lactoferrin_Concentration    R0 R2e.04 R4e.04 R8e.04 R0.0016 R0.0032 R0.0062 R0.0124 R0.025 R0.048  R0.1  R0.2
#    1                     0.0000 0.729  0.636  0.186  0.397      NA   0.678   0.748   0.768  0.748  0.720 0.746 0.762
#    2                     0.3900 0.577     NA     NA     NA      NA   0.649   0.661   0.733  0.681  0.677 0.717 0.738
#    3                     0.7800 0.564     NA     NA     NA      NA   0.560   0.618   0.671  0.693  0.694 0.712 0.701
#    4                     1.5625 0.608     NA     NA     NA      NA   0.634   0.678   0.712  0.601  0.751 0.661 0.734
#    5                     3.1250 0.605     NA     NA     NA      NA   0.612   0.615   0.633  0.640  0.727 0.732 0.725
#    6                     6.2500 0.608     NA     NA     NA      NA   0.747   0.732   0.721  0.703  0.738 0.687 0.692
#    7                    12.5000 0.675     NA     NA     NA      NA   0.709   0.751   0.735  0.705  0.762 0.747 0.769
#    8                    25.0000 0.722  0.662     NA  0.780      NA   0.803   0.770   0.803  0.730  0.784 0.786 0.857
#    9                    50.0000 0.636     NA  0.701     NA      NA      NA      NA      NA     NA     NA    NA    NA
#    10                   75.0000 0.770     NA  0.793  0.704   0.783      NA      NA      NA     NA     NA    NA    NA
#    11                  100.0000 0.705     NA     NA     NA   0.627   0.789      NA      NA     NA     NA    NA    NA
#    12                  125.0000 0.771     NA     NA     NA      NA   0.722   0.588      NA     NA     NA    NA    NA
#    13                  150.0000 0.765     NA     NA     NA      NA      NA   0.731   0.705     NA     NA    NA    NA
#    14                  175.0000 0.610     NA     NA     NA      NA      NA      NA   0.759  0.791     NA    NA    NA
#    15                  200.0000 0.816     NA     NA     NA      NA      NA      NA      NA  0.768     NA    NA    NA


##############
## RF Score ##
##############

load("intermediate_data/rf_scores_field_1999B_2020A.Rdata")

rf_scores_treatment <- rf_scores_field_1999B_2020A %>%
    dplyr::filter(!(Compound %in% c("BLANK", "PC", "NC"))) %>%
    dplyr::mutate(
        Remdesivir_Concentration = signif(Remdesivir_Concentration, 3),
        Lactoferin_Concentration = signif(Lactoferrin_Concentration, 3)) %>%
    dplyr::group_by(Remdesivir_Concentration, Lactoferrin_Concentration) %>%
    dplyr::summarize(
        rf_score_treatment_mean = mean(1-infectivity_probpos_field),
        rf_score_treatment_sem = sd(infectivity_probpos_field)/sqrt(dplyr::n()),
        rf_score_treatment_low = pmax(0, rf_score_treatment_mean-rf_score_treatment_sem),
        rf_score_treatment_high = pmin(rf_score_treatment_mean+rf_score_treatment_sem, 1)) %>%
    dplyr::ungroup()



####################
## Combine Scores ##
####################

scores_treatment <- meta_scores %>%
    dplyr::left_join(
        viral_intensity_scores_treatment %>%
        dplyr::select(
            Condition,
            Remdesivir_Concentration,
            Lactoferrin_Concentration,
            cell_count_treatment_mean,
            cell_count_treatment_sem,
            nucleoli_count_treatment_mean,
            nucleoli_count_treatment_sem,
            infected_count_treatment_mean,
            infected_count_treatment_sem),
        by=c("Condition", "Remdesivir_Concentration", "Lactoferrin_Concentration")) %>%
    dplyr::left_join(
        stratominer_scores_treatment %>%
            dplyr::select(
                Condition,
                Remdesivir_Concentration,
                Lactoferrin_Concentration,
                stratominer_score_treatment),
        by=c("Condition", "Remdesivir_Concentration", "Lactoferrin_Concentration")) %>%
    dplyr::left_join(
        rf_scores_treatment %>%
            dplyr::select(
                Remdesivir_Concentration,
                Lactoferrin_Concentration,
                rf_score_treatment_mean,
                rf_score_treatment_low,
                rf_score_treatment_high),        
        by=c("Remdesivir_Concentration", "Lactoferrin_Concentration")) %>%
    dplyr::mutate(
        rem_label = rem_label %>% as.factor() %>% reorder(Remdesivir_Concentration)) %>%
    dplyr::filter(Condition == "Treatment")



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
                 ymin=cell_count_treatment_mean-cell_count_treatment_sem,
                 ymax=cell_count_treatment_mean+cell_count_treatment_sem,
                 group=Remdesivir_Concentration),
               color=score_types['cell_count'],
               width=.1) +
            ggplot2::geom_line(
               mapping=ggplot2::aes(
                   x=Lactoferrin_Concentration*10^-6,
                   y=cell_count_treatment_mean,
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
                 ymin=nucleoli_count_treatment_mean-nucleoli_count_treatment_sem,
                 ymax=nucleoli_count_treatment_mean+nucleoli_count_treatment_sem,
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
                 ymin=infected_count_treatment_mean-infected_count_treatmetn_sem,
                 ymax=infected_count_treatment_mean+infected_count_treatment_sem),
               color=score_types['viral_intensity'],
               width=.1) +
            ggplot2::geom_line(
               mapping=ggplot2::aes(
                   x=Lactoferrin_Concentration*10^-6,
                   y=infected_count_treatment_mean),
               color=score_types['viral_intensity'],
               size=1)
    }
    if("stratominer_score" %in% names(score_types)){
        cat("Adding stratominer score layer with color ", score_types['stratominer_score'], "\n", sep="")        
        plot <- plot +
            ggplot2::geom_line(
               mapping=ggplot2::aes(
                   x=Lactoferrin_Concentration*10^-6,
                   y=stratominer_score_treatment),
               color=score_types['stratominer_score'],
               size=1)
    }
    if("rf_score" %in% names(score_types)){
        cat("Adding rf score layer with color ", score_types['rf_score'], "\n", sep="")        
        plot <- plot +
            ggplot2::geom_errorbar(
               mapping=ggplot2::aes(
                 x=Lactoferrin_Concentration*10^-6,
                 ymin=rf_score_treatment_low,
                 ymax=rf_score_treatment_high),
               color=score_types['rf_score'],
               width=.1) +
            ggplot2::geom_line(
               mapping=ggplot2::aes(
                   x=Lactoferrin_Concentration*10^-6,
                   y=rf_score_treatment_mean),
               color=score_types['rf_score'],
               size=1)
    }
    plot <- plot +
        facet_wrap(~rem_label)+
        ggplot2::ggtitle(
           label="Dose Response for Lactoferrin vs. Remdesivir",
           subtitle=subtitle) +
        ggplot2::scale_x_log10(
           "Lactoferrin Dose (µg/mL)",
           breaks=c(0, 0.39, 0.78, 1.56, 3.12, 6.25, 12.5, 25, 50, 75, 100, 125, 150, 175, 200)*10^-6,
           labels=c("0", "0.39", "0.78", "1.56", "3.12", "6.25", "12.5", "25", "50", "75", "100", "125", "150", "175", "200")) +
        ggplot2::scale_y_continuous(
           paste0(score_types, ": ", names(score_types), collapse="   "))
}


p <- make_dose_response_plot(
    scores = scores_treatment %>% dplyr::filter(Remdesivir_Concentration == 0),
    subtitle="Plate ",
    score_types=c(
        cell_count="black",
        #nucleoli_count="green",
        rf_score='red'))
#        stratominer_score="red"))
ggplot2::ggsave(
    "product/figures/dose_response/plate_1999B_2020A/Lf_rf_score_200517.pdf",
    width=4, height=4)

p <- make_dose_response_plot(
    scores = scores_treatment,
    subtitle="Plate 1999B and 2020A",
    score_types=c(
        cell_count="black",
        #nucleoli_count="green",
        rf_score='red'))
#        stratominer_score="red"))
ggplot2::ggsave(
    "product/figures/dose_response/plate_1999B_2020A/rf_score_200517.pdf",
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
