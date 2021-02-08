

fv_lf <- readr::read_csv(
  "~/Downloads/Fluvoxamine-Lactoferrin-2_Well.csv")

well_scores <- fv_lf %>%
  dplyr::filter(Condition != "PC") %>%
  dplyr::select(
    dose_1 = `Lactoferrin_Concentration (ug/mL)`,
    dose_2 = `Fluvoxamine_Concentration (uM)`,
    score = Raw_Percent_Infected) %>%
  dplyr::group_by(dose_1, dose_2) %>%
  dplyr::summarize(
    score = median(score)) %>%
  dplyr::ungroup()

MPStatas::plot_checkerboard_score_by_dose(
  well_scores = well_scores,
  title = "Fluvoxamine vs. Lactoferrin Infectivity",
  subtitle = "Huh-7, SARS-CoV-2 MOI: 10, Readout: NP fluor.",
  score_label = "Median percent cells Infected per well\n(Lower is better)",
  drug_label_1 = "Lactoferrin (ug/mL)",
  drug_label_2 = "Fluvoxamine (uM)")


#ggplot2::ggplot(data = data_fv_lf) + 
#  ggplot2::theme_bw() +
#  ggplot2::theme(
#    legend.position = "bottom",
#    panel.grid = element_blank(),
#    panel.background=element_rect(fill="grey40", colour="grey40")) +
#  ggplot2::geom_tile(
#    mapping = ggplot2::aes(
#      x = log10(d1),
#      y = log10(d2),
#      fill = Ed)) +
#  ggplot2::geom_contour(
#    mapping = ggplot2::aes(
#      x = log10(d1),
#      y = log10(d2),
#      z = Ed),
#    color ="darkorange") +
#  ggplot2::coord_fixed() +
#  ggplot2::ggtitle(
#    label="Fluvoxamine vs Lactoferrin Infectivity",
#    subtitle = "Huh-7, SARS-CoV-2 MOI: 10, Readout: NP fluor.") +
#  ggplot2::scale_fill_continuous(
#    "Median percent cells Infected per well\n(Lower is better)",
#    guide=guide_colourbar(reverse = TRUE)) +
#  ggplot2::scale_x_continuous(
#    "Lactoferrin (ug/mL)",
#    breaks = c(1.52, 3.12, 6.25, 12.5, 25, 50, 100, 200, 400) %>% log10(),
#    labels = c("0", "3.12", "6.25", "12.5", "25", "50", "100", "200", "400"),
#    expand=c(0, 0)) +
#  ggplot2::scale_y_continuous(
#    "Fluvoxamine (uM)",
#    breaks = c(.09, 0.163, 0.294, 0.529, 0.953, 1.71, 3.09, 5.56, 10) %>% log10(),
#    labels = c("0", "0.163", "0.294", "0.529", "0.953", "1.71", "3.09", "5.56", "10"),
#    expand=c(0, 0))

ggplot2::ggsave(
  filename="/tmp/Fv_vs_Lf_round_2_infectivity_20200930.pdf",
  width = 5,
  height = 5)

ggplot2::ggsave(
  filename="/tmp/Fv_vs_Lf_round_2_infectivity_20200930.png",
  width = 5,
  height = 5)



######################33

#############################

fv_lf <- readr::read_csv(
  "~/Downloads/Fluvoxamine-Lactoferrin-2_Well.csv")

data_fv_lf <- fv_lf %>%
  dplyr::filter(Condition != "PC") %>%
  dplyr::select(
    d1 = `Lactoferrin_Concentration (ug/mL)`,
    d2 = `Fluvoxamine_Concentration (uM)`,
    cell_count = Cell_Count) %>%
  dplyr::mutate(
    d1 = ifelse(d1 == 0, 1.52, d1),
    d2 = ifelse(d2 == 0, .09, d2)) %>%
  dplyr::group_by(d1, d2) %>%
  dplyr::summarize(
    cell_count = median(cell_count)) %>%
  dplyr::ungroup()


ggplot2::ggplot(data = data_fv_lf) + 
  ggplot2::theme_bw() +
  ggplot2::theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    panel.background=element_rect(fill="grey40", colour="grey40")) +
  ggplot2::geom_tile(
    mapping = ggplot2::aes(
      x = log10(d1),
      y = log10(d2),
      fill = log10(cell_count))) +
  ggplot2::geom_contour(
    mapping = ggplot2::aes(
      x = log10(d1),
      y = log10(d2),
      z = log10(cell_count)),
    color ="darkorange") +
  ggplot2::coord_fixed() +
  ggplot2::ggtitle(
    label="Fluvoxamine vs Lactoferrin Cell Count",
    subtitle="Huh-7 cells, SARS-CoV-2 MOI: 10") +
  ggplot2::scale_fill_continuous(
    "Median cell count per well",
    breaks = c(4000, 7000, 12000) %>% log10(),
    labels = c("4k", "7k", "12k")) +
  ggplot2::scale_x_continuous(
    "Lactoferrin (ug/mL)",
    breaks = c(1.52, 3.12, 6.25, 12.5, 25, 50, 100, 200, 400) %>% log10(),
    labels = c("0", "3.12", "6.25", "12.5", "25", "50", "100", "200", "400"),
    expand=c(0, 0)) +
  ggplot2::scale_y_continuous(
    "Fluvoxamine (uM)",
    breaks = c(.09, 0.163, 0.294, 0.529, 0.953, 1.71, 3.09, 5.56, 10) %>% log10(),
    labels = c("0", "0.163", "0.294", "0.529", "0.953", "1.71", "3.09", "5.56", "10"),
    expand=c(0, 0))

ggplot2::ggsave(
  filename="/tmp/Fv_vs_Lf_round_2_cell_count_20200930.pdf",
  width = 5,
  height = 5)

ggplot2::ggsave(
  filename="/tmp/Fv_vs_Lf_round_2_cell_count_20200930.png",
  width = 5,
  height = 5)
