#######################################
# Lactoferrin vs Fluvoxamine BRMS fit #
# Fit the full hill function          #
#######################################


fv_lf <- readr::read_csv(
  "~/Downloads/Fluvoxamine-Lactoferrin-2_Well.csv")

data <- fv_lf %>%
  dplyr::filter(Condition != "PC") %>%
  dplyr::select(
    d1 = `Lactoferrin_Concentration (ug/mL)`,
    d2 = `Fluvoxamine_Concentration (uM)`,
    n_positive = Infected_Cell_Count,
    cell_count = Cell_Count,
    Ed = Raw_Percent_Infected) %>%
  dplyr::mutate(
    d1 = ifelse(d1 == 0, .00001, d1),
    d1_scale = ifelse(d1 == 0, 1.52, d1),
    d2 = ifelse(d2 == 0, 0.09, d2)) %>%
  dplyr::filter(
    d2 < .1 | (d2 > 3 & d2 < 3.1)) %>%
  dplyr::mutate(
    log10_d1 = log10(d1),
    log10_d1_scale = log10(d1_scale))

# normalize well cell count to the median NC cell count
baseline_cell_count <- fv_lf %>%
  dplyr::filter(Condition == "NC") %>%
  dplyr::summarize(
    x = median(Cell_Count)) %>%
  magrittr::extract2("x")

data_median <- data %>%
  dplyr::group_by(log10_d1_scale, d2) %>%
  dplyr::summarize(
    cell_count_median = median(cell_count / baseline_cell_count)*100,
    cell_count_low = quantile(x = cell_count / baseline_cell_count, probs = c(.025))*100,
    cell_count_high = quantile(x = cell_count / baseline_cell_count, probs = c(0.975))*100,
    Ed_median = median(Ed),
    Ed_low = quantile(x = Ed, probs = c(.025)),
    Ed_high = quantile(x = Ed, probs = c(0.975))) %>%
  dplyr::ungroup()


# group the data by fv treatment level
grouped_data <- data %>%
  dplyr::group_by(d2) %>%
  tidyr::nest()

model_top_bottom_ic50_hill <- brms::brm_multiple(
  formula = brms::brmsformula(
    n_positive | trials(cell_count) ~ bottom + (top-bottom) * inv_logit(hill*4/(top-bottom)*(log10_d1 - ic50)),
    top + bottom + hill + ic50 ~ 1,
    nl=TRUE),
  data=grouped_data$data,
  prior = c(
    brms::prior(normal(.39, .1), nlpar="top", lb=0.3, ub=0.8),
    brms::prior(normal(0.03, .1), nlpar="bottom", lb=0, ub=0.3),
    brms::prior(normal(1.3,  .9), nlpar="ic50"),
    brms::prior(normal(normal(1.3,  .9),  .5), nlpar="hill", ub=-0.01)),
  family=binomial("identity"),
  inits=function(){
    list(
      top=as.array(runif(1, 0.3, 1)),
      bottom=as.array(runif(1, 0.0, 0.1)),
      ic50=as.array(rnorm(1, log10(20), .9)),
      hill=as.array(rnorm(1, -1.3*(0.4-.03)/4, .5)))},
  iter=8000,
  control=list(
    adapt_delta=0.99,
    max_treedepth=12),
  cores = 2,
  stan_model_args = list(
    verbose = TRUE),
  combine = FALSE)

model <- model_top_bottom_ic50_hill


indicator_data <- tibble::tibble(
  d2 = c(0,3),
  log10_ic50 = c(1.52, 0.84),
  log10_ic50_low = c(1.47, 0.80),
  log10_ic50_high = c(1.55, 0.88)) %>%
  dplyr::mutate(
    ic50 = 10^log10_ic50,
    ic50_low = 10^log10_ic50_low,
    ic50_high = 10^log10_ic50_high) %>%
  dplyr::mutate(
    indicator = paste0(
      "Fv=", d2, ": ",
      "IC50=", format(signif(ic50, 2), nsmall=1), " ug/mL ",
      "95%=[", format(signif(ic50_low, 1), nsmall=1), ", ", format(signif(ic50_high, 1), nsmall=1), "]"))


fit_data <- expand.grid(
  log10_d1 = seq(log10(3.12), log10(400), length.out=200),
  d2 = c(0.09, 3.09)) %>%
  dplyr::mutate(
    top = dplyr::case_when(
      d2 == 0.09 ~ 0.3855906,
      d2 == 3.09 ~ 0.3994271),
    bottom = dplyr::case_when(
      d2 == 0.09 ~ 0.005852685,
      d2 == 3.09 ~ 0.0515397),
    ic50 = dplyr::case_when(
      d2 == 0.09 ~ 1.518819,
      d2 == 3.09 ~ 0.8359207),
    hill = dplyr::case_when(
      d2 == 0.09 ~ -0.1302968,
      d2 == 3.09 ~ -0.1226358),
    pred_value = (bottom + (top-bottom) * brms::inv_logit_scaled(hill*4/(top-bottom)*(log10_d1-ic50))) * 100)



ggplot2::ggplot() +
  ggplot2::theme_bw() +
  ggplot2::theme(
    legend.position = "bottom") +
  ######
  # viability
  ######
  # ggplot2::geom_errorbar(
  #   data = data_median %>% dplyr::filter(d2 == 0.09),
  #   mapping = ggplot2::aes(
  #     x = log10_d1_scale,
  #     ymin = cell_count_low,
  #     ymax = cell_count_high),
  #   color = scales::muted("darkblue")) +
  # ggplot2::geom_errorbar(
  #   data = data_median %>% dplyr::filter(d2 == 3.09),
  #   mapping = ggplot2::aes(
  #     x = log10_d1_scale,
  #     ymin = cell_count_low,
  #     ymax = cell_count_high),
  #   color = scales::muted("gold")) +
  # ggplot2::geom_line(
  #   data = data_median %>% dplyr::filter(d2 == 0.09),
  #   mapping = ggplot2::aes(
  #     x = log10_d1_scale,
  #     y = cell_count_median),
  #   color = scales::muted("darkblue")) +
  # ggplot2::geom_line(
  #   data = data_median %>%
  #     dplyr::filter(d2 == 3.09) %>%
  #     dplyr::filter(log10_d1_scale > -1.),
  #   mapping = ggplot2::aes(
  #     x = log10_d1_scale,
  #     y = cell_count_median),
  #   color = scales::muted("gold")) +
  # ggplot2::geom_point(
  #   data = data_median %>% dplyr::filter(d2 == 0.09),
  #   mapping = ggplot2::aes(
  #     x = log10_d1_scale,
  #     y = cell_count_median),
  #   color = scales::muted("darkblue")) +
  # ggplot2::geom_point(
  #   data = data_median %>% dplyr::filter(d2 == 3.09),
  #   mapping = ggplot2::aes(
  #     x = log10_d1_scale,
  #     y = cell_count_median),
  #   color = scales::muted("gold")) +
  #####
  # infectivity
  #####
  ggplot2::geom_errorbar(
    data = data_median,
    mapping = ggplot2::aes(
      x = log10_d1_scale,
      ymin = Ed_low,
      ymax = Ed_high,
      color = log10(d2)),
    height = 0,
    width = .03,
    size = 1) +
  ggplot2::geom_point(
    data = data_median,
    mapping = ggplot2::aes(
      x = log10_d1_scale,
      y = Ed_median,
      color = log10(d2)),
    size = 1.6) +
  ggplot2::geom_line(
    data = fit_data,
    mapping = ggplot2::aes(
      x = log10_d1,
      y = pred_value,
      color = log10(d2),
      group = log10(d2)),
    size = 1.6) +
  MPStats::geom_indicator(
    data = indicator_data,
    mapping = ggplot2::aes(
      indicator = indicator,
      group = d2,
      xpos = "right",
      ypos = "top")) +
  ggplot2::scale_x_continuous(
    "Lactoferrin (ug/mL)",
    breaks = c(1.52, 3.12, 6.25, 12.5, 25, 50, 100, 200, 400) %>% log10(),
    labels = c("0", "3.12", "6.25", "12.5", "25", "50", "100", "200", "400")) +
  ggplot2::scale_y_continuous(
    "Percent infected cells per well",
    limits = c(0, 60)) +
  ggplot2::scale_color_gradient(
    "Fluvoxamine",
    low = scales::muted("darkblue"),
    high = "gold",
    breaks = c(0.09, 3.09) %>% log10(),
    labels = c("0 uM", "3 uM"),
    limits = c(0.09, 3.09) %>% log10()) +
  ggplot2::ggtitle(
    label = "Lactoferrin + Fluvoxamine antiviral activity",
    subtitle = "Huh-7 cells, SARS-CoV-2 MOI: 10, Readout: NP fluorecences")

ggplot2::ggsave(
  filename = "/tmp/lactoferrin_by_fluvoxamine_curves_dose_response_brms_full_20200930.pdf",
  width = 5,
  height =4)

ggplot2::ggsave(
  filename = "/tmp/lactoferrin_by_fluvoxamine_curves_dose_response_brms_full_20200930.png",
  width = 5,
  height =4)
