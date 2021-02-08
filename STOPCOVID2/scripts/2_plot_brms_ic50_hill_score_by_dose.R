#######################################
# Lactoferrin vs Fluvoxamine BRMS fit #
# Just fit the IC50                   #
#######################################

data <- fv_lf %>%
  dplyr::filter(Condition != "PC") %>%
  dplyr::select(
    d1 = `Lactoferrin_Concentration (ug/mL)`,
    d2 = `Fluvoxamine_Concentration (uM)`,
    n_positive = Infected_Cell_Count,
    cell_count = Cell_Count) %>%
  dplyr::mutate(
    d1 = ifelse(d1 == 0, 1.52, d1),
    d2 = ifelse(d2 == 0, 0.09, d2)) %>%
  dplyr::filter(
    d2 < .1 | (d2 > 3 & d2 < 3.1)) %>%
  dplyr::mutate(
    log_d1 = log10(d1))

data_median <- data %>%
  dplyr::group_by(d1, d2) %>%
  dplyr::summarize(
    Ed_median = median(Ed),
    Ed_low = quantile(x = Ed, probs = c(.025)),
    Ed_high = quantile(x = Ed, probs = c(0.975))) %>%
  dplyr::ungroup()


Ed_NC <- data_median %>%
  dplyr::filter(d1 == 1.52, d2 == 0.09) %>%
  magrittr::extract2("Ed_median")

grouped_data <- data %>%
  dplyr::filter(d1 != 1.52) %>%
  dplyr::group_by(d2) %>%
  tidyr::nest()




model_ic50_hill <- brms::brm_multiple(
  formula = brms::brmsformula(
    n_positive | trials(cell_count) ~ .4 * inv_logit(hill*(log_d1 - ic50)),
    hill + ic50 ~ 1,
    nl=TRUE),
  data=grouped_data$data,
  prior = c(
    brms::prior(normal(1.3,  1), nlpar="ic50"),
    brms::prior(normal(-1.3,  1), nlpar="hill")),
  family=binomial("identity"),
  inits=function(){
    list(
      ic50=as.array(rnorm(1, log10(20), 1)),
      hill=as.array(rnorm(1, -1.3, 1)))},
  iter=8000,
  control=list(
    adapt_delta=0.99,
    max_treedepth=12),
  cores = 2,
  stan_model_args = list(
    verbose = TRUE),
  combine = FALSE)




model <- model_ic50_hill


indicator_data <- tibble::tibble(
  d2 = c(3, 0),
  log10_ic50 = c(1.06, 1.48),
  log10_ic50_low = c(1.04, 1.47),
  log10_ic50_high = c(1.08, 1.49)) %>%
  dplyr::mutate(
    ic50 = 10^log10_ic50,
    ic50_low = 10^log10_ic50_low,
    ic50_high = 10^log10_ic50_high) %>%
  dplyr::mutate(
    indicator = paste0(
      "Fv=", d2, ": ",
      "IC50=", signif(ic50, 2), " ug/mL ",
      "95%=[", format(round(ic50_low, 1), nsmall=1), ", ", format(round(ic50_high, 1), nsmall=1), "]"))


hill_model_draws <- model %>%
  brms:::posterior_samples(pars=c("Intercept")) %>%
  dplyr::sample_n(50) %>%
  dplyr::mutate(draw_id = dplyr::row_number()) %>%
  purrr::pmap_dfr(function(b_ic50_Intercept, draw_id){
    tibble::tibble(
      compound = compound,
        log_dose = seq(log10(3.12), log10(400), length.out=2),
        Ed = 0.4 * brms::inv_logit_scaled(-1.13*(log_dose-b_ic50_Intercept)),
        draw_id = draw_id)
    })

fit_data <- expand.grid(
  log_d1 = seq(log10(3.12), log10(400), length.out=200),
  d2 = c(3.09, 0.09)) %>%
  dplyr::mutate(
    ic50 = dplyr::case_when(
      d2 == 3.09 ~ 1.06,
      d2 == 0.09 ~ 1.48),
    hill = dplyr::case_when(
      d2 == 3.09 ~ -1.01,
      d2 == 0.09 ~ -1.3),
    pred_value = 0.4 * brms::inv_logit_scaled(hill*(log_d1-ic50)) * 100)



ggplot2::ggplot() +
  ggplot2::theme_bw() +
  ggplot2::theme(
    legend.position = "bottom") +
  ggplot2::geom_errorbar(
    data = data_median,
    mapping = ggplot2::aes(
      x = log10(d1),
      ymin = Ed_low,
      ymax = Ed_high,
      color = log10(d2)),
    height = 0,
    width = .03,
    size = 1) +
  ggplot2::geom_point(
    data = data_median,
    mapping = ggplot2::aes(
      x = log10(d1),
      y = Ed_median,
      color = log10(d2)),
    size = 1.6) +
  ggplot2::geom_line(
    data = fit_data,
    mapping = ggplot2::aes(
      x = log_d1,
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
  ggplot2::scale_color_continuous(
    "Fluvoxamine",
    breaks = c(0.09, 3.09) %>% log10(),
    labels = c("0 uM", "3 uM"),
    limits = c(0.09, 3.09) %>% log10()) +
  ggplot2::ggtitle(
    label = "Lactoferrin + Fluvoxamine antiviral activity",
    subtitle = "Huh-7 cells, SARS-CoV-2 MOI: 10, Readout: NP fluorecences")

ggplot2::ggsave(
  filename = "/tmp/lactoferrin_by_fluvoxamine_curves_dose_response_brms_ic50_hill_20200930.pdf",
  width = 5,
  height =4)

ggplot2::ggsave(
  filename = "/tmp/lactoferrin_by_fluvoxamine_curves_dose_response_brms_ic50_hill_20200930.png",
  width = 5,
  height =4)
