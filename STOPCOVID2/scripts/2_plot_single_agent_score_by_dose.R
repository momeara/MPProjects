
###################################################

fv_single_agent <- fv_lf %>%
  dplyr::filter(Condition != "PC") %>%
  dplyr::filter(
    `Lactoferrin_Concentration (ug/mL)` == 0) %>%
  dplyr::select(
    dose = `Fluvoxamine_Concentration (uM)`,
    Ed = Raw_Percent_Infected) %>%
  dplyr::mutate(
    dose = ifelse(dose == 0, .09, dose))

fv_single_agent_median <- fv_single_agent %>%
  dplyr::group_by(dose) %>%
  dplyr::summarize(Ed = median(Ed)) %>%
  dplyr::ungroup()

ggplot2::ggplot() +
  ggplot2::theme_bw() +
  ggplot2::geom_jitter(
    data = fv_single_agent,
    mapping = ggplot2::aes(
      x = log10(dose),
      y = Ed),
    color = "grey40",
    height = 0,
    width = .03,
    size = 1) +
  ggplot2::geom_point(
    data = fv_single_agent_median,
    mapping = ggplot2::aes(
      x = log10(dose),
      y = Ed),
    color = "darkblue",
    size = 1.6) +
  ggplot2::geom_line(
    data = fv_single_agent_median %>% dplyr::filter(dose != .09),
    mapping = ggplot2::aes(
      x = log10(dose),
      y = Ed),
    color = "darkblue",
    size = 1.6) +
  ggplot2::scale_x_continuous(
    "Fluvoxamine (uM)",
    breaks = c(.09, 0.163, 0.294, 0.529, 0.953, 1.71, 3.09, 5.56, 10) %>% log10(),
    labels = c("0", "0.163", "0.294", "0.529", "0.953", "1.71", "3.09", "5.56", "10")) +
  ggplot2::scale_y_continuous(
    "Percent infected cells per well",
    limits = c(0, 60)) +
  ggplot2::ggtitle(
    label = "Fluvoxamine antiviral activity",
    subtitle = "Huh-7 cells, SARS-CoV-2 MOI: 10, Readout: NP fluorecences")

ggplot2::ggsave(
  filename = "/tmp/fluvoxamine_dose_response_20200930.pdf",
  width = 5,
  height =4)

ggplot2::ggsave(
  filename = "/tmp/fluvoxamine_dose_response_20200930.png",
  width = 5,
  height =4)

#############

lf_single_agent <- fv_lf %>%
  dplyr::filter(Condition != "PC") %>%
  dplyr::filter(
    `Fluvoxamine_Concentration (uM)` == 0) %>%
  dplyr::select(
    n_positive = Infected_Cell_Count,
    cell_count = Cell_Count,
    dose = `Lactoferrin_Concentration (ug/mL)`,
    Ed = Raw_Percent_Infected) %>%
  dplyr::mutate(
    dose = ifelse(dose == 0, 1.52, dose),
    log_dose = log10(dose))

lf_single_agent_median <- lf_single_agent %>%
  dplyr::group_by(dose) %>%
  dplyr::summarize(Ed = median(Ed)) %>%
  dplyr::ungroup()

Ed_NC <- lf_single_agent_median %>%
  dplyr::filter(dose == 1.52) %>%
  magrittr::extract2("Ed")

Ed_PC <- fv_lf %>%
  dplyr::filter(Condition == "PC") %>%
  dplyr::summarize(Ed_PC = median(Raw_Percent_Infected)) %>%
  magrittr::extract2("Ed_PC")
  
fit <- drc::drm(
  formula=Ed ~ log_dose,
  data=lf_single_agent %>% dplyr::filter(dose != 1.52),
  fct=drc::L.4(fixed=c(NA, NA, Ed_NC, NA)))
log_dose <- seq(log10(3.12), log10(400), length.out=100)
pred_value <- predict(fit, expand.grid(log_dose, 1))
fit_data <- data.frame(
  log_dose = log_dose,
  pred_value = pred_value,
  slope=fit$coefficients[1],
  bottom=fit$coefficients[2],
  top=Ed_NC,
  ic50=fit$coefficients[3])
  


lf_model <- brms::brm(
  formula = brms::brmsformula(
    n_positive | trials(cell_count) ~ bottom + (top-bottom) * inv_logit(hill*4/(top-bottom)*(log_dose - ic50)),
    top + bottom + hill + ic50 ~ 1,
    nl=TRUE),
  data=lf_single_agent %>% dplyr::mutate(dose = ifelse(dose==1.52, .000001, dose)),
  prior = c(
    brms::prior(normal(.39, .1), nlpar="top", lb=0.3, ub=0.8),
    brms::prior(normal(0.03, .1), nlpar="bottom", lb=0, ub=0.3),
    brms::prior(normal(1.3,  .9), nlpar="ic50"),
    brms::prior(normal(-1.3*(0.4-.03)/4,  .5), nlpar="hill", ub=-0.01)),
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
    verbose = TRUE))



ggplot2::ggplot() +
  ggplot2::theme_bw() +
  ggplot2::geom_jitter(
    data = lf_single_agent,
    mapping = ggplot2::aes(
      x = log10(dose),
      y = Ed),
    color = "grey40",
    height = 0,OB
    width = .03,
    size = 1) +
  ggplot2::geom_point(
    data = lf_single_agent_median,
    mapping = ggplot2::aes(
      x = log10(dose),
      y = Ed),
    color = "darkblue",
    size = 1.6) +
  #ggplot2::geom_line(
  #  data = lf_single_agent_median %>% dplyr::filter(dose != 1.52),
  #  mapping = ggplot2::aes(
  #    x = log10(dose),
  #    y = Ed),
  #  color = "darkblue",
  #  size = 1.6) +
  ggplot2::geom_line(
    data = fit_data,
    mapping = ggplot2::aes(
      x = log_dose,
      y = pred_value),
    color = "darkblue",
    size = 1.6) +
  ggplot2::scale_x_continuous(
    "Lactoferrin (ug/mL)",
    breaks = c(1.52, 3.12, 6.25, 12.5, 25, 50, 100, 200, 400) %>% log10(),
    labels = c("0", "3.12", "6.25", "12.5", "25", "50", "100", "200", "400")) +
  ggplot2::scale_y_continuous(
    "Percent infected cells per well",
    limits = c(0, 60)) +
  MPStats::geom_indicator(
    data = data.frame(indicator = paste0(expression("IC50 = ", 10^(fit_data$ic50[1]) %>% signif(2), " ug/mL"))),
    mapping = ggplot2::aes(
      indicator = indicator,
      xpos = "right",
      ypos = "top")) +
  ggplot2::ggtitle(
    label = "Lactoferrin antiviral activity",
    subtitle = "Huh-7 cells, SARS-CoV-2 MOI: 10, Readout: NP fluorecences")

ggplot2::ggsave(
  filename = "/tmp/lactoferrin_dose_response_20200930.pdf",
  width = 5,
  height =4)

ggplot2::ggsave(
  filename = "/tmp/lactoferrin_dose_response_20200930.png",
  width = 5,
  height =4)

