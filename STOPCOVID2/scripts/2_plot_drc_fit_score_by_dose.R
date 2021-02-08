

###########################################
# Lactoferrin vs Fluvoxamine with DRC Fit #
###########################################

data <- fv_lf %>%
  dplyr::filter(Condition != "PC") %>%
  dplyr::select(
    d1 = `Lactoferrin_Concentration (ug/mL)`,
    d2 = `Fluvoxamine_Concentration (uM)`,
    Ed = Raw_Percent_Infected,) %>%
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

fit_data <- data %>%
  plyr::ddply(c("d2"), function(df){
    cat("Fitting curve for d2=", df$d2[1], " ...\n", sep = "")
    fit <- drc::drm(
      formula = Ed ~ log_d1,
      data = df %>% dplyr::filter(d1 != 1.52),
      fct = drc::L.4(fixed=c(NA, NA, Ed_NC, NA)))
    log_d1 <- seq(log10(3.12), log10(400), length.out=100)
    pred_value <- predict(fit, expand.grid(log_d1, 1))
    fit_data <- data.frame(
      log_d1 = log_d1,
      pred_value = pred_value,
      slope=fit$coefficients[1],
      bottom=fit$coefficients[2],
      top=Ed_NC,
      ic50=fit$coefficients[3])
  })

indicator_data <- fit_data %>%
  dplyr::distinct(d2, .keep_all=TRUE) %>%
  dplyr::transmute(
    d2_label = signif(ifelse(d2 < .1, 0, d2), 1),
    indicator = paste0("Fv=", d2_label, ": IC50 = ", 10^(ic50) %>% signif(2), " ug/mL"))


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
      group = d2_label,
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
    breaks = c(.09, 3.09) %>% log10(),
    labels = c("0 uM", "3 uM"),
    limits = c(.09, 3.09) %>% log10()) +
  ggplot2::ggtitle(
    label = "Lactoferrin + Fluvoxamine antiviral activity",
    subtitle = "Huh-7 cells, SARS-CoV-2 MOI: 10, Readout: NP fluorecences")

ggplot2::ggsave(
  filename = "/tmp/lactoferrin_by_fluvoxamine_curves_dose_response_20200930.pdf",
  width = 5,
  height =4)

ggplot2::ggsave(
  filename = "/tmp/lactoferrin_by_fluvoxamine_curves_dose_response_20200930.png",
  width = 5,
  height =4)
