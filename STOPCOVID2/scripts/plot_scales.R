
library(ggplot2)
library(scales)


infectivity_y_scales <- list(
  ggplot2::scale_y_continuous(
    name = "% Infeected cells per well",
    labels = scales::percent_format(accuracy = 1)))

lf_single_agent_scales <- list(
  ggplot2::scale_x_continuous(
    name = "Lactoferrin (μg/mL)",
    breaks = c(1.52, 3.12, 6.25, 12.5, 25, 50, 100, 200, 400) %>% log10(),
    labels = c("0", "3.12", "6.25", "12.5", "25", "50", "100", "200", "400")))

fv_single_agent_scales <- list(
  ggplot2::scale_x_continuous(
    name = "Fluvoxamine (μM)",
    breaks = c(.09, 0.163, 0.294, 0.529, 0.953, 1.71, 3.09, 5.56, 10) %>% log10(),
    labels = c("0", "0.163", "0.294", "0.529", "0.953", "1.71", "3.09", "5.56", "10")))

checkerboard_scales <- list(
  ggplot2::scale_x_continuous(
    name = "Lactoferrin (μg/mL)",
    breaks = c(1.52, 3.12, 6.25, 12.5, 25, 50, 100, 200, 400) %>% log10(),
    labels = c("0", "3.12", "6.25", "12.5", "25", "50", "100", "200", "400"),
    expand = c(0, 0)),
  ggplot2::scale_y_continuous(
    name = "Fluvoxamine (μM)",
    breaks = c(.09, 0.163, 0.294, 0.529, 0.953, 1.71, 3.09, 5.56, 10) %>% log10(),
    labels = c("0", "0.163", "0.294", "0.529", "0.953", "1.71", "3.09", "5.56", "10"),
    expand = c(0, 0)))
