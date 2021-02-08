
library(plyr)
library(tidyverse)
library(tcpl)

source("scripts/geom_indicator.R")
load("intermediate_data/well_scores.Rdata")

dose_response_data <- well_scores %>%
  dplyr::filter(!control) %>%
  dplyr::transmute(
    compound = compound,
    log_dose = log10(concentration*1000),
    value = prob_positive)

max_value <- dose_response_data %>%
  dplyr::group_by(compound) %>%
  dplyr::summarize(max_value=max(value)) %>%
  dplyr::arrange(max_value) %>%
  data.frame

dose_response_data <- dose_response_data %>%
  dplyr::semi_join(max_value %>% dplyr::filter(max_value > .5), by="compound")



fits <- dose_response_data %>%
  plyr::ddply(c("compound"), function(curve_data){
    tryCatch({
      fit <- tcpl::tcplFit(
        logc=curve_data$log_dose,
        resp=curve_data$value,
        bmad=0,
        force.fit=TRUE,
        bidirectiona=FALSE,
        verbose=TRUE
      )
      log_dose <- seq(min(curve_data$log_dose), max(curve_data$log_dose), length.out=100)
      data.frame(log_dose=log_dose) %>%
        dplyr::mutate(
          const_aic = signif(fit$cnst_aic, 2),
          const_pred = mean(curve_data$value),
          hill_aic = signif(fit$hill_aic, 2),
          hill_pred = fit$hill_tp/(1 + 10^((fit$hill_ga - log_dose)*fit$hill_gw)),
          gnls_aic = signif(fit$gnls_aic, 2),
          gnls_pred = fit$gnls_tp * 
            (1/(1 + 10^((fit$gnls_ga - log_dose)*fit$gnls_gw))) *
            (1/(1 + 10^((log_dose - fit$gnls_la)*fit$gnls_lw)))) %>%
        return()
    }, error=function(e){
      cat("Failed to fit curve for compound: ",curve_data$compound[1], "\n", sep="")
      return(data.frame())
    })
  })

p <- ggplot2::ggplot(data=dose_response_data) +
  ggplot2::theme_bw() +
  ggplot2::geom_line(data=fits, mapping=aes(x=log_dose, y=hill_pred), color="blue") +
  ggplot2::geom_line(data=fits, mapping=aes(x=log_dose, y=gnls_pred), color="brown") +
  ggplot2::geom_line(data=fits, mapping=aes(x=log_dose, y=const_pred), color="black") +
  ggplot2::geom_point(
    mapping=ggplot2::aes(
      x=log_dose,
      y=value)) +
   geom_indicator(
     data=fits %>% dplyr::distinct(compound, hill_aic),
     mapping=ggplot2::aes(
       indicator=paste0("hill_aic: ", hill_aic)),
     size=3,
     xpos="left",
     ypos="top",
     group=1) +
  geom_indicator(
    data=fits %>% dplyr::distinct(compound, gnls_aic),
    mapping=ggplot2::aes(
      indicator=paste0("gnls_aic: ", gnls_aic)),
    size=3,
    xpos="left",
    ypos="top",
    group=2) +
  geom_indicator(
    data=fits %>% dplyr::distinct(compound, const_aic),
    mapping=ggplot2::aes(
      indicator=paste0("const_aic: ", const_aic)),
    size=3,
    xpos="left",
    ypos="top",
    group=3) +
  ggplot2::ggtitle("BBBC021 screen: Dose-response (Const, Hill, Gain-Loss)") +
  ggplot2::scale_x_continuous(
    "log[Compound dose]",
    labels=function(breaks){breaks-9}) +
  ggplot2::scale_y_continuous(
    "Well-score (% taxol-like)",
    limits=c(0,1),
    labels=scales::percent_format()) +
  ggplot2::facet_wrap(~compound, scales="free_x")

ggplot2::ggsave(
  filename=paste0("product/dose_response_curves_tcpl_", MPStats::date_code(), ".pdf"),
  height=15,
  width=15)
ggplot2::ggsave(
  filename=paste0("product/dose_response_curves_tcpl_", MPStats::date_code(), ".png"),
  height=15,
  width=15)
