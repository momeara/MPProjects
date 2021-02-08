
library(plyr)
library(tidyverse)
library(MPStats)
library(brms)
library(tidybayes)


source("scripts/plot_scales.R")


data_raw <- readr::read_csv(
  "raw_data/Fluvoxamine-Lactoferrin-2_Well_20200930.csv")


# MuSyC params:
#   E0 <- median NC score
#
#   Lactoferrin params
#   d1 <- dose ug/mL scaled so max(d1) == 1
#   C1 <- scaled IC50
#   E1 <- maximum effect for 
#   s1 <- slope at IC50 if 

E0 <- data_raw %>%
  dplyr::filter(Condition == "NC") %>%
  dplyr::summarize(
    E0 = median(Raw_Percent_Infected/100)) %>%
  magrittr::extract2("E0")

data <- data_raw %>%
  dplyr::filter(Condition != "PC") %>%
  dplyr::transmute(
    
    # plotting dose 0 on the log scale has to be handled specially,
    # set the dose to be lower than the lowest dose with a bit of a gap
    dose_1 = ifelse(d1 == 0, 1.52, `Lactoferrin_Concentration (ug/mL)`),
    dose_2 = ifelse(d2 == 0, 0.09, `Fluvoxamine_Concentration (uM)`),
    
    # Numerically, the MuSyC model is more stable if doses fall in [0-1]
    d1_scale_factor = max(`Lactoferrin_Concentration (ug/mL)`),
    d2_scale_factor = max(`Fluvoxamine_Concentration (uM)`),
    d1 = `Lactoferrin_Concentration (ug/mL)` / d1_scale_factor,
    d2 = `Fluvoxamine_Concentration (uM)` / d2_scale_factor,
    
    # model the response as bernoulli trials to down-weight
    # doses where there are very few cells
    n_positive = Infected_Cell_Count,
    cell_count = Cell_Count,
    Ed = Raw_Percent_Infected/100) %>%
  dplyr::mutate(
    # no drug params
    E0 = E0,
    # lactoferrin params
    C1 = 33 / d1_scale_factor,
    E1 = .01,
    s1 = 1,
    # fluvoxamine params
    C2 = 1 / d2_scale_factor,
    E2 = E0,
    s2 = 0.5,
    # transformed params
    h1 = s1 * (4 * C1) / (E0 + E1),
    h2 = s2 * (4 * C2) / (E0 + E2))
#   s1 = h1 * (E0 + E1) / (4 * C1),
#   s2 = h2 * (E0 + E2) / (4 * C2))

#   h1 = s1 * (4 * C1) / (E0 + E1)
#   h2 = s2 * (4 * C2) / (E0 + E2)

# setting h1 = 1 then
#   s1 = (4 * C1) / (E0 + E1) 

MuSyC_Ed_formula <- brms::brmsformula(
  n_positive | trials(cell_count) ~ (
    C1^h1 * C2^h2 * E0 +
    d1^h1 * C2^h2 * E1 +
    C1^h1 * d2^h2 * E2 +
    d1^h1 * d2^h2 * E3 * alpha
  ) / (
    C1^h1 * C2^h2 +
    d1^h1 * C2^h2 +
    C1^h1 * d2^h2 +
    d1^h1 * d2^h2 * alpha),
  brms::nlf(h1 ~ s1 * (4 * C1) / (E0 + E1)),
  brms::nlf(h2 ~ s2 * (4 * C2) / (E0 + E2)),
  C1 + E1 + s1 + C2 + E2 + s2 + alpha + E3 ~ 1,
  nl = TRUE)

MuSyC_Ed_prior <- c(
  # drug 1 params
  b_C1_Intercept = brms::prior(normal(.5, .5), nlpar = "C1", lb = 0),
  b_E1_Intercept = brms::prior(normal(.1, .2), nlpar = "E1", lb=0, ub=1),
  b_s1_Intercept = brms::prior(normal(1, 3), nlpar = "s1", lb = .1),

  # drug 2 params
  b_C2_Intercept = brms::prior(normal(.5, .5), nlpar = "C2", lb = 0),
  b_E2_Intercept = brms::prior(normal(.39, .2), nlpar = "E2", lb=0, ub=1),
  b_s2_Intercept = brms::prior(normal(1, 3), nlpar = "s2", lb = .1),

  # synergy params
  b_alpha_Intercept = brms::prior(normal(1, 3), nlpar = "alpha", lb = 0),
  b_E3_Intercept = brms::prior_string(paste0("normal(", E1, ", .3)"), nlpar="E3", lb=0))

MuSyC_Ed_inits <- function(){list(
  # drug 1 params
  b_C1 = function(){as.array(normal(1, .5, .5))},
  b_E1 = function(){as.array(.1)},
  b_s1 = function(){as.array(runif(1, .1, 1.5))},
  
  # drug 2 params
  b_C2 = function(){as.array(normal(1, .5, .5))},
  b_E2 = function(){as.array(.39)},
  b_s2 = function(){as.array(runif(1, .5, 1.5))},
  
  # synergy params
  b_alpha = function(){as.array(runif(1, 0, 3))},
  b_E3 = function(){as.array(runif(1, 0, .5))})}


brms::get_prior(
  MuSyC_Ed_formula,
  data = list(data),
  family=binomial("identity"))

model_code <- brms::make_stancode(
  formula = MuSyC_Ed_formula,
  data = list(data),
  family=binomial("identity"),
  prior = MuSyC_Ed_prior)


model <- brms::brm_multiple(
  formula = MuSyC_Ed_formula,
  prior = MuSyC_Ed_prior,
  data = list(data),
  family=binomial("identity"),
  inits = MuSyC_Ed_inits,
  iter = 8000,
  cores = 4,
  stan_model_args = list(verbose = TRUE),
  control = list(
    adapt_delta = .99,
    max_treedepth = 12),
  combine = FALSE,
  save_model = "/tmp/model.stan")

# get posterior intervals for parameters
estimated_parameters <- model[[1]] %>%
  tidybayes::spread_draws(
    b_C1_Intercept,
    b_E1_Intercept,
    b_s1_Intercept,
    b_C2_Intercept,
    b_E2_Intercept,
    b_s2_Intercept,
    b_alpha_Intercept,
    b_E3_Intercept) %>%
    tidybayes::median_qi()

model[[1]] %>% shinystan::launch_shinystan()

  
data_fitted <- data %>%
  tidybayes::add_fitted_draws(
    model = model[[1]],
    value = "fitted_Ed",
    n = 10) %>%
  dplyr::ungroup()


generate_MuSyC_effects <- function(
  d1,
  d2,
  E0,
  h1, C1, E1,
  h2, C2, E2,
  alpha,
  E3 = NULL,
  beta = NULL) {
  
  if(!is.null(beta)){
    E3 <- min(E1, E2) - beta * min(E1, E2)
  } else if(is.null(E3)){
    stop("either E3 or beta must be non-null")
  }
  
  if(d1 == -Inf & d2 == -Inf){
    response <- E0
  } else if(d2 == -Inf){
    if(d1 == Inf){
      response <- E1
    } else {
      response <-
        (C1^h1 * E0 + d1^h1 * E1) / 
        (C1^h1      + d1^h1)
    }
  } else if(d1 == -Inf){
    if(d2 == Inf){
      response <- E2
    } else {
      response <-
        (C2^h2 * E0 + d2^h2 * E2) /
        (C2^h2      + d2^h2)
    }
  } else {
    if(d1 == Inf & d2 == Inf) {
      response <- E3
    } else {
      response <-
        (
          C1^h1 * C2^h2 * E0 +
          d1^h1 * C2^h2 * E1 +
          C1^h1 * d2^h2 * E2 +
          d1^h1 * d2^h2 * E3 * alpha
        ) / (
          C1^h1 * C2^h2 +
          d1^h1 * C2^h2 +
          C1^h1 * d2^h2 +
          d1^h1 * d2^h2 * alpha)
    }
  }
  return(response)
}

data_fitted <- data %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    fitted_Ed = generate_MuSyC_effects(
      d1 = d1,
      d2 = d2,
      E0 = E0,
      C1 = estimated_parameters$b_C1_Intercept[1],
      E1 = estimated_parameters$b_E1_Intercept[1],
      h1 = estimated_parameters$b_s1_Intercept[1] * (4 * C1) / (E0 + E1),
      C2 = estimated_parameters$b_C2_Intercept[1],
      E2 = estimated_parameters$b_E2_Intercept[1],
      h2 = estimated_parameters$b_s2_Intercept[1] * (4 * C2) / (E0 + E2),
      alpha = estimated_parameters$b_alpha_Intercept[1],
      E3 = estimated_parameters$b_E3_Intercept[1]),
    .draw = 1) %>%
  dplyr::ungroup()



ggplot2::ggplot() +
  ggplot2::theme_bw() +
  ggplot2::geom_line(
    data = data_fitted,
    mapping = ggplot2::aes(
      x=log10(d1*d1_scale_factor),
      y=fitted_Ed,
      group = paste(.draw,d2*d2_scale_factor),
      color = d2*d2_scale_factor),
    size = .5) +
  lf_single_agent_scales +
  infectivity_y_scales

ggplot2::ggplot() +
  ggplot2::theme_bw() +
  ggplot2::geom_line(
    data = data_fitted %>%
      dplyr::filter(d1 == 0),
    mapping = ggplot2::aes(
      x=log10(d2),
      y=fitted_Ed)) +
  fv_single_agent_scales +
  infectivity_y_scales



MPStats::plot_checkerboard_score_by_dose(
  well_scores = data %>%
    dplyr::group_by(dose_1, dose_2) %>%
    dplyr::summarize(
      score = median(Ed)) %>%
    dplyr::ungroup()) +
  ggplot2::geom_contour(
    data = data_fitted,
    mapping = ggplot2::aes(
      x = log10(dose_1),
      y = log10(dose_2),
      z = fitted_Ed,
      group = .draw),
    color = "purple",
    size = .3) +
   ggplot2::ggtitle("Fluvoxamine vs. Lactoferrin Infectivity") +
  checkerboard_scales
  

