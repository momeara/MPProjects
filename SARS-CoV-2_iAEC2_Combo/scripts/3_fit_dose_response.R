
library(plyr)
library(tidyverse)

data_raw <- dplyr::bind_rows(
  readr::read_csv("raw_data/Plate1_20201014.csv") %>% dplyr::mutate(plate_id = "Plate1"),
  readr::read_csv("raw_data/Plate2_20201014.csv") %>% dplyr::mutate(plate_id = "Plate2"),
  readr::read_csv("raw_data/Plate3_20201014.csv") %>% dplyr::mutate(plate_id = "Plate3"),
  readr::read_csv("raw_data/Plate4_20201014.csv") %>% dplyr::mutate(plate_id = "Plate4"),
  readr::read_csv("raw_data/Plate5_20201014.csv") %>% dplyr::mutate(plate_id = "Plate5"),
  readr::read_csv("raw_data/Plate6_20201014.csv") %>% dplyr::mutate(plate_id = "Plate6")) %>%
  dplyr::rename(
    well_id = 1,
    n_positive = `Infected Cells`,
    cell_count = `Total Cells`) %>%
  dplyr::mutate(
    ray = dplyr::case_when(
      Treatment == "Negative" ~ "Negative",
      Treatment == "Positive" ~ "Positive",
      drug_concentration_1 == 0 ~ "Drug2",
      drug_concentration_2 == 0 ~ "Drug1",
      TRUE ~ "D1:D2"))


# infectivity rate in negative controls across plates
plate_nc <- data_raw %>%
  dplyr::filter(Treatment == "Negative") %>%
  dplyr::group_by(plate_id) %>%
  dplyr::summarize(
    plate_NC_percent_infected = mean(n_positive / cell_count),
    .groups = "drop")

# compute normalized percent infectivity
data <- data_raw %>%
  dplyr::left_join(plate_nc, by = "plate_id") %>%
  dplyr::mutate(
    normalized_percent_infected = (n_positive / cell_count) / plate_NC_percent_infected)

##############################
# Single agent dose response #
##############################

data_single <- data %>%
  dplyr::filter(ray != "D1:D2")

#
# normalized percent infectivity by dose
ggplot2::ggplot() +
  ggplot2::geom_hline(yintercept = 1, size=.1) +
  ggplot2::geom_point(
    data = data_single %>% dplyr::filter(ray == "Drug2"),
    mapping = ggplot2::aes(
      x = log10(drug_concentration_2) - 6,
      y = normalized_percent_infected)) +
  ggplot2::ggtitle(label = "Drug dose response iAEC2") +
  ggplot2::scale_y_continuous(
    "Normalized percent infected per well",
    labels = scales::percent_format()) +
  ggplot2::scale_x_continuous("log[Drug]") +
  ggplot2::facet_wrap(facets = dplyr::vars(drug_2))
ggplot2::ggsave(
  filename = "product/normalize_percent_infected_single_20201015.pdf",
  width = 8,
  height = 8)


data_combo <- dplyr::bind_rows(
  data %>%
    dplyr::filter(ray == "D1:D2"),
  data %>%
    dplyr::filter(ray == "D1:D2") %>%
    dplyr::rename(
      drug_1 = drug_2,
      drug_2 = drug_1,
      drug_concentration_1 = drug_concentration_2,
      drug_concentration_2 = drug_concentration_1) %>%
    dplyr::mutate(
      ray_ratio = factor(signif(drug_concentration_1 / drug_concentration_2, 2))))

################
# dose layout #
###############
ggplot2::ggplot() + 
  ggplot2::geom_point(
    data = data_combo %>% dplyr::filter(drug_1 < drug_2),
    mapping = ggplot2::aes(
      x = log10(drug_concentration_2) - 6,
      y = log10(drug_concentration_1) - 6)) +
  ggplot2::facet_grid(
    rows = dplyr::vars(drug_1),
    cols = dplyr::vars(drug_2),
    scales = "free") +
  ggplot2::scale_x_continuous("log[Row Drug]") +
  ggplot2::scale_y_continuous("log[Column Drug]")
ggplot2::ggsave(
  filename = "product/dose_layout_20201015.pdf",
  width = 15, 
  height = 15)



p <- ggplot2::ggplot() +
  ggplot2::geom_hline(yintercept = 0, size=.1) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "bottom") +
  MPStats::geom_indicator(
    data = data_combo %>% dplyr::distinct(plate_id, drug_1, drug_2),
    mapping = ggplot2::aes(
      indicator = plate_id,
      group = 1),
    size=3,
    color="grey50",
    xpos="left",
    ypos="bottom") +
  ggplot2::geom_point(
    data = data_combo,
    mapping = ggplot2::aes(
      x = log10(drug_concentration_2) - 6,
      y = log10(normalized_percent_infected),
      color = log10(drug_concentration_1) - 6)) +
  ggplot2::ggtitle(
    label = "Drug Combination in inhibition of SARS-CoV-2 iAEC2 infection") +
  ggplot2::scale_y_continuous(
    "log[Normalized percent infected per well]") +
    #labels = scales::percent_format()) +
  ggplot2::scale_x_continuous("log[Column Drug]") + 
  ggplot2::scale_color_continuous("log[Row Drug]") +
  ggplot2::facet_grid(
    rows = dplyr::vars(drug_1),
    cols = dplyr::vars(drug_2))
ggplot2::ggsave(
  filename = "product/normalize_percent_infected_combo_20201015.pdf",
  plot = p,
  width = 15,
  height = 13)

#
data %>%
  dplyr::group_by(
    Treatment, drug_2) %>%
  dplyr::do({
    data_combo <- .
    treatment <- data_combo$Treatment[1]
    drug_1 <- data_combo$drug_1[1]
    drug_2 <- data_combo$drug_2[1]
    output_dir <- paste0("product/drug_combos/", treatment, "_", drug_2)
    dir.create(output_dir, recursive=TRUE)
    cat("plotting to ", output_dir, "\n", sep = "")
    ggplot2::ggplot() +
      ggplot2::geom_point(
        data = data_combo,
        mapping = ggplot2::aes(
          x = log10(drug_concentration_2),
          y = normalized_percent_infected,
          color = ray_ratio)) +
      ggplot2::ggtitle(
        label = "Combo drug dose response iAEC2",
        subtitle = paste0(drug_1, " + ", drug_2)) +
      ggplot2::scale_y_continuous(
        "Normalized percent infected per well",
        labels = scales::percent_format()) +
      ggplot2::scale_x_log10("Drug Dose 2 (uM)") + 
      ggplot2::facet_wrap(facets = dplyr::vars(drug_2))
    ggplot2::ggsave(
      filename = paste0(output_dir, "/normalize_percent_infected_20201015.pdf"),
      width = 5,
      height = 4)
  })


data_combo %>%
  dplyr::filter(drug_1 < drug_2 ) %>%
  
    
  
E0 <- data_raw %>%
  dplyr::filter(Treatment == "Negative") %>%
  dplyr::summarize(
    E0 = median(Raw_Percent_Infected/100)) %>%
  magrittr::extract2("E0")




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


