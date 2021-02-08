
                

ggplot2::ggplot(data = data_observed) + 
  ggplot2::theme_bw() +
  ggplot2::theme(
    legend.position = "bottom") +
  ggplot2::geom_raster(
    mapping = ggplot2::aes(
      x = log10_d1,
      y = log10_d2,
      fill = Ed)) +
  ggplot2::geom_contour(
    mapping = ggplot2::aes(
      x = log10_d1,
      y = log10_d2,
      z = response),
    color ="darkorange") +
  ggplot2::coord_fixed() +
  ggplot2::scale_x_continuous("log(Drug 1)", expand=c(0, 0)) +
  ggplot2::scale_y_continuous("log(Drug 2)", expand=c(0, 0))

# just fit alpha with known E3
model <- brms::brm_multiple(
  formula = brms::brmsformula(
    Ed ~ (
      C1^h1 * C2^h2 * E0 +
        d1^h1 * C2^h2 * E1 +
        C1^h1 * d2^h2 * E2 +
        d1^h1 * d2^h2 * alpha * E3
    ) / (
      C1^h1 * C2^h2 +
        d1^h1 * C2^h2 +
        C1^h1 * d2^h2 +
        d1^h1 * d2^h2 * alpha),
    alpha ~ 1,
    nl = TRUE),
  data = list(data_observed),
  prior = c(
    brms::prior(uniform(0, 2), nlpar="alpha", lb=0),
    brms::prior(student_t(3, 0, 5), class = sigma)),
  iter = 8000,
  inits = function(){
    list(
      b_alpha = as.array(runif(1, 0.5, 2)))},
  combine = FALSE,
  save_model = "/tmp/model.stan",
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 12),
  cores = 2,
  stan_model_args = list(
    verbose = TRUE))

# fit alpha and E3
model <- brms::brm_multiple(
  formula = brms::brmsformula(
    log10_Ed ~ log10(
        C1 * C2 * E0 +
        d1 * C2 * E1 +
        C1 * d2 * E2 +
        d1 * d2 * E3 * alpha
      ) - log10(
        C1 * C2+
        d1 * C2+
        C1 * d2 +
        d1 * d2 * alpha),
    alpha + E3 ~ 1,
    nl = TRUE),
  data = list(data_observed),
  prior = c(
    brms::prior(normal(1, 3), nlpar="alpha", lb=0),
    brms::prior(normal(25, 20), nlpar="E3", lb=0),
    brms::prior(student_t(3, 5, 2), class = sigma)),
  family = gaussian(),
  inits = function(){
    list(
      b_alpha = as.array(runif(1, 0, 3)),
      b_E3 = as.array(runif(1, 0, 30)))},
  iter = 8000,
  cores = 2,
  #stan_model_args = list(verbose = TRUE),
  control = list(
    adapt_delta = .99,
    max_treedepth = 12),
  combine = FALSE,
  save_model = "/tmp/model.stan")




  
# simplifying to h1 = h2 = 1
# fixing E3 = .005

get_prior(
  formula = brms::brmsformula(
    Ed ~ (
      C1^h1 * C2^h2 * E0 +
        d1^h1 * C2^h2 * E1 +
        C1^h1 * d2^h2 * E2 +
        d1^h1 * d2^h2 * alpha * E3
    ) / (
      C1^h1 * C2^h2 +
        d1^h1 * C2^h2 +
        C1^h1 * d2^h2 +
        d1^h1 * d2^h2 * alpha),
    alpha + E3 ~ 1,
    nl = TRUE),
  data = list(data_observed))


fit_data <- expand.grid(
  log10_d1 = seq(-2,3, length.out = 12),
  log10_d2 = seq(-2,3, length.out = 12)) %>%
  dplyr::mutate(
    d1 = 10^log10_d1,
    d2 = 10^log10_d2,
    h1 = 1,
    h2 = 1,
    C1 = 6.97,
    C2 = 0.93,
    E0 = 100,
    E1 = 30,
    E2 = 20,
    E3 = 9.66,
    beta = (min(E1, E2) - E3)/min(E1, E2),
    alpha = 3.94,
    response = generate_effect(
      d1 = d1, d2 = d2,
      h1 = h1, h2 = h2,
      C1 = C1, C2 = C2,
      E0 = E0, E1 = E1, E2 = E2,
      beta = beta, alpha = alpha),
    Ed = response + rnorm(dplyr::n(), 0, 0.11))

ggplot2::ggplot(data = data_observed) + 
  ggplot2::theme_bw() +
  ggplot2::theme(
    legend.position = "bottom") +
  ggplot2::geom_raster(
    mapping = ggplot2::aes(
      x = log10_d1,
      y = log10_d2,
      fill = Ed)) +
  ggplot2::geom_contour(
    mapping = ggplot2::aes(
      x = log10_d1,
      y = log10_d2,
      z = response),
    color ="darkorange") +
  ggplot2::geom_contour(
    data = fit_data,
    mapping = ggplot2::aes(
      x = log10_d1,
      y = log10_d2,
      z = response),
    color ="purple") +
  ggplot2::coord_fixed() +
  ggplot2::ggtitle(
    label="Toy synergy model",
    subtitle="Alpha=1.8 beta=.5; orange=sample, purple=fit") +
  ggplot2::scale_x_continuous("log(Drug 1)", expand=c(0, 0)) +
  ggplot2::scale_y_continuous("log(Drug 2)", expand=c(0, 0))


##################################################################################

# fit alpha, E3, C1, C2
model <- brms::brm_multiple(
  formula = brms::brmsformula(
    log10_Ed ~ log10(
      C1 * C2 * E0 +
        d1 * C2 * E1 +
        C1 * d2 * E2 +
        d1 * d2 * E3 * alpha
    ) - log10(
      C1 * C2+
        d1 * C2+
        C1 * d2 +
        d1 * d2 * alpha),
    C1 + C2 + alpha + E3 ~ 1,
    nl = TRUE),
  data = list(data_observed),
  prior = c(
    brms::prior(normal(1, 3), nlpar="alpha", lb=0),
    brms::prior(normal(25, 20), nlpar="E3", lb=0),
    brms::prior(normal(0.5, 5), nlpar="C1"),
    brms::prior(normal(0.5, 5), nlpar="C2"),
    brms::prior(student_t(3, 5, 2), class = sigma)),
  family = gaussian(),
  inits = function(){
    list(
      b_alpha = as.array(runif(1, 0, 3)),
      b_E3 = as.array(runif(1, 0, 30)),
      b_C1 = as.array(runif(1, -2, 3)),
      b_C2 = as.array(runif(1, -2, 3)))},
  iter = 8000,
  cores = 2,
  #stan_model_args = list(verbose = TRUE),
  control = list(
    adapt_delta = .99,
    max_treedepth = 12),
  combine = FALSE,
  save_model = "/tmp/model.stan")

#######################################
fv_lf <- readxl::read_excel(
  "~/Downloads/LFxFluvoxamine_Processed_Huh7.xlsx",
  sheet = 3)

data_fv_lf <- fv_lf %>%
  dplyr::filter(Condition != "PC") %>%
  dplyr::select(
    d1 = `Lactoferrin_Concentration (ug/mL)`,
    d2 = `Fluvoxamine_Concentration (uM)`,
    Ed = Percent_Infection) %>%
  dplyr::mutate(
    d1 = ifelse(d1 == 0, 1.52, d1),
    d2 = ifelse(d2 == 0, .09, d2)) %>%
  dplyr::group_by(d1, d2) %>%
  dplyr::summarize(
    Ed = median(Ed)) %>%
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
      fill = Ed)) +
  ggplot2::geom_contour(
    mapping = ggplot2::aes(
      x = log10(d1),
      y = log10(d2),
      z = Ed),
    color ="darkorange") +
  ggplot2::coord_fixed() +
  ggplot2::ggtitle(
    label="Fv vs Lf") +
  ggplot2::scale_x_continuous("log(Lf (ug/mL))", expand=c(0, 0)) +
  ggplot2::scale_y_continuous("log(Fv(uM))", expand=c(0, 0))

#############################

























SYNERGY















