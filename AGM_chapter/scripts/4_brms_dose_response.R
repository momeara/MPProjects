library(MPStats)


library(plyr)                   # data manipulation
library(tidyverse)              # data manipulation
library(rstan)                  # interface to stan
library(brms)                   # wrapper for rstan for
library(bayesplot)              # plotting bayesian models
library(broom.mixed)            # tidying bayesian models
library(tidybayes)              # tidying bayesian models

options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = TRUE)


load("intermediate_data/well_scores.Rdata")


### response model
# The target response is prob_positive, the fraction of treatment-like cells in the well.
#
#   prob_positive ~ sigmoid(log_dose)
#
# A challenge is that when there are few cells in the well, the estimate becomes unreliable.
# How to model that uncertainty? Binomial regression: The number of successes n in a series
# of N independent Bernoulli trials, where each trial has probability of success theta
#
#   binomial(n | well_cell_count, prop_positive) ~ sigmoid(log_dose)
#
# to model the sigmoid function we rescale the function inv_logit(x) = 1/(1+exp(-x)),
# using brms's formula syntax, with parameters top01, slope, and ec50:
#
#   n | trials(well_cell_count) ~ top01 * inv_logit(slope(log_dose - ec50))
#
# With this parameterization, the slope changes when the top changes.
# To improve sampling we can decouple with slope=hill*4/top. This comes from evaluating derivative at the ec50
#   
#               hill := d sigmoid / d log_dose |_{log_dose = ec50}
#                     = slope * top01 * exp(log_dose - ec50) / (1 + exp(slope * (log_dose - ec50))^2
#                     = slope * top01 * 1 / (1 + 1)^2
#                     = slope * top01 / 4
#      hill * 4 / top01 = slope
#
# since the success probability is [0, 1], the top parameter should be bounded to [0, 1] as well.
# We can do this by again using the inv_logit function
#
#          logit(top01) = top
#                 top01 = inv_logit(top)
#
# together this gives
#
#   n | trials(well_cell_count) ~ inv_logit(top) * inv_logit(hill * 4 / inv_logit(top) * (log_dose - ec50))
#

flat_model <- MPStats::fit_brms_flat_score_by_dose(
  well_scores=well_scores %>% dplyr::filter(!is_control))
save(flat_model, file="intermediate_data/flat_model.Rdata")

hill_model <- MPStats::fit_brms_hill_score_by_dose(
  well_scores=well_scores %>% dplyr::filter(!is_control))

save(hill_model, file="intermediate_data/hill_model.Rdata")

gainloss_model <- MPStats::fit_brms_gainloss_score_by_dose(well_scores=well_scores)
save(gainloss_model, file="intermediate_data/gainloss_model.Rdata")



plot_all <- MPStats::plot_brms_score_by_dose(
  well_scores=well_scores,
  flat_model=flat_model,
  hill_model=hill_model,
  subtitle="Human MCF7 cells -- compound-profiling experiment (BBBC021v1)")
ggplot2::ggsave(
  plot=plot_all,
  filename=paste0("product/brms_score_by_dose_", MPStats::date_code(), ".pdf"),
  height=20,
  width=20)
ggplot2::ggsave(
  plot=plot_all,
  filename=paste0("product/brms_score_by_dose_", MPStats::date_code(), ".png"),
  height=20,
  width=20)
#######



dose_response_data <- well_scores %>%
  dplyr::filter(!control) %>%
  dplyr::mutate(
    log_dose = log10(concentration*1000),
    n_positive = as.integer(prob_positive * well_cell_count))

# brms can fit a model against multiple data at once
dose_response_data_by_compound <- dose_response_data %>%
  dplyr::group_by(compound) %>%
  tidyr::nest()

flat_model <- brms::brm_multiple(
  formula= n_positive | trials(well_cell_count) ~ 0 + Intercept,
  data=dose_response_data_by_compound$data,
  family = binomial("logit"),
  prior = c(brms::prior(normal(0, 10), coef="Intercept")),
  combine=FALSE)
loo_flat <- add_criterion(flat_model, c("loo", "bayes_R2"))
save("flat_model", file="intermediate_data/flat_model.Rdata")

loo_flat <- add_criterion(flat_model, c("loo", "bayes_R2"))


hill_model_bounds <- brms::brm_multiple(
  formula = brms::brmsformula(
    n_positive | trials(well_cell_count) ~ top * inv_logit(hill*4/top*(log_dose - ic50)),
    top + ic50 + hill ~ 1,
    nl=TRUE),
  data=dose_response_data_by_compound$data,
  prior = c(
    brms::prior(uniform(0, 1), nlpar="top", lb=0, ub=1),
    brms::prior(normal(2, 5), nlpar="ic50"),
    brms::prior(normal(1, 5), nlpar="hill")),
  family=binomial("identity"),
  inits=function(){
    list(
      top=as.array(runif(1, 0, 1)),
      ic50=as.array(rnorm(1, 2, 5)),
      hill=as.array(rnorm(1, 1, 5)))},
  iter=4000,
  control=list(
    adapt_delta=0.99,
    max_treedepth=12),
  combine=FALSE)
save("hill_model_bounds", file="intermediate_data/hill_model_bounds.Rdata")

# summarize fits


fit_summary <- plyr::ldply(1:length(dose_response_data_by_compound$compound), function(i){
  cat("summarizing model for compound ", i, "\n", sep="")
  compound <- dose_response_data_by_compound$compound[[i]]
  flat_fit <- flat_model[[i]] %>%
    brms::add_criterion(criterion=c("loo"), model_name="flat_model")
  
  flat_summary <- flat_fit %>%
    tidybayes::spread_draws(b_Intercept) %>%
    tidybayes::median_qih(b_Intercept)
  
  hill_fit <- hill_model[[i]] %>%
    brms::add_criterion(criterion=c("loo"), model_name="hill_model")

  hill_summary <- hill_fit %>%
    tidybayes::spread_draws(b_top_Intercept, b_ic50_Intercept, b_hill_Intercept) %>%
    tidybayes::median_qih(b_top_Intercept, b_ic50_Intercept, b_hill_Intercept)
  
  tibble::tibble(
    compound=compound,
    flat_elpd_loo = flat_fit$criteria$loo$estimates[1,1],
    flat_se_elpd_loo = flat_fit$criteria$loo$estimates[1,2],
    flat_p_loo = flat_fit$criteria$loo$estimates[2,1],
    flat_se_p_loo = flat_fit$criteria$loo$estimates[2,2],
    flat_value_median_hqi = flat_summary$b_Intercept %>% brms::inv_logit_scaled(),
    flat_value_median_hqi_low95 = flat_summary$.lower %>% brms::inv_logit_scaled(),
    flat_value_median_hqi_high95 = flat_summary$.upper %>% brms::inv_logit_scaled(),
    hill_elpd_loo = hill_fit$criteria$loo$estimates[1,1],
    hill_se_elpd_loo = hill_fit$criteria$loo$estimates[1,2],
    hill_p_loo = hill_fit$criteria$loo$estimates[2,1],
    hill_se_p_loo = hill_fit$criteria$loo$estimates[2,2],
    hill_top_median_qih = hill_summary$b_top_Intercept %>% brms::inv_logit_scaled(),
    hill_top_median_qih_low95 = hill_summary$b_top_Intercept.lower %>% brms::inv_logit_scaled(),
    hill_top_median_qih_high95 = hill_summary$b_top_Intercept.upper %>% brms::inv_logit_scaled(),
    hill_ic50_median_qih = hill_summary$b_ic50_Intercept,
    hill_ic50_median_qih_low95 = hill_summary$b_ic50_Intercept.lower,
    hill_ic50_median_qih_high95 = hill_summary$b_ic50_Intercept.upper,
    hill_hill_median_qih = hill_summary$b_hill_Intercept,
    hill_hill_median_qih_low95 = hill_summary$b_hill_Intercept.lower,
    hill_hill_median_qih_high95 = hill_summary$b_hill_Intercept.upper)
})

fit_summary %>%
    readr::write_tsv(
      paste0("product/dose_response_bayesian_fits_", MPStats::date_code(), ".tsv"))

z <- hill_model[[50]] %>% brms::conditional_effects("log_dose")



