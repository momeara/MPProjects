library(tidyverse)
library(brms)



stan_code <- "
data {
  int Y[10];  // response variable
  int trials[10];  // number of trials
}
parameters {
  vector[1] b_top;  // population-level effects
}
model {
  vector[10] mu;
  for (n in 1:10) {
    // compute non-linear predictor values
    mu[n] = b_top;
  }
  target += binomial_lpmf(Y | trials, mu);
}
"

model <- rstan::stan_model(
  model_name="flat_model",
  model_code=stan_code)

fit <- rstan::sampling(
  object=model,
  data=list(Y=rbinom(10, 100, .5), t=rep(100, 10)))


z <- brms::brm(
  formula=y|trials(t) ~ 1,
  data=data.frame(y=rbinom(10, 100, .5), t=rep(100, 10)),
  family=binomial("logit"),
  inits=list(
    list(Intercept=.5),
    list(Intercept=.5),
    list(Intercept=.5),
    list(Intercept=.5)))


z <- brms::brm(
  formula=brms::brmsformula(
    y|trials(t) ~ top,
    top ~ 1,
    nl=TRUE),
  data=data.frame(y=rbinom(10, 100, .5), t=rep(100, 10)),
  family=binomial("identity"),
  prior=c(
    brms::prior(uniform(0, 1), nlpar="top", lb=0, ub=1)),
  init_r=100,
  inits=list(
    list(Intercept=.5),
    list(Intercept=.5),
    list(Intercept=.5),
    list(Intercept=.5)))

# init works
z <- brms::brm(
  formula=y|trials(t) ~ 1,
  data=data.frame(y=rbinom(10, 100, .5), t=rep(100, 10)),
  family=binomial("identity"),
  inits=list(
    list(Intercept=.5),
    list(Intercept=.5),
    list(Intercept=.5),
    list(Intercept=.5)))

# 
z2 <- brms::brm(
  formula= brms::brmsformula(
    y|trials(t) ~ inv_logit(z),
    z ~ 1,
    nl=TRUE),
  data=data.frame(y=rbinom(10, 100, .5), t=rep(100, 10)),
  prior = c(brms::prior(normal(0, 11), nlpar="z")),
  family=binomial("identity"))


z2_code <- "
// generated with brms 2.11.0
functions {
}
data {
  int<lower=1> N;  // number of observations
  int Y[N];  // response variable
  int trials[N];  // number of trials
  int<lower=1> K_value;  // number of population-level effects
  matrix[N, K_value] X_value;  // population-level design matrix
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
}
parameters {
  vector[K_value] b_value;  // population-level effects
}
transformed parameters {
}
model {
  // initialize linear predictor term
  vector[N] nlp_value = X_value * b_value;
  // initialize non-linear predictor term
  vector[N] mu;
  for (n in 1:N) {
    // compute non-linear predictor values
    mu[n] = nlp_value[n];
  }
  // priors including all constants
  target += uniform_lpdf(b_value | 0, 1);
  // likelihood including all constants
  if (!prior_only) {
    target += binomial_lpmf(Y | trials, mu);
  }
}
generated quantities {
}
"

z2_model <- rstan::stan_model(
  model_name="z2_model",
  model_code=z2_code)

fit <- rstan::sampling(
  object=z2_model,
  data=list(
    N=10,
    Y=rbinom(10, 100, .5),
    trials=rep(100, 10),
    K_value=1,
    X_value=matrix(rep(1, 10)),
    prior_only=0),
  init=list(
    chain1=list(b_value=list(.3)),
    chain2=list(b_value=list(.3)),
    chain3=list(b_value=list(.3)),
    chain4=list(b_value=list(.3))))
    

####

z3_code <- "
data {
  int Y[10];  // response variable
  int trials[10];  // number of trials
}
parameters {
  vector[1] b_value;
}
model {
  vector[10] mu;
  for (n in 1:10) {
    mu[n] = b_value[1];
  }
  target += binomial_lpmf(Y | trials, mu);
}
"
z3_model <- rstan::stan_model(
  model_name="z3_model",
  model_code=z3_code)

z3_fit <- rstan::sampling(
  object=z3_model,
  data=list(
    Y=rbinom(10, 100, .5),
    trials=rep(100, 10)),
  chains=1,
  init_r=100,
  init=list(
    chain1=list(
      b_value=as.array(.3))))




#hyper parameters
n_trials = 100      # measurements are n bernoulli trials
n_doses = 10
log_dose = 1:10

# true parameters
# Drug A is inactive
drugA_prob <- function(log_dose){.01}

# DrugB is active and has a sigmoidal dose-response
drugB_prob <- function(log_dose){
  top <- .95
  bottom <- .01
  slope <- 1
  ec50 <- 5
  (top-bottom)/(1 + exp(-slope*(log_dose - ec50)))
}

data <- data.frame(n_trials = n_trials, log_dose=log_dose) %>%
  dplyr::mutate(
    drugA = rbinom(n=n_doses, size=n_trials, prob=drugA_prob(log_dose)),
    drugB = rbinom(n=n_doses, size=n_trials, prob=drugB_prob(log_dose))) %>%
  tidyr::pivot_longer(
    cols=tidyr::starts_with("drug"), names_to="treatment", values_to="response")

model <- brms::brm(
  formula= brms::brmsformula(
    response | trials(n_trials) ~ top * inv_logit(hill*4/top*(log_dose - ec50)),
    top + ec50 + hill ~ treatment,
    nl=TRUE),
  data=data,
  prior = c(
    brms::prior(normal(1, 5), nlpar="top"),
    brms::prior(normal(2, 5), nlpar="ec50"),
    brms::prior(horseshoe(), nlpar="hill")),
  family = binomial("identity"),
  iter=2000,
  init=c(top=1, bottom=0, ic50=2, slope=1),
  control=list(
    adapt_delta=0.99,
    max_treedepth=12))


