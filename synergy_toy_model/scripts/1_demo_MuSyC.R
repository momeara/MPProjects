
library(plyr)
library(tidyverse)
library(brms)
library(drc)


# Functional form:
# Ed = (
#        C1^h1 * C2^h2 * E0 +
#        d1^h1 * C2^h2 * E1 +
#        C1^h1 * d2^h2 * E2 +    
#        d1^h1 * d2^h2 * E3 * alpha
#      ) / (
#        C1^h1 * C2^h2 +
#        d1^h1 * C2^h2 +
#        C1^h1 * d2^h2 +    
#        d1^h1 * d2^h2 * alpha
#      )
#
# When d1=0 and d2=0:
# Ed = (
#        C1^h1 * C2^h2 * E0
#      ) / (
#        C1^h1 * C2^h2)
#    = E0
#
#
# When d1=0 and d2 -> Inf
# Ed = (C2^h2 * E0 + d2^h2 * E2) /
#      (C2^h2      + d2^h2)
#
# The terms without d2 go away:
# Ed = (d2^h2 * E2) / (d2^h2)
#    = E2
#
#
# When d1>0 and d2 -> Inf
# Ed = (C1^h1 * E2 + d1^h1 * E3 * alpha) /
#      (C1^h1 +      d1^h1      * alpha 
#
#
# When d1=0 and d2=C2
# Ed = (
#        C1^h1 * C2^h2 * E0 +
#        C1^h1 * C2^h2 * E2    
#      ) / (
#        C1^h1 * C2^h2 +
#        C1^h1 * C2^h2
#      )
#    = (E0 + E2) / 2
#
# When d1 > 0 what is the 


#
# When d1=0, what is the slope at d2=C2
# d(Ed)/d(d2)
#   =  d/d(d2)
#      (C1^h1 * C2^h2 * E0 + C1^h1 * d2^h2 * E2) /
#      (C1^h1 * C2^h2      + C1^h1 * d2^h2)
#
#Cancle the C1^h1 terms:
#   =  d/d(d2)
#      (C2^h2 * E0 + d2^h2 * E2) /
#      (C2^h2      + d2^h2)
#   
#
# distribute the derivative across the terms in the numerator
#   =  E0 * C2^h2 * [d/d(d2) 1     / (C2^h2 + d2^h2)] +
#      E2         * [d/d(d2) d2^h2 / (C2^h2 + d2^h2)]
#
#   =  E0 * C2^h2 * [h2 * d2^(h2-1) / (C2^h2 + d2^h2)^2] +
#      E2 * [C2^h2 * h2 * d2^(h2-1) / (C2^h2 + d2^h2)^2]
#
#   =  (E0 + E2) * C2^h2 * h2 * d2^(h2-1)/(C2^h2 + d2^h2)^2
#
# Evaluate at d2 = C2:
#   =  (E0 + E2) * h2 * C2^(2*h2-1) / [4*C2^(2*h2))]
#   =  h2 * (E0 + E2) / (4 * C2)   

# E


data <- expand.grid(
  log10_d1 = seq(-2,3, length.out = 12),
  log10_d2 = seq(-2,3, length.out = 12)) %>%
  dplyr::mutate(
    d1 = 10^log10_d1,
    d2 = 10^log10_d2,
    h1 = 1,
    h2 = 1,
    C1 = 5,
    C2 = 1,
    E0 = 100,
    E1 = 30,
    E2 = 20,
    beta = .5,
    E3 = min(E1, E2) - beta * min(E1, E2),
    alpha = 1.8,
    response = generate_effect(
      d1 = d1, d2 = d2,
      h1 = h1, h2 = h2,
      C1 = C1, C2 = C2,
      E0 = E0, E1 = E1, E2 = E2,
      beta = beta, alpha = alpha))

cat(
  "(",
    data$C1[1] * data$C2[1] * E0[1], " + ",
    data$C2[1] * data$E1[1], "*d1 + ",
    data$C1[1] * data$E2[1], "*d2 + ",
    data$E3[1] * data$alpha[1], "*d1*d2",
  ")/(",
    data$C1[1] * data$C2[1], " + ",
    data$C2[1], "*d1 + ",
    data$C1[1], "*d2 + ",
    data$alpha[1], "*d1*d2)\n", sep = "")


data_observed <- data %>%
  dplyr::mutate(
    Ed = response + rnorm(dplyr::n(), 0, 5),
    log10_Ed = log10(Ed))
