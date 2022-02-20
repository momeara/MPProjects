
library(tidyverse)
library(BayesPharma)

source("parameters.R")


activities <- readr::read_tsv(
    file = "raw_data/activities_20210323.tsv")



# Evaluate Retigabine Baseline
activities %>%
    dplyr::filter(substance_name == "retigabine") %>%
    

