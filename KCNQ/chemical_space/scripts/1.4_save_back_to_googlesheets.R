
library(plyr)
library(tidyverse)
library(googlesheets4)

source("parameters.R")


substances <- readr::read_tsv("intermediate_data/substances_sanitized_20210323.tsv")


substances %>%
    googlesheets4::write_sheet(
        ss = parameters$project_data_googlesheets_id,
        sheet = "Substances")
