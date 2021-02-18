
library(plyr)
library(tidyverse)

library(googlesheets4)

source("parameters.R")
source("scripts/get_dataset.R")


# gather the dataset summary and ids from google drive
# HLOs/datasets
datasets <- googlesheets4::read_sheet(
    ss = parameters$datasets_googlesheet_id)

datasets %>%
    readr::write_tsv("raw_data/datasets_20210208.tsv")
