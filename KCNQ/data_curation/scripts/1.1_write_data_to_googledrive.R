

library(plyr)
library(tidyverse)
library(googlesheets4)


date_code <- "10210221"

activities_summary <- readr::read_tsv(
    paste0("raw_data/activities_summary_", date_code, ".tsv"))
# careful writing back to make sure it's formatted correctly


substances %>% googlesheets4::write_sheet(
    ss = parameters$project_data_googlesheets_id,
    sheet = "Substances")

