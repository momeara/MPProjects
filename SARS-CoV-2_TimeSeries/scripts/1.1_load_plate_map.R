library(plyr)
library(magrittr)
options(tidyverse.quiet = TRUE)
library(tidyverse, quietly=TRUE)
library(readxl)
library(arrow, quietly=TRUE, warn.conflicts = FALSE)

source("parameters.R")

######################
# Time Series 202006 #
######################
cat("Loading time series plate 202006 ...\n")
plate_map_TS <-  readxl::read_excel(
        path = parameters$plate_map_fname,
        sheet = "Time Series 202006") %>%
    dplyr::rename(Compound = Compound_Name) %>%
    dplyr::mutate(
        master_plate_id = Plate_Name,
        plate_id = Plate_Name %>%
            stringr::str_replace("SARS_", ""),
        Concentration = ifelse(is.na(Concentration), 0, Concentration),
        row = Well_ID %>%
            stringr::str_extract("^[A-Z]") %>%
            purrr::map_int(~which(LETTERS==., arr.ind=T)),
        column = Well_ID %>%
            stringr::str_extract("[0-9]+$") %>%
            as.integer(),
        is_control = `Condition` %in% c("PC", "NC"),
        time_point = `20200616T154655`)
save(plate_map_TS, file = "intermediate_data/plate_map_TS.Rdata")


######################
# Time Series 202008 #
######################
cat("Loading time series plate 202008 ...\n")
plate_map_TS_202008 <-  readxl::read_excel(
        path = parameters$plate_map_fname,
        sheet = "Time Series 202008") %>%
    dplyr::rename(Compound = Compound_Name) %>%
    dplyr::mutate(
        master_plate_id = Plate_Name,
        plate_id = Plate_Name %>%
            stringr::str_replace("SARS_", ""),
        row = Well_ID %>%
            stringr::str_extract("^[A-Z]") %>%
            purrr::map_int(~which(LETTERS == ., arr.ind = T)),
        column = Well_ID %>%
            stringr::str_extract("[0-9]+$") %>%
            as.integer(),
        is_control = `Compound` %in% c("PC", "NC"),
        time_point = `Condition`)

save(plate_map_TS_202008, file = "intermediate_data/plate_map_TS_202008.Rdata")
plate_map_TS_202008 %>% arrow::write_parquet("product/plate_map_TS_202008.parquet")
