library(plyr)
library(magrittr)
library(DBI)
options(tidyverse.quiet = TRUE)
library(tidyverse, quietly=TRUE)
library(readxl)
library(RMySQL)
library(magrittr)
library(tictoc)
library(arrow, quietly=TRUE, warn.conflicts = FALSE)


plate_map_fname <- "raw_data/plate_map_200826.xlsx"
plate_map_1999B_fname <- "raw_data/plate_map_1999B.xlsx"
plate_map_2021A_fname <- "raw_data/plate_map_2021A.xlsx"


cat("Loading plate 999A ...\n")
plate_map_999A <- readxl::read_excel(
    path=plate_map_fname,
    sheet="999A_Metadata") %>%
    dplyr::mutate(
        master_plate_id = '999',
        plate_id = '0999A',
        Plate_Name = 'SARS_0999A') %>%
    dplyr::rename(Well_ID=WellID) %>%
    dplyr::mutate(
        row = Well_ID %>%
            stringr::str_extract("^[A-Z]") %>%
            purrr::map_int(~which(LETTERS==., arr.ind=T)),
        column = Well_ID %>%
            stringr::str_extract("[0-9]+$") %>%
            as.integer(),
        dose_nM = Concentration * 1000) %>%
    dplyr::mutate(
        COND = ifelse(COND == "JQ1", "Treatment", COND))
plate_map_999A %>% save(file="intermediate_data/plate_map_999A.Rdata")

cat("Loading top 140 hits from primary screen ...\n")
hits_10XX <- readxl::read_excel(
    path = plate_map_fname,
    sheet = "140_Compound_Hit_List") %>%
    dplyr::rename(
        master_plate_id = FDA_384_Barcode,
        Well_ID = FDA_384_Well_ID,
        Compound = Compound_Name) %>%
    dplyr::mutate(
        row = Well_ID %>%
            stringr::str_extract("^[A-Z]") %>%
            purrr::map_int(~ifelse(is.na(.), NA_integer_, which(LETTERS == ., arr.ind = T))),
        column = Well_ID %>%
            stringr::str_extract("[0-9]+$") %>%
            as.integer())
hits_10XX %>% save(file = "intermediate_data/hits_10XX.Rdata")

cat("Loading 2000 series dose response plates ...\n")
plate_map_20XX <- dplyr::bind_rows(
    readxl::read_excel(
        path = plate_map_fname,
        sheet = "2006A-2010A_Metadata") %>%
        dplyr::rename(Compound = Compound_Name),    
    readxl::read_excel(
        path = plate_map_fname,
        sheet = "2011A-2015A_Metadata") %>%
        dplyr::rename(Compound = Compound_Name),
    readxl::read_excel(
        path = plate_map_fname,
        sheet = "2016A-2019A_Metadata")) %>%
    dplyr::mutate(
        master_plate_id = Plate_Name %>%
            stringr::str_extract(".....$") %>%
            stringr::str_replace(".$", "") %>%
            as.numeric(),
        plate_id = Plate_Name %>%
            stringr::str_extract(".....$"),
        Concentration = ifelse(is.na(Concentration), 0, Concentration),
        row = Well_ID %>%
            stringr::str_extract("^[A-Z]") %>%
            purrr::map_int(~which(LETTERS==., arr.ind=T)),
        column = Well_ID %>%
            stringr::str_extract("[0-9]+$") %>%
            as.integer(),
        is_control = Condition %in% c("PC", "NC"),
        dose_nM=Concentration*1000)
save(plate_map_20XX, file="intermediate_data/plate_map_20XX.Rdata")



############################################
## Lactoferrin and Remdesivir Combo 1999B ##
############################################
cat("Loading lactoferrin Remdesivir Combo plate 1999B ... \n")
plate_map_1999B <- readxl::read_excel(plate_map_1999B_fname) %>%
    select(-tidyselect::starts_with("...")) %>%
    dplyr::mutate(
        plate_id = Plate_ID,
        row = Well_ID %>%
            stringr::str_extract("^[A-Z]") %>%
            purrr::map_int(~which(LETTERS==., arr.ind=T)),
        column = Well_ID %>%
            stringr::str_extract("[0-9]+$") %>%
            as.integer(),
        Compound = dplyr::case_when(
            Compound == "NC" ~ "Negative Control",
            Compound == "PC" ~ "Positive Control",
            TRUE ~ Compound),
        is_control = Compound %in% c("Negative Control", "Positive Control"))

# label the factor levels
rem_levels <- plate_map_1999B %>%
    dplyr::arrange(Remdesivir_Concentration) %>%
    dplyr::distinct(Remdesivir_Concentration) %>%
    magrittr::extract2("Remdesivir_Concentration")

rem_labels <- c(
    "Remdesivir: 0 nM",
    "Remdesivir: 3.2 nM",
    "Remdesivir: 6.2 nM",
    "Remdesivir: 12.4 nM",
    "Remdesivir: 25 nM",
    "Remdesivir: 48 nM",
    "Remdesivir: 100 nM",
    "Remdesivir: 200 nM")

# label the factor levels
lf_levels <- plate_map_1999B %>%
    dplyr::arrange(Lactoferrin_Concentration) %>%
    dplyr::distinct(Lactoferrin_Concentration) %>%
    magrittr::extract2("Lactoferrin_Concentration")

lf_labels <- c(
    "Lactoferrin: 0 µg/mL",
    "Lactoferrin: 0.39 µg/mL",
    "Lactoferrin: 0.78 µg/mL",
    "Lactoferrin: 1.56 µg/mL",
    "Lactoferrin: 3.12 µg/mL",
    "Lactoferrin: 6.25 µg/mL",
    "Lactoferrin: 12.5 µg/mL",
    "Lactoferrin: 25 µg/mL")   

plate_map_1999B <- plate_map_1999B %>%
    dplyr::mutate(
        rem_label = factor(
            Remdesivir_Concentration,
            levels=rem_levels,
            labels=rem_labels),
        lf_label = factor(
            Lactoferrin_Concentration,
            levels=lf_levels,
            labels=lf_labels))

save(plate_map_1999B, file="intermediate_data/plate_map_1999B.Rdata")



############################################
## Lactoferrin and Remdesivir Combo 2020A ##
############################################
cat("Loading lactoferrin remdesivir combo plate 2020A ...\n")
plate_map_2020A <- readxl::read_excel(
    path = plate_map_fname,
    sheet = "2020A_Metadata") %>%
    dplyr::mutate(
        plate_id = "2020A",
        row = Well_ID %>%
            stringr::str_extract("^[A-Z]") %>%
            purrr::map_int(~which(LETTERS==., arr.ind=T)),
        column = Well_ID %>%
            stringr::str_extract("[0-9]+$") %>%
            as.integer(),
        Compound = dplyr::case_when(
            Compound == "NC" ~ "Negative Control",
            Compound == "PC" ~ "Positive Control",
            TRUE ~ Compound),
        is_control = Compound %in% c("Negative Control", "Positive Control"))

# label the factor levels
rem_levels <- plate_map_2020A %>%
    dplyr::arrange(Remdesivir_Concentration) %>%
    dplyr::distinct(Remdesivir_Concentration) %>%
    magrittr::extract2("Remdesivir_Concentration")

rem_labels <- c(
    "Remdesivir: 0 nM",
    "Remdesivir: 0.2 nM",
    "Remdesivir: 0.4 nM",
    "Remdesivir: 0.8 nM",
    "Remdesivir: 1.6 nM",
    "Remdesivir: 3.2 nM",
    "Remdesivir: 6.2 nM",
    "Remdesivir: 12.4 nM",    
    "Remdesivir: 25 nM")

# label the factor levels
lf_levels <- plate_map_2020A %>%
    dplyr::arrange(Lactoferrin_Concentration) %>%
    dplyr::distinct(Lactoferrin_Concentration) %>%
    magrittr::extract2("Lactoferrin_Concentration")

lf_labels <- c(
    "Lactoferrin: 0 µg/mL",
    "Lactoferrin: 25 µg/mL",
    "Lactoferrin: 50 µg/mL",
    "Lactoferrin: 75 µg/mL",
    "Lactoferrin: 100 µg/mL",
    "Lactoferrin: 125 µg/mL",
    "Lactoferrin: 150 µg/mL",
    "Lactoferrin: 175 µg/mL",
    "Lactoferrin: 200 µg/mL")

plate_map_2020A <- plate_map_2020A %>%
    dplyr::mutate(
        rem_label = factor(
            Remdesivir_Concentration,
            levels=rem_levels,
            labels=rem_labels),
        lf_label = factor(
            Lactoferrin_Concentration,
            levels=lf_levels,
            labels=lf_labels))

save(plate_map_2020A, file="intermediate_data/plate_map_2020A.Rdata")


#############################################
## Lactoferrin and Chloroquine Combo 2021A ##
#############################################
cat("Loading lactoferrin chloroquine combo plate 2021A ...\n")
plate_map_2021A <- readxl::read_excel(
    path = plate_map_2021A_fname,
    sheet = "2021A_Metadata") %>%
    dplyr::mutate(
        Plate_Name = `Plate ID`,
        plate_id = "2021A",
        row = Well_ID %>%
            stringr::str_extract("^[A-Z]") %>%
            purrr::map_int(~which(LETTERS == ., arr.ind = T)),
        column = Well_ID %>%
            stringr::str_extract("[0-9]+$") %>%
            as.integer(),
        Compound = dplyr::case_when(
            Compound == "NC" ~ "Negative Control",
            Compound == "PC" ~ "Positive Control",
            TRUE ~ Compound),
        is_control = Compound %in% c("Negative Control", "Positive Control"))

# label the factor levels
hcq_levels <- plate_map_2021A %>%
    dplyr::arrange(Hydroxychloroquine_Concentration) %>%
    dplyr::distinct(Hydroxychloroquine_Concentration) %>%
    magrittr::extract2("Hydroxychloroquine_Concentration")

hcq_labels <- c(
    "HCQ: 0 nM",
    "HCQ: 80 nM",
    "HCQ: 15 nM",
    "HCQ: 310 nM",
    "HCQ: 620 nM",
    "HCQ: 1.3 µM",
    "HCQ: 2.5 µM",
    "HCQ: 5 µM",    
    "HCQ: 10 µM")

# label the factor levels
lf_levels <- plate_map_2021A %>%
    dplyr::arrange(Lactoferrin_Concentration) %>%
    dplyr::distinct(Lactoferrin_Concentration) %>%
    magrittr::extract2("Lactoferrin_Concentration")

lf_labels <- c(
    "Lactoferrin: 0 µg/mL",
    "Lactoferrin: 1.6 µg/mL",
    "Lactoferrin: 3.1 µg/mL",
    "Lactoferrin: 6.2 µg/mL",
    "Lactoferrin: 12 µg/mL",
    "Lactoferrin: 25 µg/mL",
    "Lactoferrin: 50 µg/mL",
    "Lactoferrin: 100 µg/mL",
    "Lactoferrin: 200 µg/mL")

plate_map_2021A <- plate_map_2021A %>%
    dplyr::mutate(
        hcq_label = factor(
            Hydroxychloroquine_Concentration,
            levels = hcq_levels,
            labels = hcq_labels),
        lf_label = factor(
            Lactoferrin_Concentration,
            levels = lf_levels,
            labels = lf_labels))

save(plate_map_2021A, file = "intermediate_data/plate_map_2021A.Rdata")


######################
# Time Series 202006 #
######################
cat("Loading time series plate 202006 ...\n")
plate_map_TS <-  readxl::read_excel(
        path = plate_map_fname,
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
        path = plate_map_fname,
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
