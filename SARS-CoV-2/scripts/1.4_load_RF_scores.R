

library(plyr)
library(tidyverse)
library(ggplot2)
library(readxl)
library(MPStats)
library(arrow)

rf_scores_well_10XX <- dplyr::bind_rows(
    readr::read_csv("raw_data/RF_scores_well_200508/1001_ProbPos_Field.csv") %>%
        dplyr::mutate(Concentration = as.numeric(Concentration)) %>%
        dplyr::rename(Well_ID = `First(WellID)`),
    readr::read_csv("raw_data/RF_scores_well_200508/1002_ProbPos_Field.csv") %>%
        dplyr::mutate(Concentration = as.numeric(Concentration)) %>%
        dplyr::rename(Well_ID = `First(FDA_384_Well_ID)`),
    readr::read_csv("raw_data/RF_scores_well_200508/1003_ProbPos_Field.csv") %>%
        dplyr::mutate(Concentration = as.numeric(Concentration)) %>%
        dplyr::rename(Well_ID = `First(FDA_384_Well_ID)`),
    readr::read_csv("raw_data/RF_scores_well_200508/1004B_ProbPos_Field.csv") %>%
        dplyr::mutate(Concentration = as.numeric(Concentration)) %>%
        dplyr::rename(Well_ID = `First(FDA_384_Well_ID)`),
    readr::read_csv("raw_data/RF_scores_well_200508/1005B_ProbPos_Field.csv") %>%
        dplyr::mutate(Concentration = as.numeric(Concentration)) %>%
        dplyr::rename(Well_ID = `First(FDA_384_Well_ID)`),
    readr::read_csv("raw_data/RF_scores_well_200508/1005C_ProbPos_Field.csv") %>%
        dplyr::mutate(Concentration = as.numeric(Concentration)) %>%
        dplyr::rename(Well_ID = `First(FDA_384_Well_ID)`)) %>%
    dplyr::rename(
        Plate_Name = Barcode,
        cell_count_well = `Sum(Image_Count_Cells)`,
        infectivity_probpos_well_median = `Median(Image_Classify_Positive_PctObjectsPerBin)`,
        infectivity_probpos_well_sd = `Standard deviation(Image_Classify_Positive_PctObjectsPerBin)`,
        infectivity_probpos_well_PC_median = `Median(P (COND=PC))`,
        infectivity_probpos_well_PC_sd = `Standard deviation(P (COND=PC))`,
        infectivity_probpos_well_cov = `median-stddev`,
        infectivity_probpos_well_NC_median = `Median(P (COND=NC))`) %>%
    dplyr::mutate(
        master_plate_id = Plate_Name %>%
            stringr::str_extract("SARS_....") %>%
            stringr::str_replace("SARS_", "") %>%
            as.numeric(),
        dose_nM = Concentration,
        row = Well_ID %>%
            stringr::str_extract("^[A-Z]") %>%
            purrr::map_int(~which(LETTERS==., arr.ind=T)),
        column = Well_ID %>%
            stringr::str_extract("[0-9]+$") %>%
            as.integer()) 
rf_scores_well_10XX %>% save(file="intermediate_data/rf_scores_well_10XX.Rdata")


rf_scores_field_10XX <- dplyr::bind_rows(
    readr::read_csv("raw_data/RF_scores_field_200509/1001A_ProbPos_Field_Writeout.csv") %>%
        dplyr::mutate(Concentration = as.numeric(Concentration)),
    readr::read_csv("raw_data/RF_scores_field_200509/1002A_ProbPos_Field_Writeout.csv") %>%
        dplyr::mutate(Concentration = as.numeric(Concentration)),
    readr::read_csv("raw_data/RF_scores_field_200509/1003A_ProbPos_Field_Writeout.csv") %>%
        dplyr::mutate(Concentration = as.numeric(Concentration)),
    readr::read_csv("raw_data/RF_scores_field_200509/1004B_ProbPos_Field_Writeout.csv") %>%
        dplyr::mutate(Concentration = as.numeric(Concentration)),
    readr::read_csv("raw_data/RF_scores_field_200509/1005B_ProbPos_Field_Writeout.csv") %>%
        dplyr::mutate(Concentration = as.numeric(Concentration)),
    readr::read_csv("raw_data/RF_scores_field_200509/1005C_ProbPos_Field_Writeout.csv") %>%
        dplyr::mutate(Concentration = as.numeric(Concentration))) %>%
    dplyr::rename(
        Plate_Name = Barcode,
        Well_ID = WellID,
        cell_count_field = Image_Count_Cells,
        infectivity_probpos_field = `P (COND=PC)`) %>%
    dplyr::mutate(
        master_plate_id = Plate_Name %>%
            stringr::str_extract("SARS_....") %>%
            stringr::str_replace("SARS_", "") %>%
            as.numeric(),
        dose_nM = Concentration,
        row = Well_ID %>%
            stringr::str_extract("^[A-Z]") %>%
            purrr::map_int(~which(LETTERS == ., arr.ind = T)),
        column = Well_ID %>%
            stringr::str_extract("[0-9]+$") %>%
            as.integer())
rf_scores_field_10XX %>%
    save(file = "intermediate_data/rf_scores_field_10XX.Rdata")

###############################
## RF Scores for 20XX Series ##
###############################

rf_scores_well_20XX <- dplyr::bind_rows(
    readr::read_csv("raw_data/RF_scores_well_200508/2006_2010_ProbPos_No2008_Field.csv") %>%
        dplyr::transmute(
            Plate_Name = Barcode, 
            Well_ID = WellID,            
            Compound=Compound_Name,
            dose_nM = Concentration*1000,
            cell_count_well = Count,
            infectivity_probpos_well_PC_median = ProbPos),
    readr::read_csv("raw_data/RF_scores_well_200508/2011_2015_ProbPos_Field.csv") %>%
        dplyr::transmute(
            Plate_Name = Barcode,
            Well_ID = WellID,
            Compound=Compound_Name,
            dose_nM = Concentration*1000,   
            cell_count_well = `Sum(Image_Count_Cells)`,
            infectivity_probpos_well_PC_median = `Median(P (COND=PC))`)) %>%
    dplyr::mutate(
        master_plate_id = Plate_Name %>%
            stringr::str_extract("^....") %>%
            as.numeric(),   
        row = Well_ID %>%
            stringr::str_extract("^[A-Z]") %>%
            purrr::map_int(~which(LETTERS==., arr.ind=T)),
        column = Well_ID %>%
            stringr::str_extract("[0-9]+$") %>%
            as.integer()) 
rf_scores_well_20XX %>% save(file="intermediate_data/rf_scores_well_20XX.Rdata")

rf_scores_field_20XX <- dplyr::bind_rows(
    readr::read_csv("raw_data/RF_scores_field_200509/2006_2010_ProbPos_Field_Writeout.csv") %>%
       dplyr::rename(Compound = Compound_Name),
    readr::read_csv("raw_data/RF_scores_field_200509/2011_2015_ProbPos_Field_Writeout.csv") %>%
       dplyr::rename(Compound = Compound_Name),
    readr::read_csv("raw_data/RF_scores_field_200509/2016_2019_ProbPos_Field_Writeout.csv")) %>%
        dplyr::transmute(
            Plate_Name = Barcode,
            Well_ID = WellID,
            Field_ID = Image_Metadata_Field,
            Compound,
            dose_nM = Concentration*1000,   
            cell_count_field = Image_Count_Cells,
            infectivity_probpos_field = `P (COND=PC)`) %>%
    dplyr::mutate(
        master_plate_id = Plate_Name %>%
            stringr::str_extract("^....") %>%
            as.numeric(),   
        row = Well_ID %>%
            stringr::str_extract("^[A-Z]") %>%
            purrr::map_int(~which(LETTERS==., arr.ind=T)),
        column = Well_ID %>%
            stringr::str_extract("[0-9]+$") %>%
            as.integer())
rf_scores_field_20XX %>% save(file="intermediate_data/rf_scores_field_20XX.Rdata")


##########################################
## RF Scores for 1999B and 2020A on CQ1 ##
##########################################

rf_scores_field_1999B_2020A <- dplyr::bind_rows(
    readr::read_csv("raw_data/RF_scores_field_200509/1999B_cq1_ProbPos_Field_Writeout.csv"),
    readr::read_csv("raw_data/RF_scores_field_200509/2020A_cq1_ProbPos_Field_Writeout.csv")) %>%
    dplyr::transmute(
        Plate_Name = Barcode,
        Well_ID = WellID,
        Field_ID = Image_Metadata_FieldID,
        Compound,
        Lactoferrin_Concentration,
        Remdesivir_Concentration,
        cell_count_field = Image_Count_Cells,
        infectivity_probpos_field = `P (COND=PC)`,
        master_plate_id = Plate_Name %>%
            stringr::str_extract("^....") %>%
            as.numeric(),   
        row = Well_ID %>%
            stringr::str_extract("^[A-Z]") %>%
            purrr::map_int(~which(LETTERS==., arr.ind=T)),
        column = Well_ID %>%
            stringr::str_extract("[0-9]+$") %>%
            as.integer())
rf_scores_field_1999B_2020A %>% save(file="intermediate_data/rf_scores_field_1999B_2020A.Rdata")



#############################
## RF Scores for 2021A CX5 ##
#############################

rf_scores_field_2021A <- readr::read_csv(
       file="raw_data/RF_scores_field_200509/2021A_cx5_ProbPos_Field_Writeout.csv") %>%
        dplyr::transmute(
            Plate_Name = Barcode,
            master_plate_id = Plate_Name %>% stringr::str_extract("[0-9]+"),
            Well_ID = WellID,
            Field_ID = Image_Metadata_Field,
            Compound,
            Hydroxychloroquine_Concentration,
            Lactoferrin_Concentration,
            cell_count_field = Image_Count_Cells,
            infectivity_probpos_field = `P (COND=PC)`) %>%
    dplyr::mutate(
        master_plate_id = Plate_Name %>%
            stringr::str_extract("^....") %>%
            as.numeric(),   
        row = Well_ID %>%
            stringr::str_extract("^[A-Z]") %>%
            purrr::map_int(~which(LETTERS==., arr.ind=T)),
        column = Well_ID %>%
            stringr::str_extract("[0-9]+$") %>%
            as.integer())
rf_scores_field_2021A %>% save(file="intermediate_data/rf_scores_field_2021A.Rdata")


