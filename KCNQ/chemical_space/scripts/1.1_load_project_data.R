
library(plyr)
library(tidyverse)
library(googlesheets4)


source("parameters.R")

activities_summary <- googlesheets4::read_sheet(
    ss = parameters$project_data_googlesheets_id,
    sheet = "Activities Summary") %>%
    dplyr::transmute(
        substance_class = `substance class`,
        substance_source = `substance source`,
        substance_name = `substance name`,
        substance_smiles = `substance smiles`,
        substance_zinc_id = `substance zinc_id`,
        substance_iupac = `substance IUPAC Name`,
        substance_binding_site = `substance binding site`,
        measurement,
        activity_KCNQ1_uM = `KCNQ1 activity (uM)`,
        activity_KCNQ1_KCNE1_uM = `KCNQ1/KCNE1 activity (uM)`,
        activity_KCNQ2_uM = `KCNQ2 activity (uM)`,
        activity_KCNQ3_uM = `KCNQ3 activity (uM)`,
        activity_KCNQ2_KCNQ3_uM = `KCNQ2/3 activity (uM)`,
        activity_KCNQ4_uM = `KCNQ4 activity (uM)`,
        activity_KCNQ5_uM = `KCNQ5 activity (uM)`,
        activity_KCNQ3_KCNQ5_uM = `KCNQ3/5 activity (uM)`,
        activity_hERG_uM = `hERG activity (uM)`,
        activity_GABAA_uM = `GABA(A) activity (uM)`,
        reference = reference,
        status,
        comment) %>%
    dplyr::mutate(
        dplyr::across(
            tidyselect::ends_with("uM"),
            function(x) {ifelse(x == "NULL", NA, as.character(x))}))
activities_summary %>% readr::write_tsv(
    file = paste0("raw_data/activities_summary_", parameters$date_code, ".tsv"))

 
# check for unique substance names
# it appears that in some papers, several compounds are reported multiple times
activities_summary %>%
    dplyr::distinct(substance_name, .keep_all=TRUE) %>%
    dplyr::count(substance_smiles) %>%
    dplyr::filter(n>1) %>%
    dplyr::left_join(activities_summary %>% dplyr::select(substance_name, substance_smiles, by="substance_name")) %>%
    data.frame()

### SUBSTANCES ###
substances <- googlesheets4::read_sheet(
    ss = parameters$project_data_googlesheets_id,
    sheet = "Substances")

substances <- activities_summary %>%
    dplyr::distinct(substance_smiles, .keep_all = TRUE) %>%
    dplyr::select(
        substance_class,
        substance_source,
        substance_name,
        substance_smiles,
        substance_zinc_id,
        substance_iupac,
        substance_binding_site) %>%
    dplyr::mutate(
        substance_dock_id = substance_name %>%
            stringr::str_replace_all(" ", "_") %>%
            stringr::str_replace_all("′", "") %>%
            stringr::str_sub(-16, -1) %>%
            stringr::str_trim())

substances %>% readr::write_tsv(
    file = paste0("raw_data/substances_", parameters$date_code, ".tsv"))



### RECEPTORS ###
receptors <- googlesheets4::read_sheet(
    ss = parameters$project_data_googlesheets_id,
    sheet = "Receptors") %>%
    dplyr::rename(
        gene_name = `Gene Name`,
        species = Species,
        mutations = Mutations,
        mutation_comment = `Mutation Comment`,
        reference = Reference,
        measurements = `Measurements`)
receptors %>% readr::write_tsv(
    file = paste0("raw_data/receptors_", parameters$date_code, ".tsv"))

### REFERENCES ###
references <- googlesheets4::read_sheet(
    ss = parameters$project_data_googlesheets_id,
    sheet = "References") %>%
    dplyr::transmute(
        reference = Reference,
        year = as.numeric(Year),
        lab = Lab,
        institution = Institution,
        todo = `To-Do`)
references %>% readr::write_tsv(
    file = paste0("raw_data/references_", parameters$date_code, ".tsv"))


### ACTIVITIES ###
activities <- googlesheets4::read_sheet(
    ss = parameters$project_data_googlesheets_id,
    sheet = "Activity Annotations") %>%
    dplyr::transmute(
        reference = Reference,
        document_location = `Document Location`,
        substance_name = `Substance Name`,
        receptor = `Receptor`,
        receptor_variant = `Receptor Variant`,
        cell_line = `Cell Line`,
        assay_type = `Assay Type`,
        assay_platform = `Assay Platform`,
        assay_version = `Assay Version`,
        assay_parameters = `Assay Parameters`,
        assay_notes = `Assay Notes`,
        fixed_dose = `Fixed Dose uM`,
        measured_value = as.character(`Measured Value`),
        measured_units = `Measured Units`,
        measurement_notes = `Measurement Notes`,
        retigabine_site = `Active at Retigabine Site`,
        retigabine_site_gold = `Retigabine Site Gold`,
        ICA27243_site = as.character(`Active at ICA-27243 Site`))

activities %>% readr::write_tsv(
    file = paste0("raw_data/activities_", parameters$date_code, ".tsv"))

