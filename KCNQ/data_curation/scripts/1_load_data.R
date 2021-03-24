
library(plyr)
library(tidyverse)
library(googlesheets4)
library(dm)


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



activities_summary %>%
    dplyr::distinct(substance_name, .keep_all = TRUE) %>%
    dplyr::count(substance_smiles) %>%
    dplyr::filter(n > 1) %>%
    dplyr::left_join(
        activities_summary %>% dplyr::select(
            substance_name,
            substance_smiles,
            by = "substance_name")) %>%
    data.frame()


substances <- activities_summary %>%
    dplyr::distinct(substance_smiles, .keep_all = TRUE) %>%
    dplyr::select(
        substance_class,
        substance_source,
        substance_name,
        substance_smiles,
        substance_zinc_id,
        substance_iupac,
        substance_binding_site)

substances %>%
    googlesheets4::write_sheet(
        ss = parameters$project_data_googlesheets_id,
        sheet = "Substances")


substances %>% readr::write_tsv(
    file = paste0("raw_data/substances_", parameters$date_code, ".tsv"))


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


activities <- googlesheets4::read_sheet(
    ss = parameters$project_data_googlesheets_id,
    sheet = "Activity Annotations") %>%
    dplyr::transmute(
        reference,
        substance_name = `substance name`,
        receptor,
        cell_line = `cell line`,
        measurement,
        assay_params = `assay params`,
        value = as.character(value),       # todo, move non-numeric scores to info
        units,
        significant,
        info = INFO,
        retigabine_site = `Retigabine Site`,
        ICA27243_site = as.character(`ICA-27243 Site`))

activities %>% readr::write_tsv(
    file = paste0("raw_data/activities_", parameters$date_code, ".tsv"))


# check primary keys
library(dm)

kcnq_data <- dm::dm(
    activities_summary,
    receptors,
    references,
    substances,
    activities) %>%
    dm::dm_add_pk(table = activities_summary, columns = substance_name) %>%
    dm::dm_add_pk(table = references, columns = reference) %>%
    dm::dm_add_fk(table = activities_summary, columns = reference, ref_table = references) %>%
    dm::dm_add_fk(table = activities, columns = reference, ref_table = references)


kcnq_data %>%
    dm::dm_examine_constraints()

kcnq_database <- DBI::dbConnect(
    RSQLite::SQLite(),
    dbname = "intermediate_data/kcnq_database.sqlite3")

# this fails because of contraint fails
kcnq_data %>%
    dm::copy_dm_to(
        dest = kcnq_database,
        dm = kcnq_data,
        temporary = FALSE)
