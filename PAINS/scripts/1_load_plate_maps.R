library(plyr)
library(tidyverse)
library(arrow)

# one row per assay plate
plate_tracking <- readxl::read_xlsx(
    "raw_data/Plate\ Maps/Cell\ Painting\ Compound\ and\ Assay\ Plate\ Tracking.xlsx") %>%
    dplyr::filter(!is.na(`Assay Plate`)) %>%
    tidyr::fill(
        `Compound Plate`,
        `Compound Plate Nickname`,
        `Plate Map`) %>%
    dplyr::select(
        -`Imaging Check`,
        -`Size (GB)`,
        -`Images exported?`,
        -`Transferred to AWS?`,
        -`Processed on AWS?`,
        -`Archived on GCP?`)
plate_tracking %>%
    arrow::write_parquet(
        "intermediate_data/plate_tracking.parquet")


# one row per compound plate well
plate_map <- dplyr::bind_rows(
    readxl::read_xls("raw_data/Plate\ Maps/JLD012020-Ono20.xls"),
    readxl::read_xls("raw_data/Plate\ Maps/JLD012020-Ono21.xls"),
    readxl::read_xls("raw_data/Plate\ Maps/JLD012020-Ono22.xls"),
    readxl::read_xls("raw_data/Plate\ Maps/JLD012020-Ono23.xls"),
    readxl::read_xls("raw_data/Plate\ Maps/JLD012020-Ono24.xls"),
    readxl::read_xls("raw_data/Plate\ Maps/JLD012020-Ono25.xls"),
    readxl::read_xls("raw_data/Plate\ Maps/JLD012020-Ono26.xls"),
    readxl::read_xls("raw_data/Plate\ Maps/JLD012020-Ono27.xls"))
plate_map %>%
    arrow::write_parquet(
        "intermediate_data/plate_map.parquet")


compound_map_fname <- "raw_data/Plate\ Maps/NextGenKATinhibitors-V9-SuppData1.xlsx"
compound_map_sheets <- readxl::excel_sheets(compound_map_fname)

compound_map <- plyr::ldply(compound_map_sheets, function(sheet){
    data <- readxl::read_xlsx("raw_data/Plate\ Maps/NextGenKATinhibitors-V9-SuppData1.xlsx", sheet = sheet) %>%
        dplyr::mutate(compound_set = sheet) %>%
        data.frame()
})

compound_map <- compound_map %>%
    dplyr::mutate(
        Manuscript.number = ifelse(
            is.na(Manuscript.number),
            Manscript.number,
            Manuscript.number),
        Institutional.ID = ifelse(
            compound_set %in% c("Annexin data", "Caspases data", "CytoTox data"),
            Compound.name.1,
            Institutional.ID),
        Alias = ifelse(
            compound_set %in% c("Annexin data", "Caspases data", "CytoTox data"),
            Compound.name.2,
            Alias),
        Institutional.ID = ifelse(
            !is.na(Alias) & (Alias == "DMSO"),
            "DMSO",
            Institutional.ID),
        Normalized = ifelse(
            is.na(Normalized),
            Normalization,
            Normalized)) %>%
    dplyr::rename_at(
        vars(tidyselect::matches("X[0-9]+")),
        ~ paste0("measurement_", stringr::str_replace(., "X", ""))) %>%
    dplyr::select(
        -tidyselect::starts_with("..."),
        -Manscript.number,
        -Normalization,
        -Compound.name.1,
        -Compound.name.2,
        compound_set,
        manuscript_number = Manuscript.number,
        alias = Alias,
        class = Class,
        quality = Quality,
        inference_moa = InterferenceMOA,
        compound_id = Institutional.ID,
        purity = Purity....,
        purity_comments = Purity.comments,
        smiles = SMILES,
        cellular_injury_category = Cellular.injury.category,
        grade = Grade,
        comments = Comments,
        reactive_category = Reactive.category,
        inactive_control = Inactive.control,
        category = Category,
        correlation = Correlation,
        dose = Dose,
        plate = Plate,
        replicate = Replicate,
        assay = Assay,
        glutathione_uM = X.Glutathione.species...uM,
        glutathione_percent_control = X.Glutathione.species...Percent.control,
        cell_viablity_percent_control = Cell.viability..Percent.control,
        thiol = Thiol,
        reaction_time_min = Reaction.time..min.,
        compound_uM = Compound.concentration..uM.,
        unreacted_thiols_percent = Unreacted.thiols....,
        measurement = Measurement,
        normalized = Normalized) %>%
    dplyr::filter(!is.na(compound_id))


compound_map %>%
    arrow::write_parquet(
        "intermediate_data/compound_map.parquet")


# extract plate barcodes from the cpdata.h5 files
dataset_ids <- readr::read_tsv("raw_data/dataset_ids.tsv") %>%
    plyr::adply(1, function(df) {
        dataset_id <- df$dataset_id[1]
        cat("Getting plate barcode for dataset '", dataset_id, "' ...\n", sep = "")
        dataset <- hdf5r::H5File$new(paste0("raw_data/", dataset_id, "/cpdata.h5"), "r")
        plate_barcode <- dataset[['meta/experiment/assay_plate_barcode']][]
        data.frame(plate_barcode = plate_barcode)
    })
dataset_ids %>% readr::write_tsv("raw_data/dataset_ids.tsv")
