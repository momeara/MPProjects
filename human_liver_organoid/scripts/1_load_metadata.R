
library(plyr)
library(tidyverse)

library(googlesheets4)

source("parameters.R")
source("scripts/get_dataset.R")


# gather the dataset summary and ids from google drive
# HLOs/datasets
datasets <- googlesheets4::read_sheet(
    ss = parameters$datasets_googlesheet_id,
    sheet = parameters$datasets_googlesheet_sheet)

datasets %>%
    readr::write_tsv(parameters$datasets_fname)


# Gather marker genes from google drive
# HLOs/datasets
marker_genes <- googlesheets4::read_sheet(
    ss = parameters$marker_genes_googlesheet_id,
    sheet = parameters$marker_genes_googlesheet_sheet) %>%
    dplyr::rename(
        gene_set = `Gene Set`,
        gene = Gene,
        tissue_type = `Tissue Type`)

marker_genes %>%
    readr::write_tsv(parameters$marker_genes_fname)


