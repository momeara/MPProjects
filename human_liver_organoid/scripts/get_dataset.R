
library(plyr)
library(tidyverse)
library(monocle3)


# Load the output of the 10x genomics CellRanger pipeline stored in turbo
# https://cole-trapnell-lab.github.io/monocle3/docs/starting/#10x-output

get_dataset <- function(sample_id, dataset_path = NULL) {
    if (is.null(dataset_path)) {
        source("parameters.R")
        datasets <- readr::read_tsv(parameters$datasets_fname, col_types = readr::cols())
        dataset_path <- paste(
            parameters$data_base_dir,
            datasets %>%
                dplyr::filter(sample_id == !!sample_id) %>%
                dplyr::pull(data_path),
            sep = "/")
    } 
    cat("Loading cell ranger datset into Monocle3 with path '", dataset_path, "'\n", sep = "")
    monocle3::load_cellranger_data(dataset_path)
}
