
library(plyr)
library(dplyr)
library(monocle3)
library(MPStats)


source("parameters.R")
marker_genes <- readr::read_tsv(
    parameters$marker_genes_fname)


MPStats::make_scatterboard

