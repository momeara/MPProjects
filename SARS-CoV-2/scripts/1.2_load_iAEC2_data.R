
library(plyr)
library(tidyverse)
library(RMySQL)
library(magrittr)
library(tictoc)
library(arrow)

# Charles Zhang, May 8th, 2020
# iAEC2 cells
# infected with MOIs of 0, 10, and 25
# For each MOI, they were treated with
#    vehicle control
#    remdesivir
#    lactoferrin
#    rem+LF
# Cells were stained by
#    hoechst
#    acetylated tubulin
#    nucleocapsid antibody
#    a population of cells expresses tomato
cell_features <- readr::read_csv(
    file = "~/bucket/Reid_FileTransfer/Charles_Files/features_iAE2C_MOI_0_10_25_rem_lf_pruned_200508.csv")

cell_features %>% arrow::write_parquet(
    sink = "product/features_iAE2C_MOI_0_10_25_rem_lf_pruned_200508.parquet")
