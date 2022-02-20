
library(plyr)
library(tidyverse)


decoys <- readr::read_delim(
    file = "../docking/databases/project_decoys_20210323/decoys_final/decoys.smi",
    delim = " ",
    col_names = c("substance_smiles", "substance_zinc_id"))
decoys %>%
    readr::write_tsv("intermediate_data/project_decoys_20210323.tsv")
