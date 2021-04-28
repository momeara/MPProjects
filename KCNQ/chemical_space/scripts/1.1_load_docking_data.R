# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(tidyverse)

source("../docking/scripts/gather_docking_features.R")

source("parameters.R")

date_code <- "20210323"
tag <- paste0("project_substances_", date_code)



docking_features <- load_mol2_features(
    mol2_fname = "raw_data/KCNQ2_7CR2_retigabine_AB_20201028,project_20210323,,20210402/poses.mol2",
    verbose = TRUE)

docking_features %>%
    readr::write_tsv(
        file = "intermediate_data/KCNQ2_7CR2_retigabine_AB_20201028,project_20210323,,20210402/docking_features.tsv")


# generate fingerprints from docked substances
cmd <- paste0(
parameters$featurize_substances_program, " \\
  --library_type mol2 \\
  --library_path raw_data/KCNQ2_7CR2_retigabine_AB_20201028,project_20210323,,20210402/poses.mol2 \\
  --library_fields _Name \\
  --output_path intermediate_data/KCNQ2_7CR2_retigabine_AB_20201028,project_20210323,,20210402 \\
  --fingerprint_type ECFP4 \\
  --fingerprint_n_bits 4096 \\
  --verbose
")
cat(cmd, sep = "")
system(cmd)
