
library(plyr)
library(tidyverse)

source("parameters.R")

date_code <- "20210501"


## APDP Fingerprints ##

cmd <- paste0(
parameters$featurize_substances_program, " \\
  --library_path raw_data/AID2239.sdf.gz \\
  --substance_id_field 'PUBCHEM_SUBSTANCE_ID' \\
  --fingerprint_type 'APDP' \\
  --fingerprint_n_bits 5000 \\
  --output_path intermediate_data/AID2239_APDP_", date_code, " \\
  --verbose
")
cat(cmd, sep = "")
system(cmd)


cmd <- paste0(
parameters$featurize_substances_program, " \\
  --library_path raw_data/UPCMLD_1.sdf \\
  --substance_id_field 'CdId' \\
  --fingerprint_type 'APDP' \\
  --fingerprint_n_bits 5000 \\
  --output_path intermediate_data/UPCMLD_1_APDP_", date_code, " \\
  --verbose
")
cat(cmd, sep = "")
system(cmd)


cmd <- paste0(
parameters$featurize_substances_program, " \\
  --library_path raw_data/UPCMLD_2.sdf \\
  --substance_id_field 'CdId' \\
  --fingerprint_type 'APDP' \\
  --fingerprint_n_bits 5000 \\
  --output_path intermediate_data/UPCMLD_2_APDP_", date_code, " \\
  --verbose
")
cat(cmd, sep = "")
system(cmd)

