

library(plyr)
library(tidyverse)

source("parameters.R")

date_code <- "20210323"
tag <- paste0("project_substances_", date_code)

system(paste0(
parameters$featurize_substances_program, " \\
  --library_path ../data_curation/raw_data/substances_", date_code, ".tsv \\
  --substance_id_field 'substance_name' \\
  --smiles_field 'substance_smiles' \\
  --output_path intermediate_data/", tag, " \\
  --verbose
"))



expand.grid(
    a = c(1),
    b = c(1)) %>%
    dplyr::rowwise() %>%
    dplyr::do({
        params <- .
        alt_tag <- paste0("project_substances_a=", params$a, ",b=", params$b, "_", date_code)
        cat("Computing embedding for: ", alt_tag, "\n", sep = "")
        command <- paste0("
            ", parameters$embed_umap_program, " \\
              --dataset intermediate_data/", tag, "/fingerprints.parquet \\
              --feature_columns intermediate_data/", tag, "/fingerprint_feature_columns.tsv \\
              --tag ", alt_tag, " \\
              --umap_a ", params$a, " \\
              --umap_b ", params$b, " \\
              --verbose
        ", sep = "")
        cat(command, sep = "")
        system(command)
        data.frame()
    })
