

library(plyr)
library(tidyverse)

source("parameters.R")

date_code <- "20210323"
tag <- paste0("project_substances_", date_code)

## ECFP4 Fingerprints ##
cmd <- paste0(
parameters$featurize_substances_program, " \\
  --library_path intermediate_data/substances_sanitized_", date_code, ".tsv \\
  --substance_id_field 'substance_name' \\
  --smiles_field 'substance_smiles_rdkit' \\
  --output_path intermediate_data/project_substances_", date_code, " \\
  --verbose
")
cat(cmd, sep = "")
system(cmd)


cmd <- paste0(
parameters$featurize_substances_program, " \\
  --library_path intermediate_data/project_decoys_20210323.tsv \\
  --substance_id_field 'substance_zinc_id' \\
  --smiles_field 'substance_smiles' \\
  --output_path intermediate_data/project_decoys_20210323 \\
  --verbose
")
cat(cmd, sep = "")
system(cmd)

cmd <- paste0(
parameters$featurize_substances_program, " \\
  --library_path /scratch/maom_root/maom99/maom/chembl25_substances.tsv \\
  --substance_id_field 'zinc_id' \\
  --smiles_field 'substance.smiles' \\
  --output_path intermediate_data/chembl25_substances_20210323 \\
  --verbose
")
cat(cmd, "\n", sep = "")
system(cmd)

cmd <- paste0(
parameters$featurize_substances_program, " \\
  --library_path raw_data/chembl27_substances_20210323.tsv \\
  --substance_id_field 'zinc_id' \\
  --smiles_field 'smiles' \\
  --output_path intermediate_data/chembl27_substances_ecfp4_20210323 \\
  --verbose
")
cat(cmd, "\n", sep = "")
system(cmd)


tag <- paste0("project_substances_", parameters$date_code)
expand.grid(
    a = c(1),
    b = c(1.8),
    n_neighbors = c(100),
    metric = c("cosine")) %>%
    dplyr::rowwise() %>%
    dplyr::do({
        params <- .
#        alt_tag <- paste0("project_substances_decoys_a=", params$a, ",b=", params$b, ",n_neighbors=", params$n_neighbors, ",metric=", params$metric, "_", date_code)
        alt_tag <- paste0("project_substances_decoys_chembl25_275k_a=", params$a, ",b=", params$b, ",n_neighbors=", params$n_neighbors, ",metric=", params$metric, "_", date_code)
        cat("Computing embedding for: ", alt_tag, "\n", sep = "")
        command <- paste0(
            parameters$embed_umap_program, " ",
            "--dataset ",
            "intermediate_data/project_substances_20210323/fingerprints.parquet ",
            "intermediate_data/project_decoys_20210323/fingerprints.parquet ",
            "intermediate_data/chembl25_substances_20210323/fingerprints.parquet ",
            "--pca_n_components 100 ",
            "--feature_columns intermediate_data/", tag, "/fingerprint_feature_columns.tsv ",
            "--tag ", alt_tag, " ",
            "--umap_a ", params$a, " ",
            "--umap_b ", params$b, " ",
            "--umap_n_neighbors ", params$n_neighbors, " ",
            "--umap_metric ", params$metric, " ",
            "--umap_n_epochs 2000 ",
            "--verbose",
        sep = "")
        cat(command, sep = "")
        system(command)
        data.frame()
    })




## APDP Fingerprints ##
cmd <- paste0(
parameters$featurize_substances_program, " \\
  --library_path intermediate_data/substances_sanitized_", date_code, ".tsv \\
  --substance_id_field 'substance_name' \\
  --smiles_field 'substance_smiles_rdkit' \\
  --fingerprint_type 'APDP' \\
  --fingerprint_n_bits 5000 \\
  --output_path intermediate_data/project_substances_APDP_", date_code, " \\
  --verbose
")
cat(cmd, sep = "")
system(cmd)


cmd <- paste0(
parameters$featurize_substances_program, " \\
  --library_path intermediate_data/project_decoys_20210323.tsv \\
  --substance_id_field 'substance_zinc_id' \\
  --smiles_field 'substance_smiles' \\
  --fingerprint_type 'APDP' \\
  --fingerprint_n_bits 5000 \\
  --output_path intermediate_data/project_decoys_APDP_", date_code, " \\
  --verbose
")
cat(cmd, sep = "")
system(cmd)


cmd <- paste0(
parameters$featurize_substances_program, " \\
  --library_path /scratch/maom_root/maom99/maom/chembl25_substances.tsv \\
  --substance_id_field 'zinc_id' \\
  --smiles_field 'substance.smiles' \\
  --fingerprint_type 'APDP' \\
  --fingerprint_n_bits 5000 \\
  --output_path intermediate_data/chembl25_substances_APDP_20210323 \\
  --verbose
")
cat(cmd, sep = "")
system(cmd)


tag <- paste0("project_substances_APDP_", parameters$date_code)
expand.grid(
    a = c(1),
    b = c(1.8),
    metric = c("cosine")) %>%
    dplyr::rowwise() %>%
    dplyr::do({
        params <- .
        alt_tag <- paste0(
            "project_substances_decoys_APDP_",
            "a=", params$a, ",",
            "b=", params$b, ",",
            "metric=", params$metric, "_",
            date_code)
        cat("Computing embedding for: ", alt_tag, "\n", sep = "")
        command <- paste0(
            parameters$embed_umap_program, " ",
            "--dataset ",
            "intermediate_data/project_substances_APDP_20210323/fingerprints.parquet ",
            "intermediate_data/project_decoys_APDP_20210323/fingerprints.parquet ",
#            "intermediate_data/chembl25_substances_APDP_20210323/fingerprints.parquet ",
            "--feature_columns intermediate_data/", tag, "/fingerprint_feature_columns.tsv ",
            "--tag ", alt_tag, " ",
            "--umap_a ", params$a, " ",
            "--umap_b ", params$b, " ",
            "--umap_metric ", params$metric, " ",
            "--verbose",
        sep = "")
        cat(command, sep = "")
        system(command)
        data.frame()
    })
