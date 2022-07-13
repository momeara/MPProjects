library(tidyverse)

source("parameters.R")

date_code <- "20220713"





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




# library_substances_APDP
command <- paste0(
    parameters$embed_umap_program, " ",
    "--dataset ",
    "intermediate_data/AID2239_APDP_20210501/fingerprints.parquet ",
    "intermediate_data/UPCMLD_1_APDP_20210501/fingerprints.parquet ",
    "intermediate_data/UPCMLD_2_APDP_20210501/fingerprints.parquet ",
    "--feature_columns intermediate_data/AID2239_APDP_20210501/fingerprint_feature_columns.tsv ",
    "--tag library_substances_APDP ",
    "--umap_a 1 ",
    "--umap_b 1.8 ",
    "--umap_metric cosine ",
    "--verbose",
    sep = "")
cat(command, sep = "")
system(command)



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



date_code <- "20220711"
tag <- paste0("project_substances_decoys_ChemGPT-1.2B")
expand.grid(
    a = c(1),
    b = c(2),
    n_neighbors = c(100),
    metric = c("cosine")) %>%
    dplyr::rowwise() %>%
    dplyr::do({
        params <- .
        alt_tag <- paste0("project_substances_decoys_ChemGPT-1.2B_a=", params$a, ",b=", params$b, ",n_neighbors=", params$n_neighbors, ",metric=", params$metric, "_", date_code)
        cat("Computing embedding for: ", alt_tag, "\n", sep = "")
        command <- paste0(
            parameters$embed_umap_program, " ",
            "--dataset ",
            "intermediate_data/project_substances_ChemGPT-1.2B_20220711/fingerprints.parquet ",
            "intermediate_data/project_decoys_ChemGPT-1.2B_20220711/fingerprints.parquet ",
            "--pca_n_components 100 ",
            "--feature_columns intermediate_data/project_substances_ChemGPT-1.2B_20220711/fingerprint_feature_columns.tsv ",
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

