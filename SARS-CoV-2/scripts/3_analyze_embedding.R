
library(plyr)
library(tidyverse)
library(arrow)
library(ggplot2)
library(MPStats)
library(caret)
library(gbm)

library(doParallel)
cl <- makePSOCKcluster(30)
registerDoParallel(cl)

feature_columns <- readr::read_tsv("raw_data/cell_feature_columns.tsv")

master_plate_id <- 1001
plate_id <- "1001005A"

sample_indices <- readr::read_tsv(
    file=paste0("~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/SARS_", master_plate_id, "_Cell_umap2_2M_15_0.0/sample_indices.tsv"),
    col_names="cell_index") %>%
    dplyr::mutate(cell_index = cell_index+1)

cell_features <- dplyr::bind_rows(
    arrow::read_parquet(
        file=paste0("product/SARS_", master_plate_id, "0050A_Cell_MasterDataTable.parquet"),
        col_select=feature_columns=feature_columns$feature) %>%
        dplyr::mutate(dose_nM = 50),
    arrow::read_parquet(
        file=paste0("product/SARS_", master_plate_id, "0250A_Cell_MasterDataTable.parquet"),
        col_select=feature_columns=feature_columns$feature) %>%
        dplyr::mutate(dose_nM = 250),
    arrow::read_parquet(
        file=paste0("product/SARS_", master_plate_id, "0500A_Cell_MasterDataTable.parquet"),
        col_select=feature_columns=feature_columns$feature) %>%
        dplyr::mutate(dose_nM = 500),
    arrow::read_parquet(
        file=paste0("product/SARS_", master_plate_id, "1000A_Cell_MasterDataTable.parquet"),
        col_select=feature_columns=feature_columns$feature) %>%
        dplyr::mutate(dose_nM = 1000),
    arrow::read_parquet(
        file=paste0("product/SARS_", master_plate_id, "2000A_Cell_MasterDataTable.parquet"),
        col_select=feature_columns=feature_columns$feature) %>%
        dplyr::mutate(dose_nM = 2000)) %>%
    dplyr::mutate(cell_index = dplyr::row_number())

umap_embedding <- arrow::read_parquet(
    file=paste0("~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/SARS_", master_plate_id, "_Cell_umap2_2M_15_0.0/umap_embedding.parquet"))

cluster_labels <- arrow::read_parquet(
    file=paste0("~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/SARS_", master_plate_id, "_Cell_umap2_2M_15_0.0/hdbscan_clustering_min100.parquet"))

umap_features <- sample_indices %>%
    dplyr::left_join(cell_features,  by=c("cell_index")) %>%
    dplyr::bind_cols(
       umap_embedding,
       cluster_labels)

umap_features <- umap_features %>%
    dplyr::mutate(
       umap_infected = dplyr::case_when(
           master_plate_id == 1001 & UMAP_2 < -3 ~ TRUE,
           master_plate_id == 1002 & UMAP_1 < -2.5 & UMAP_2 > 14 ~ TRUE,
           master_plate_id == 1002 & UMAP_1 > 9.4 & UMAP_2 > 9 ~ TRUE,
           master_plate_id == 1003 & UMAP_2 < -3 ~ TRUE,
           TRUE ~ FALSE))
                                  

data_matrix <- umap_features %>%
    dplyr::mutate(label = dplyr::case_when(
        cluster_label == 1 ~ "infected",
        TRUE ~ "not_infected") %>% factor(
            levels=c("infected", "not_infected"), labels=c("infected", "not_infected"))) %>%
    dplyr::select(label, tidyselect::one_of(feature_columns$feature)) %>%
    tidyr::drop_na()

#ok
#data_matrix %>% caret::nearZeroVar()

# check for correlated features
#cor_matrix <- data_matrix %>% dplyr::select(-roi) %>% as.matrix() %>% cor

in_training <- caret::createDataPartition(data_matrix$label, p = .005, list = FALSE)
data_train <- data_matrix %>% dplyr::slice(in_training)
data_test  <- data_matrix %>% dplyr::slice(-in_training)

fitControl <- caret::trainControl(## 10-fold CV
    method = "repeatedcv",
    number = 2,
    classProbs=TRUE)

plsFit2 <- caret::train(label ~ ., data = data_train, 
                 method = "pls", 
                 trControl = fitControl,
                 preProc = c("center", "scale"),
                 tuneLength = 3,
                 verbose=TRUE)
predFit2 <- predict(plsFit2, newdata = data_test)

gbmFit1 <- caret::train(label ~ ., data = data_train, 
                 method = "gbm", 
                 trControl = fitControl,
                 preProc = c("center", "scale"),                 
                 tuneLength=5,
                 verbose=TRUE)

gbmFit1 %>% caret::varImp()

twoClassSummary(data_set, lev = levels(test_set$obs))
