

library(tidyverse)
library(arrow)

cell_features_TS2PL1 <- arrow::read_parquet(
  file = "product/covid19cq1_SARS_TS2PL1_Cell_MasterDataTable.parquet")

cell_feature_columns <- readr::read_tsv(
  "product/cell_feature_columns_TS_202008.tsv")

cell_metadata_columns <- readr::read_tsv(
  "product/cell_metadata_columns_TS_202008.tsv")


cell_features_rel <- cell_features_TS2PL1 %>%
  dplyr::filter(row == 1, column == 2)
feature_duplication <- cell_feature_columns %>%
  dplyr::rowwise() %>%
  dplyr::do({
    feature1 <- .$feature
    cat("Feature1:", feature1, "\n")
    cell_feature_columns %>%
      dplyr::rowwise() %>%
      dplyr::do({
        feature2 <- .$feature
        n_equal = sum(cell_features_rel[[feature1]] == cell_features_rel[[feature2]])
        if( (feature1 != feature2) && n_equal > 1000){
          cat("Feature1:", feature1, " Feature2:", feature2, "n_equal: ", n_equal, "/", nrow(cell_features_rel), "\n")
        }
        z <- data.frame(
          feature1 = feature1,
          feature2 = feature2,
          n_equal = n_equal)
      })
  })
