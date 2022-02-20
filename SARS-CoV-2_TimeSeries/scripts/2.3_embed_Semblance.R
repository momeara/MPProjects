


library(Semplance)

cell_feature_columns_TS_202008 <- readr::read_tsv(
    "product/cell_feature_columns_TS_202008.tsv")

cell_features <- arrow::read_parquet(
    file = "intermediate_data/covid19cq1_SARS_TS2PL1_Cell_MasterDataTable.parquet",
    col_select = cell_feature_columns_TS_202008$feature)



for (sample_size in c(100, 1000, 10000, 100000, 1000000)) {
    cat("Computing Semblance for ", sample_size, " features\n", sep = "")
    tictoc::tic()
    s <- cell_features %>%
        dplyr::sample_n(!!sample_size) %>%
        as.matrix() %>%
        Semblance::ranksem()
    tictoc::toc()
    x <- sort( sapply(ls(),function(x){object.size(get(x))}))
    x <- x[names(x) == "s"]
    cat("sample_size ", sample_size, ", bytes: ", x, "\n", sep ="")
}
