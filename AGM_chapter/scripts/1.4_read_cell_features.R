
library(MPStats)

cat("Reading in cell feature data\n")

cell_features <- MPStats::read_cell_features("raw_data/AGM_FeatureSet_RAW55.csv")
save(cell_features, file="intermediate_data/cell_features.Rdata")
arrow::write_parquet("intermediate_data/cell_features.parquet")
