

library(Zr)
library(arrow)

# I cancled this after collecting 275k substances
chembl25_substances <- Zr::catalog_items(
    catalog_short_name = "chembl25",
    output_fields = c("zinc_id", "supplier_code", "substance.preferred_name", "substance.smiles"),
    result_batch_size = 5000,
    temp_file_base = "/scratch/maom_root/maom99/maom/Zr",
    verbose = TRUE)

command <- paste0("cd /scratch/maom_root/maom99/maom && awk '(NR == 1) || (FNR > 1)' Zr_*.csv | sed 's/,/\t/g' > chembl25_substances.tsv")
cat(command, "\n", sep = "")
system(command)

chembl25_fingerprints <- arrow::read_parquet("intermediate_data/chembl25_substances_20210323/fingerprints.parquet")

chembl25_fingerprints %>%
    dplyr::sample_n(5000) %>%
    arrow::write_parquet("intermediate_data/chembl25_substances_20210323/fingerprints_5k.parquet")

chembl25_fingerprints %>%
    dplyr::sample_n(10000) %>%
    arrow::write_parquet("intermediate_data/chembl25_substances_20210323/fingerprints_10k.parquet")

chembl25_fingerprints %>%
    dplyr::sample_n(50000) %>%
    arrow::write_parquet("intermediate_data/chembl25_substances_20210323/fingerprints_50k.parquet")


load("raw_data/chembl27_substances_20210323.Rdata")

full_chembl_data %>%
    readr::write_tsv("raw_data/chembl27_substances_20210323.tsv")
