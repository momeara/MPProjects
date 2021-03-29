

library(Zr)


# I cancled this after collecting 275k substances
chembl25_substances <- Zr::catalog_items(
    catalog_short_name = "chembl25",
    output_fields = c("zinc_id", "supplier_code", "substance.preferred_name", "substance.smiles"),
    result_batch_size = 5000,
    temp_file_base = "/scratch/maom_root/maom99/maom/Zr",
    verbose = TRUE)

command <- paste0("cd /scratch/maom_root/maom99/maom && awk '(NR == 1) || (FNR > 1)' Zr_*.csv | sed 's/,/\t/g' | head -n 2000 > chembl25_substances.tsv")
cat(command, "\n", sep = "")
system(command)


