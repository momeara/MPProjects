








library(plyr)
library(dplyr)
library(monocle3)
library(MPStats)


load("intermediate_data/cds_CZ_x.Rdata")

# by inspecting the distribution of marker genes across the UMAP clusters
# decide that the clusters 1, 4, 7, 9 and 11 are likely enriched for hepatocyte cells
cds_CZ_x_hepatocyte <- cds_CZ_x[,
    which(
        monocle3::clusters(cds_CZ_x, reduction_method = "UMAP") %in% c(1, 4, 7, 9, 11),
        arr.ind = TRUE)]
save(cds_CZ_x_hepatocyte, file = "intermediate_data/cds_CZ_x_hepatocyte.Rdata")

col_data <- cds_CZ_x_hepatocyte %>%
    SummarizedExperiment::colData() %>%
    data.frame()

col_data %>% dplyr::count(sample)
