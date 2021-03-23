#TODO

## what are size factors?

run when creating a cell_data_se

  cds <- estimate_size_factors(cds)
            method = c("mean-geometric-mean-total", "mean-geometric-mean-log-total")

  estimate_sf_sparse
     cell_total <- Matrix::colSums(counts)
     sfs <- cell_total / exp(mean(log(cell_total)))

cds <- estimateDispersions(cds)    
when combining datasets, what is going on with 'conf' conflicts?


## Filtering out bad cells



## Count Normalization

Compute pearson residuals for each gene
Jan Lause, Philipp Berens, and Dmitry Kobak
Analytic Pearson residuals for normalization of single-cell RNA-seq UMI data

