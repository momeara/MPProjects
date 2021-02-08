#' Load a feature table from CellProfiler HDF5 data export
#'
#' Usage:
#'   h5_dataset <- hdf5r::H5File$new("DefaultOUT.h5", "r")
#'   h5_dataset$ls(recursive = TRUE)
#'   # look up the table names
#'   cell_features <- h5_dataset %>%
#'      load_features_h5(
#'        table_name = "Measurements/2020-10-15-15-10-50/Nuclei",
#'        exclude_features = c("ObjectNumber"))
#'
#' Watch out that some of the data associated with a table may not
#' have the same number of values. In that case either exclude a few
#' features or just include what you need with the include_features or
#' exclude_features arguments.
#'
#' Returns a tibble::tibble with columns as features, rows as objects
load_features_h5 <- function(
    h5_dataset,
    table_name,
    exclude_features = NULL,
    include_features = NULL,
    verbose = TRUE) {
    environment <- h5_dataset[[table_name]]
    if (!is.null(include_features)) {
        features <- include_features
        if (!is.null(exclude_features)) {
            cat("ERROR: Please use as most one of 'include_features' and 'exclude_features'\n")
            return(NULL)
        }
    } else if (!is.null(exclude_features)) {
        if (!is.null(include_features)) {
            cat("ERROR: Please use as most one of 'include_features' and 'exclude_features'\n")
            return(NULL)
        }        
        features <- environment %>%
            names() %>%
            tibble::tibble(feature = .) %>%
            dplyr::filter(
                !(feature %in% exclude_features)) %>%
            magrittr::extract2("feature")
    } else {
        features <- environment %>%
            names()
    }
    
    if (verbose) {
        cat("Extracting '", length(features), "' for table '", table_name, "'\n", sep = "")
    }
    data <- purrr::map_dfc(features, function(feature) {
            values <- environment[[feature]][["data"]][]
            if (verbose) {
                cat("getting ", length(values), " values for feature '", feature, "' for table '", table_name, "'\n", sep = "")
            }
            z <- tibble::tibble(value =  values)
            names(z)[1] <- feature
            z
        })
}
