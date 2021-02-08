
library(plyr)
library(tidyverse)
library(spatstat)
library(testthat)

populate_point_process <- function(
    data,
    coordinate_columns,
    mark_columns = NULL,
    window = NULL) {

    testthat::test_that(
        "coordinate_columns must be of length 2",
        coordinate_columns %>% length == 2)
    testthat::test_that(
        "coordinate_columns exist in the input data",
        all(coordinate_columns %in% names(data)))
    
    if (!is.null(mark_columns)) {
        testthat::test_that(
            "If specified, mark_columns must exist in the input data",
            all(coordinate_columns %in% names(data)))
        marks <- data[mark_columns]
    } else {
        marks <- NULL
    }
    if (is.null(window)) {
        window <- spatstat::owin(
            xrange = data[[coordinate_columns[1]]] %>% range,
            yrange = data[[coordinate_columns[2]]] %>% range)
    }
    
    spatstat::ppp(
        x = data[[coordinate_columns[1]]],
        y = data[[coordinate_columns[2]]],
        window = window,
        marks = marks)
}

clustering_density <- function(
    point_process,
    correction="none") {
    l_function <- spatstat::Lest(
        X = point_process,
        correction = correction)
    mean(l_function$un)
}

