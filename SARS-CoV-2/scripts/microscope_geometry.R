
get_image_width_height <- function() {
    # get the image height/width for CQ1
    source("scripts/database.R")
    con <- get_primary_database_connection("covid19cq1")
    plate_id <- "2019A"
    image_table_name <- paste0("SARS_", plate_id, "_Per_Image")
    image_features <- con %>%
        dplyr::tbl(image_table_name) %>%
        dplyr::collect(n = Inf)
    image_width <- image_features$Image_Width_Virus[1]
    image_height <- image_features$Image_Height_Virus[1]
    
    assertthat::assert_that(image_width == 2560)
    assertthat::assert_that(image_height == 2160)
    c(image_width, image_height)
}



cq1_pixel_to_mm <- function(coordinate){
    # in Jiji -> image/Properties, pixel width = 0.0104167
    # There are 25.4 mm per inch
    coordinate * 0.0104167 * 25.4
}

# for the CQ1 there are 9 fields layed out like this:
# 0 1 2
# 3 4 5
# 6 7 8

cq1_field_pixel_to_well_mm_x <- function(
    coordinate_x,
    field_index,
    image_width = 2560) {
    field_x <- field_index %% 3
    (coordinate_x + image_width * field_x) * 0.0104167 * 25.4
}

cq1_field_pixel_to_well_mm_y <- function(
    coordinate_y,
    field_index,
    image_height = 2160) {
    field_y <- floor(field_index / 3)
    (coordinate_y + image_height * field_y) * 0.0104167 * 25.4
}
