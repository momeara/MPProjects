library(tidyverse)
library(MPStats)

cat("Reading in well scores\n")

# load well scores
well_scores <- MPStats::read_well_scores(input="raw_data/ProbPos_CellCount_55Plates.csv")
save(well_scores, file="intermediate_data/well_scores.Rdata")

# show plate layout for plate 4 on week 9
well_scores %>%
  MPStats::plate_layout(week=9, plate=4) %>%
  readr::write_tsv(paste("product/example_plate_map_", MPStats::date_code(), ".tsv"))
