
library(MPStats)

cat("Reading in Zprime by plate data\n")

Zprime_by_plate <- MPStats::read_Zprime_by_plate("raw_data/Master55.csv")
save(Zprime_by_plate, file="intermediate_data/Zprime_by_plate.Rdata")
