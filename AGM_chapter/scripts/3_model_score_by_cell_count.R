library(plyr)
library(tidyverse)
library(ggplot2)
library(MPStats)

cat("Modeling score by cell count\n")


load("intermediate_data/well_scores.Rdata")

plot <- MPStats::plot_score_by_cell_count(
  well_scores=well_scores,
  subtitle="Human MCF7 cells -- compound-profiling experiment (BBBC021v1)")
ggplot2::ggsave(
  plot=plot,
  filename=paste0("product/score_by_cell_count_", MPStats::date_code(), ".pdf"),
  width=8, height=5)
ggplot2::ggsave(
  plot=plot,
  filename=paste0("product/score_by_cell_count_", MPStats::date_code(), ".png"),
  width=8, height=5)
