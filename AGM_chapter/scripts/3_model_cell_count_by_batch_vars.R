library(plyr)
library(tidyverse)
library(MPStats)
library(ggplot2)

cat("Modeling cell count by batch variables\n")

load("intermediate_data/well_scores.Rdata")

model <- MPStats::model_cell_count_by_batch_vars_lm(well_scores=well_scores)
summary(model) %>%
  capture.output(file=paste0("product/cell_count_by_batch_vars_", MPStats::date_code(), ".txt"))


plot <- MPStats::plot_cell_count_by_batch_vars_scatter(
  well_scores=well_scores,
  subtitle="Human MCF7 cells -- compound-profiling experiment (BBBC021v1)")
ggplot2::ggsave(
  plot=plot,
  filename=paste0("product/cell_count_by_batch_vars_scatter_", MPStats::date_code(), ".pdf"),
  width=8, height=5)
ggplot2::ggsave(
  plot=plot,
  filename=paste0("product/cell_count_by_batch_vars_scatter_", MPStats::date_code(), ".png"),
  width=8, height=5)


plot <- MPStats::plot_cell_count_by_batch_vars_density(
  well_scores=well_scores,
  subtitle="Human MCF7 cells -- compound-profiling experiment (BBBC021v1)")
ggplot2::ggsave(
  plot=plot,
  filename=paste0("product/cell_count_by_batch_vars_density_", MPStats::date_code(), ".pdf"),
  width=8, height=5)
ggplot2::ggsave(
  plot=plot,
  filename=paste0("product/cell_count_by_batch_vars_density_", MPStats::date_code(), ".png"),
  width=8, height=5)


