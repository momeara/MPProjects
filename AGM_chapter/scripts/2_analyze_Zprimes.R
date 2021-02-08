library(plyr)
library(tidyverse)
library(MPStats)
library(ggplot2)
library(arrow)
library(patchwork)

cat("Analyzing Zprime scores by features and dyes\n")

load("intermediate_data/Zprime_by_plate.R")
cell_features <- arrow::read_parquet(file="intermediate_data/cell_features.parquet")

dye_plot <- MPStats::plot_Zprime_by_dye(Zprime_by_plate)
ggplot2::ggsave(
  plot=dye_plot,
  filename=paste0("product/Zprime_by_dye_", MPStats::date_code(), ".pdf"),
  width=4, height=6,
  useDingbats=FALSE)
ggplot2::ggsave(
  plot=dye_plot,
  filename=paste0("product/Zprime_by_dye_", MPStats::date_code(), ".png"),
  width=4, height=6)

split_violin_plot <- MPStats::plot_top_features_by_control(cell_features)
ggplot2::ggsave(
  plot=split_violin_plot,
  filename=paste0("product/top_features_by_control_", MPStats::date_code(), ".pdf"),
  width=5, height=7,
  useDingbats=FALSE)
ggplot2::ggsave(
  plot=split_violin_plot,
  filename=paste0("product/top_features_by_control_", MPStats::date_code(), ".png"),
  width=5, height=7)


Zprime_plot <- MPStats::plot_Zprime_by_top_feature(cell_features)
ggplot2::ggsave(
  plot=Zprime_plot,
  filename=paste0("product/Zprime_by_top_feature_", MPStats::date_code(), ".pdf"),
  width=5, height=7,
  useDingbats=FALSE)
ggplot2::ggsave(
  plot=Zprime_plot,
  filename=paste0("product/Zprime_by_top_feature_", MPStats::date_code(), ".png"),
  width=5, height=7)


# Combine the split_violin and Zprime plots
top_plot <-
  split_violin_plot +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      legend.position="bottom") +
  Zprime_plot +
    ggplot2::theme(             
      axis.title.y=ggplot2::element_blank(),
      axis.text.y=ggplot2::element_blank(),                 
      axis.ticks.y=ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      legend.position="bottom")

ggplot2::ggsave(
  plot=top_plot,
  filename=paste0("product/top_features_", MPStats::date_code(), ".pdf"),
  width=7, height=7,
  useDingbats=FALSE)
ggplot2::ggsave(
  plot=top_plot,
  filename=paste0("product/top_features_", MPStats::date_code(), ".png"),
  width=7, height=7)


