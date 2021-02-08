



library(plyr)
library(tidyverse)
library(ggplot2)
library(readxl)
library(MPStats)

plate_ids <- readr::read_tsv("raw_data/plate_ids.tsv")
load("intermediate_data/image_scores.Rdata")
feature_columns <- readr::read_tsv("raw_data/cell_feature_columns.tsv")

plate_id <- "1999B"
cell_features <- arrow::read_parquet(
       file=paste0("product/SARS_", plate_id, "_Cell_MasterDataTable.parquet"),
       cols=c(
           feature_columns$feature,
           row, column, is_control, Compound, dose_nM, plate_id, Image_Metadata_Field))


feature_columns <- feature_columns %>%
    dplyr::mutate(
       transform = dplyr::case_when(
           # these should be log-transformed, but they have values <= 0 so use log1p
           feature == "Cells_Intensity_MinIntensityEdge_Hoe" ~ "log1p",
           feature == "Cells_Intensity_MinIntensity_Hoe" ~ "log1p",
           feature == "Nuclei_Intensity_MADIntensity_Hoe" ~ "log1p",
           feature == "Nuclei_Intensity_MinIntensityEdge_Hoe" ~ "log1p",
           feature == "Nuclei_Intensity_MinIntensity_Hoe" ~ "log1p",
           feature == "Cytoplasm_Intensity_LowerQuartileIntensity_Hoe" ~ "log1p",
           feature == "Cytoplasm_Intensity_MinIntensityEdge_Hoe" ~ "log1p",
           feature == "Cytoplasm_Intensity_MinIntensity_Hoe" ~ "log1p",
           # log transform
           feature %>% stringr::str_detect("_AreaShape_Area") ~ "log",
           feature %>% stringr::str_detect("_AreaShape_Compactness") ~ "log",
           feature %>% stringr::str_detect("_AreaShape_MajorAxisLength") ~ "log",
           feature %>% stringr::str_detect("_AreaShape_MaxFeretDiameter") ~ "log",
           feature %>% stringr::str_detect("_AreaShape_MaximumRadius") ~ "log",
           feature %>% stringr::str_detect("_AreaShape_MeanRadius") ~ "log",
           feature %>% stringr::str_detect("_AreaShape_MedianRadius") ~ "log",
           feature %>% stringr::str_detect("_AreaShape_MinFeretDiameter") ~ "log",
           feature %>% stringr::str_detect("_AreaShape_MinorAxisLength") ~ "log",
           feature %>% stringr::str_detect("_AreaShape_Perimeter") ~ "log",
           feature %>% stringr::str_detect("_Intensity_") ~ "log",
           feature %>% stringr::str_detect("_RadialDistribution_") ~ "log",
           TRUE ~ "identity"))
feature_columns %>% readr::write_tsv("raw_data/cell_feature_columns.tsv")

plate_feature_densities <- feature_columns %>%
    plyr::adply(1, function(feature_column){
        cat("Computing plate density distributions for feature '", feature_column$feature[1], "'\n", sep="")
        values <- cell_features[[feature_column$feature[1] ]]
        if(feature_column$transform[1] == "log1p"){
            cat("transforming feature '", feature_column$feature[1], "' with 'log1p'\n", sep="")
            if(sum(values+1 <= 0) > 0){
                cat("   ", sum(values <= 0), " values are <= 0\n", sep="")
            }
            values <- log(values+1)
        } else if(feature_column$transform[1] == "log"){
            cat("transforming feature '", feature_column$feature[1], "' with 'log'\n", sep="")
            if(sum(values <= 0) > 0){
                cat("   ", sum(values <= 0), " values are <= 0\n", sep="")
            }
            values <- log(values)
        }
        dens <- density(values)
        tibble::tibble(
            plate_id = plate_id,
            feature=feature_column$feature[1],
            value = dens$x,
            density = dens$y)
    })
plate_feature_densities <- plate_feature_densities %>%
    dplyr::mutate(feature_name = feature) %>%
    tidyr::separate(
        col = feature,
        into = c("object", "measure_type", "measure", "channel"),
        sep = "_")

###########################################
cat("  Plotting AreaShape feature densities ...\n")
plot <- ggplot2::ggplot() +
  ggplot2::theme_bw() +
  theme(legend.position = "bottom") +
  ggplot2::geom_line(
    data = plate_feature_densities %>% dplyr::filter(measure_type == "AreaShape"),
    mapping = ggplot2::aes(x = value, y = log10(density + 1), color = object, group = feature_name),
    size = 1,
    alpha = 1) +
  facet_wrap(~measure, scale = "free") +
  ggplot2::scale_x_continuous("Feature Value") +
  ggplot2::scale_y_continuous("Log Density") +
  ggplot2::ggtitle("AreaShape densities ", subtitle = paste0("Plate id: ", plate_id))
ggplot2::ggsave(
  filename = paste0("product/figures/feature_densities/", plate_id, "_AreaShape_feature_densities_200501.pdf"),
  plot = plot,
  width = 8,
  height = 8)


cat("  Plotting Intensity feature densities ...\n")
plot <- ggplot2::ggplot() +
  ggplot2::theme_bw() +
  theme(legend.position="bottom") +
  ggplot2::geom_line(
    data=plate_feature_densities %>% dplyr::filter(measure_type == "Intensity"),           
    mapping=ggplot2::aes(x=value, y=log10(density+1), color=object, group=feature_name),
    size=1,
    alpha=1) +
  facet_wrap(measure~channel, scale="free") +
  ggplot2::scale_x_continuous("Feature Value") +
  ggplot2::scale_y_continuous("Log Density") +
  ggplot2::ggtitle("Intensity feature distributions", subtitle=paste0("Plate id: ", plate_id))
ggplot2::ggsave(
  filename=paste0("product/figures/feature_densities/", plate_id, "_Intensity_feature_densities_200501.pdf"),
  plot=plot,
  width=15,
  height=20)

cat("  Plotting RadialDistribution row batch effects ...\n")                        
plot <- ggplot2::ggplot() +
  ggplot2::theme_bw() +
  theme(legend.position="bottom") +
  ggplot2::geom_line(
    data=plate_feature_densities %>% dplyr::filter(measure_type == "RadialDistribution"),    
    mapping=ggplot2::aes(x=value, y=log10(density+1), color=object, group=feature_name),
    size=1,
    alpha=1) +
  facet_wrap(measure~channel, scale="free") +
  ggplot2::scale_x_continuous("Feature Value") +
  ggplot2::scale_y_continuous("Log Density") +
  ggplot2::ggtitle("RadialDistribution feature distributions", subtitle=paste0("Plate id: ", plate_id))    
ggplot2::ggsave(
  filename=paste0("product/figures/feature_densities/", plate_id, "_RadialDistribution_feature_dendisities_200501.pdf"),
  plot=plot,
  width=6,
  height=8)
