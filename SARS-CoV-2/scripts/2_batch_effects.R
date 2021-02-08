
library(plyr)
library(tidyverse)
library(ggplot2)
library(readxl)
library(MPStats)

plate_ids <- readr::read_tsv("raw_data/plate_ids.tsv")
load("intermediate_data/image_scores.Rdata")
feature_columns <- readr::read_tsv("raw_data/cell_feature_columns.tsv", col_names="feature")

plate_ids %>%
    plyr::a_ply(1, function(plate){
        plate_id <- plate$plate_id[1]
        cat("Plotting batch effects for plate ", plate_id, "\n", sep="")

        cat("  Getting cell features ...\n")
        cell_features <- arrow::read_parquet(
               file=paste0("product/SARS_", plate_id, "_Cell_MasterDataTable.parquet"),
               cols=c(
                   feature_columns$feature,
                   row, column, is_control, Compound, dose_nM, plate_id, Image_Metadata_Field))
        
        cat("  Computing well summaries ...\n")
        well_features <- feature_columns$feature %>%
            plyr::ldply(function(feature){
                cat("Computing well summaries for feature '", feature, "'\n", sep="")
                cell_features %>%
                    dplyr::select(plate_id, dose_nM, row, column, is_control, Compound, value=!!feature) %>%
                    dplyr::mutate(normed_value = scale(value)) %>%
                    dplyr::group_by(plate_id, dose_nM, row, column, is_control, Compound) %>%
                    dplyr::summarize(
                         well_mean = mean(normed_value),
                         well_sd = sd(normed_value)) %>%
                    dplyr::mutate(feature=!!feature)
            })
        
        well_features <- well_features %>%
            dplyr::mutate(feature_name = feature) %>%
            tidyr::separate(col=feature, into=c("object", "measure_type", "measure", "channel"), sep="_")
        
        ###########################################
        cat("  Plotting AreaShape row batch effects ...\n")        
        plot <- ggplot2::ggplot() +
          ggplot2::theme_bw() +
          theme(legend.position="bottom") +
          ggplot2::geom_smooth(
            data=well_features %>% dplyr::filter(measure_type == "AreaShape"),
            mapping=ggplot2::aes(x=row, y=well_mean, color=object, group=feature_name),
            size=1,
            alpha=1,
            se=FALSE) +
          facet_wrap(~measure) +
          ggplot2::scale_x_continuous("Plate Row") +
          ggplot2::scale_y_continuous("Average well feature value (Standard Deviations)") +
          ggplot2::ggtitle("AreaShape row batch bffects", subtitle=paste0("Plate id: ", plate_id))
        ggplot2::ggsave(
          filename=paste0("product/figures/batch_effects/", plate_id, "_AreaShape_row_batch_effects_200430.pdf"),
          plot=plot,
          width=8,
          height=8)
        
        cat("  Plotting Intensity row batch effects ...\n")                
        plot <- ggplot2::ggplot() +
          ggplot2::theme_bw() +
          theme(legend.position="bottom") +
          ggplot2::geom_smooth(
            data=well_features %>% dplyr::filter(measure_type == "Intensity"),           
            mapping=ggplot2::aes(x=row, y=well_mean, color=object, group=feature_name),
            size=1,
            alpha=1,
            se=FALSE) +
          facet_grid(measure~channel) +
          ggplot2::scale_x_continuous("Plate Row") +
          ggplot2::scale_y_continuous("Average well feature value (Standard Deviations)") +
          ggplot2::ggtitle("Intensity row batch effects", subtitle=paste0("Plate id: ", plate_id))
        ggplot2::ggsave(
          filename=paste0("product/figures/batch_effects/", plate_id, "_Intensity_row_batch_effects_200430.pdf"),
          plot=plot,
          width=8,
          height=20)
        
        cat("  Plotting RadialDistribution row batch effects ...\n")                        
        plot <- ggplot2::ggplot() +
          ggplot2::theme_bw() +
          theme(legend.position="bottom") +
          ggplot2::geom_smooth(
            data=well_features %>% dplyr::filter(measure_type == "RadialDistribution", !is_control),    
            mapping=ggplot2::aes(x=row, y=well_mean, color=object, group=feature_name),
            size=1,
            alpha=1,
            se=FALSE) +
          facet_grid(measure~channel) +
          ggplot2::scale_x_continuous("Plate Row") +
          ggplot2::scale_y_continuous("Average well feature value (Standard Deviations)") +
          ggplot2::ggtitle("RadialDistribution row batch bffects", subtitle=paste0("Plate id: ", plate_id))    
        ggplot2::ggsave(
          filename=paste0("product/figures/batch_effects/", plate_id, "_RadialDistribution_row_batch_effects_200430.pdf"),
          plot=plot,
          width=6,
          height=8)
        
    })

                            
