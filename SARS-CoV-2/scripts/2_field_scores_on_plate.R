

#layout
#3
#   0   1   2
#2
#   3   4   5
#1
#   6   7   8
#0    1   2   3
###################

load("intermediate_data/rf_scores_field_10XX.Rdata")
rf_scores_field_10XX <- rf_scores_field_10XX %>%
    dplyr::mutate(
        xmin = mod(Image_Metadata_Field, 3),
        ymin = 2-floor(Image_Metadata_Field/3),               
        xmax = 1+mod(Image_Metadata_Field, 3),
        ymax = 3-floor(Image_Metadata_Field/3))

rf_scores_field_10XX %>%
    plyr::d_ply("Plate_Name", function(rf_scores){
        cat("Making rf_score field map plot for plate ", rf_scores$Plate_Name[1], "\n", sep="")
        plot <- ggplot2::ggplot(data=rf_scores) +
            ggplot2::theme_bw() +
            ggplot2::theme(legend.position="bottom") +
            ggplot2::geom_rect(
                mapping=ggplot2::aes(
                    xmin=xmin,
                    xmax=xmax,
                    ymin=ymin,
                    ymax=ymax,
                    fill=1-infectivity_probpos_field)) +
            ggplot2::facet_grid(row ~ column) +
        ggplot2::scale_fill_gradient(
            "Infectivity score (yellow is more infected)",
            low="#00274C",       # MPStats::umich_colors['blue'],
            high="#FF00FF",      # magenta
            expand = c(0, 0)) +
            ggplot2::ggtitle(
                label="Mean Infectivity Score by Field",
                subtitle=paste0(plate_map$Plate_Name[1]))
        ggplot2::ggsave(
            paste0("product/figures/field_scores/plate_",rf_scores$Plate_Name[1],"_rf_scores_200510.pdf"),
            width=13, height=13)
})

#################






#################

load("intermediate_data/rf_scores_field_20XX.Rdata")
rf_scores_field_20XX <- rf_scores_field_20XX %>%
    dplyr::mutate(
        xmin = mod(Field_ID, 3),
        ymin = 2-floor(Field_ID/3),               
        xmax = 1+mod(Field_ID, 3),
        ymax = 3-floor(Field_ID/3))

rf_scores_field_20XX %>%
    plyr::d_ply("Plate_Name", function(rf_scores){
        cat("Making rf_score field map plot for plate ", rf_scores$Plate_Name[1], "\n", sep="")
        plot <- ggplot2::ggplot(data=rf_scores) +
            ggplot2::theme_bw() +
            ggplot2::theme(legend.position="bottom") +
            ggplot2::geom_rect(
                mapping=ggplot2::aes(
                    xmin=xmin,
                    xmax=xmax,
                    ymin=ymin,
                    ymax=ymax,
                    fill=1-infectivity_probpos_field)) +
            ggplot2::facet_grid(row ~ column) +
        ggplot2::scale_fill_gradient(
            "Infectivity score (yellow is more infected)",
            low="#00274C",       # MPStats::umich_colors['blue'],
            high="#FF00FF",      # magenta
            expand = c(0, 0)) +
            ggplot2::ggtitle(
                label="Mean Infectivity Score by Field",
                subtitle=paste0(rf_scores$Plate_Name[1]))
        ggplot2::ggsave(
            paste0("product/figures/field_scores/plate_",rf_scores$Plate_Name[1],"_rf_scores_200510.pdf"),
            width=13, height=13)
})
