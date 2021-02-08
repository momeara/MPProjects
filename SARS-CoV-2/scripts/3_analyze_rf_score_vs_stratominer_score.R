








load("intermediate_data/rf_scores_field_20XX.Rdata")
rf_well_scores <- rf_scores_field_20XX %>%
    dplyr::mutate(Plate_Name = paste0("SARS_", Plate_Name)) %>%
    dplyr::group_by(Plate_Name, Compound, row, column) %>%
    dplyr::summarize(rf_probpos_well = mean(infectivity_probpos_field))

stratominer_well_scores <- arrow::read_parquet(
    file="intermediate_data/stratominer_well_scores_2016A-2019A_200514.parquet") %>%
    dplyr::select(
        Plate_Name, Compound, row, column, stratominer_probpos_well=probPOSITIVE)    

well_scores <- dplyr::inner_join(
    rf_well_scores,
    stratominer_well_scores,
    by=c("Plate_Name", "Compound", "row", "column"))


ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position="bottom") + 
    ggplot2::geom_point(
        data=well_scores %>%
            dplyr::filter(!(Compound %in% c("PC", "NC"))),
        mapping=ggplot2::aes(                                      
            x=rf_probpos_well,
            y=stratominer_probpos_well),
        color="black") +
    ggplot2::geom_point(
        data=well_scores %>% dplyr::filter(Compound == "NC"),
        mapping=ggplot2::aes(                                      
            x=rf_probpos_well,
            y=stratominer_probpos_well),
        color="magenta") +
    ggplot2::geom_point(
        data=well_scores %>% dplyr::filter(Compound == "PC"),
        mapping=ggplot2::aes(                                      
            x=rf_probpos_well,
            y=stratominer_probpos_well),
        color="blue") +
    ggplot2::geom_smooth(
        data=well_scores,
        mapping=ggplot2::aes(
            x=rf_probpos_well,
            y=stratominer_probpos_well),
        color="red") +
    ggplot2::scale_x_continuous("mean RF score per well") +
    ggplot2::scale_y_continuous("mean Stratominer score per well") +
    ggplot2::ggtitle(
        label="Compare RF vs Stratominer scores",
        subtitle="plates: 2016A-2019A (CQ1) (Blue: PC, Mageneta: NC)")

ggplot2::ggsave(
   file="product/figures/batch_effects/rf_vs_stratominer_score_2016A-2019A.pdf")
