xo
library(plyr)
library(tidyverse)
library(ggplot2)
library(readxl)
library(MPStats)
library(arrow)

### CX5
infectivity_score_CX5_features <- c(
    "Cells_Intensity_IntegratedIntensityEdge_Virus",
    "Cells_Intensity_MeanIntensityEdge_Virus",
    "Cells_Intensity_MaxIntensityEdge_Virus",
    "Cells_Intensity_MaxIntensity_Virus")
add_infectivity_score_CX5 <- function(cell_features){
    cell_features %>%
        dplyr::mutate(
            infectivity_score = -5.064328 +
                Cells_Intensity_IntegratedIntensityEdge_Virus * 1.487025e-01 +
                Cells_Intensity_MeanIntensityEdge_Virus * -3.840196e+01 +
                Cells_Intensity_MaxIntensityEdge_Virus * 4.270269e+01 +
                Cells_Intensity_MaxIntensity_Virus * 4.254849e+01)
}

### CQ1
infectivity_score_CQ1_features <- c(
    "Cells_Intensity_IntegratedIntensityEdge_Virus",
    "Cells_Intensity_MeanIntensityEdge_Virus",
    "Cells_Intensity_MaxIntensityEdge_Virus",
    "Cytoplasm_Intensity_StdIntensity_Virus")

add_infectivity_score_CQ1 <- function(cell_features){
    cell_features %>%
        dplyr::mutate(
            infectivity_score = -0.59661 * Cells_Intensity_IntegratedIntensityEdge_Virus +
                3985.99907 * Cells_Intensity_MeanIntensityEdge_Virus +
                -54.14521 * Cells_Intensity_MaxIntensityEdge_Virus +
                -96.29369 * Cytoplasm_Intensity_StdIntensity_Virus)
}

#### 1999B ####
infectivity_score_field <- arrow::read_parquet(
    file="product/covid19cq1_SARS_1999B_Cell_MasterDataTable.parquet",
    col_select=c(
        "plate_id",
        "Condition",
        "row",
        "column",
        "Image_Metadata_FieldID",
        "Remdesivir_Concentration",
        "Lactoferrin_Concentration",
        infectivity_score_CQ1_features)) %>%
    dplyr::mutate(
        Remdesivir_Concentration = signif(Remdesivir_Concentration, 3),
        Lactoferin_Concentration = signif(Lactoferrin_Concentration, 3)) %>%
    add_infectivity_score_CQ1() %>%
    dplyr::mutate(
        infectivity_score_nc_mean =
            ifelse(Condition %in% c("NC"), infectivity_score, NA) %>% mean(na.rm=TRUE),
        infectivity_score_pc_mean =
            ifelse(Condition %in% c("PC"), infectivity_score, NA) %>% mean(na.rm=TRUE),        
        infectivity_score_cell =
            (infectivity_score-infectivity_score_pc_mean) /
            (infectivity_score_nc_mean - infectivity_score_pc_mean)) %>%
    dplyr::group_by(
        plate_id,
        Condition,
        row,
        column,
        Image_Metadata_FieldID,
        Lactoferrin_Concentration,
        Remdesivir_Concentration) %>%
    dplyr::summarize(
    	infectivity_score_field = mean(infectivity_score_cell)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(plate_id, Remdesivir_Concentration, Lactoferrin_Concentration) %>%
    dplyr::arrange(infectivity_score_field) %>%
    dplyr::mutate(
        infectivity_score_field_rank = dplyr::row_number() / dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
        keep_field = infectivity_score_field_rank >= .66) %>%
    dplyr::left_join(
        plate_map_1999B %>%
        dplyr::mutate(
            Remdesivir_Concentration = signif(Remdesivir_Concentration, 3),
            Lactoferin_Concentration = signif(Lactoferrin_Concentration, 3)) %>%
        dplyr::distinct(
            Remdesivir_Concentration,
            Lactoferrin_Concentration,
            rem_label,
            lf_label),
        by=c("Remdesivir_Concentration", "Lactoferrin_Concentration"))

infectivity_score_field <- infectivity_score_field %>%
    dplyr::filter(Condition == "Treatment") %>%
    dplyr::group_by(Remdesivir_Concentration, Lactoferrin_Concentration) %>%
    dplyr::do({
        fit <- lm(
            formula=infectivity_score_field ~ infectivity_score_field_rank,
            data= dplyr::filter(., infectivity_score_field_rank >= 1/2))
        pred_values <- predict(fit, newdata=.)
        dplyr::mutate(., zscore = (infectivity_score_field - pred_values)/pred_values)
    }) %>%
    dplyr::ungroup() %>% 
    dplyr::mutate(
        outlier_field = infectivity_score_field_rank >=.66 & zscore > 0.6)

infectivity_score_treatment <- infectivity_score_field %>%
    dplyr::filter(infectivity_score_field_rank >=.66, !outlier_field) %>%
    dplyr::group_by(Remdesivir_Concentration, rem_label, Lactoferrin_Concentration, lf_label) %>%
    dplyr::summarize(
        infectivity_score_treatment_mean = mean(infectivity_score_field),
        infectivity_score_treatment_sem = sd(infectivity_score_field)/sqrt(dplyr::n())) %>%
    dplyr::ungroup()


# variants
infectivity_score_treatment_no_outliers <- infectivity_score_field %>%
    dplyr::filter(infectivity_score_field_rank >=.66) %>%
    dplyr::group_by(Remdesivir_Concentration, rem_label, Lactoferrin_Concentration, lf_label) %>%
    dplyr::summarize(
        infectivity_score_treatment_mean = mean(infectivity_score_field),
        infectivity_score_treatment_sem = sd(infectivity_score_field)/sqrt(dplyr::n())) %>%
    dplyr::ungroup()

infectivity_score_treatment_half <- infectivity_score_field %>%
    dplyr::filter(infectivity_score_field_rank >=.5) %>%
    dplyr::group_by(Remdesivir_Concentration, rem_label, Lactoferrin_Concentration, lf_label) %>%
    dplyr::summarize(
        infectivity_score_treatment_mean = mean(infectivity_score_field),
        infectivity_score_treatment_sem = sd(infectivity_score_field)/sqrt(dplyr::n())) %>%
    dplyr::ungroup()

infectivity_score_treatment_auc <- infectivity_score_field %>%
    dplyr::group_by(Remdesivir_Concentration, rem_label, Lactoferrin_Concentration, lf_label) %>%
    dplyr::summarize(
        infectivity_score_treatment_AUC =
            DescTools::AUC(
                x=infectivity_score_field_rank,
                y=infectivity_score_field)) %>%
    dplyr::ungroup()

infectivity_score_treatment_mean <- infectivity_score_field %>%
    dplyr::group_by(Remdesivir_Concentration, rem_label, Lactoferrin_Concentration, lf_label) %>%
    dplyr::summarize(
        infectivity_score_treatment_mean = mean(infectivity_score_field)) %>%
    dplyr::ungroup()





p <- ggplot2::ggplot(data=infectivity_score_field)+
    ggplot2::theme_bw() +
    ggplot2::geom_rect(
        xmin=-Inf,
        xmax=.66,
        ymin=-Inf,
        ymax=Inf,
        fill="grey",
        alpha=.1) +
    ggplot2::geom_point(
        data=infectivity_score_field,
        mapping=ggplot2::aes(
            x=infectivity_score_field_rank,
            y=infectivity_score_field,
            color=outlier_field),
        size=.6) +
    ggplot2::geom_smooth(
        data=infectivity_score_field %>% dplyr::filter(
            infectivity_score_field_rank >= 1/2),
        mapping=ggplot2::aes(
            x=infectivity_score_field_rank,
            y=infectivity_score_field),
        method=lm,
        color="blue",
        se=FALSE,
        size=.2)+
    ggplot2::geom_segment(
        data=infectivity_score_treatment %>% dplyr::mutate(x=.66, xend=1),
        mapping=ggplot2::aes(
            x=x,
            xend=xend,
            y=infectivity_score_treatment_mean,
            yend=infectivity_score_treatment_mean),
        size=.5) +
#    ggplot2::geom_segment(
#        data=infectivity_score_treatment_half %>% dplyr::mutate(x=.55, xend=1),
#        mapping=ggplot2::aes(
#            x=x,
#            xend=xend,
#            y=infectivity_score_treatment_mean,
#            yend=infectivity_score_treatment_mean),
#        size=.5,
#        color="darkgrey") +
#    ggplot2::geom_segment(
#        data=infectivity_score_treatment_auc %>% dplyr::mutate(x=0, xend=1),
#        mapping=ggplot2::aes(
#            x=x,
#            xend=xend,
#            y=infectivity_score_treatment_AUC,
#            yend=infectivity_score_treatment_AUC),
#        size=.5,
#        color="purple") +
#    ggplot2::geom_segment(
#        data=infectivity_score_treatment_mean %>% dplyr::mutate(x=0, xend=1),
#        mapping=ggplot2::aes(
#            x=x,
#            xend=xend,
#            y=infectivity_score_treatment_mean,
#            yend=infectivity_score_treatment_mean),
#        size=.5,
#        color="orange") +
    ggplot2::facet_grid(rem_label ~ lf_label) +
    ggplot2::scale_x_continuous(
        "Rank of mean field infectivity score",
        breaks=c(.5, 2/3, 1),
        labels=c(".5", "2/3", "1")) +
    ggplot2::scale_y_continuous(
        "Infectivity Score") +
    ggplot2::scale_color_discrete("Is outlier")

ggplot2::ggsave(
    filename="product/figures/batch_effects/infectivity_score_1999B_200517.pdf",
    width=13, height=7)


sigmoid_fits <- infectivity_score_field %>%
    dplyr::filter(Lactoferrin_Concentration > 0) %>%
    dplyr::filter(infectivity_score_field_rank >=.66, !outlier_field) %>%
    plyr::ddply(c("Remdesivir_Concentration", "rem_label"), function(df){
    df <- df %>%
        dplyr::mutate(log_dose = log10(Lactoferrin_Concentration+.2))
    tryCatch({
        fit <- drc::drm(
            formula=infectivity_score_field ~ log_dose,
            data=df,
            fct=drc::L.4(fixed=c(NA, NA, NA, NA)))
            upperl=c(35,Inf,Inf,Inf))
        log_dose <- seq(log10(0+.2), log10(25+.2), length.out=200)
        pred_value <- predict(fit, expand.grid(log_dose, 1))
        fit_data <- data.frame(log_dose, pred_value) %>%
          dplyr::mutate(
v            slope=fit$coefficients[1],
            bottom=fit$coefficient[4],
            top=fit$coefficients[2],
            ic50=fit$coefficients[3],
            Lactoferrin_Concentration = 10^(log_dose)-.2)
    }, error=function(e){
        cat("    Failed to fit curve...\n")
        data.frame()
    })
})

sigmoid_fits_auc <- infectivity_score_treatment_auc %>%
    dplyr::filter(Lactoferrin_Concentration > 0) %>%
    plyr::ddply(c("Remdesivir_Concentration", "rem_label"), function(df){
    df <- df %>%
        dplyr::mutate(log_dose = log10(Lactoferrin_Concentration+.2))
    tryCatch({
        fit <- drc::drm(
            formula=infectivity_score_treatment_AUC ~ log_dose,
            data=df,
            fct=drc::L.4(fixed=c(NA, NA, NA, NA)))
            upperl=c(20,Inf,Inf,Inf))
        log_dose <- seq(log10(0+.2), log10(25+.2), length.out=200)
        pred_value <- predict(fit, expand.grid(log_dose, 1))
        fit_data <- data.frame(log_dose, pred_value) %>%
          dplyr::mutate(
            slope=fit$coefficients[1],
            bottom=fit$coefficient[4],
            top=fit$coefficients[2],
            ic50=fit$coefficients[3],
            Lactoferrin_Concentration = 10^(log_dose)-.2)
    }, error=function(e){
        cat("    Failed to fit curve...\n")
        data.frame()
    })
})



p <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position="bottom") +
    ggplot2::geom_errorbar(
        data=infectivity_score_treatment,
        mapping=ggplot2::aes(
            x=log10(Lactoferrin_Concentration+.8),
            ymin=infectivity_score_treatment_mean-infectivity_score_treatment_sem,
            ymax=infectivity_score_treatment_mean+infectivity_score_treatment_sem))+
    ggplot2::geom_line(
        data=infectivity_score_treatment %>% dplyr::filter(Lactoferrin_Concentration > 0),
        mapping=ggplot2::aes(
            x=log10(Lactoferrin_Concentration+.8),
            y=infectivity_score_treatment_mean,
            group=Remdesivir_Concentration),
        size=.7) +
    ggplot2::geom_line(
        data=sigmoid_fits,
        mapping=ggplot2::aes(
            x=log10(Lactoferrin_Concentration+.8),
            y=pred_value,
            group=Remdesivir_Concentration),
        size=1.3) +
    ggplot2::geom_point(
        data=infectivity_score_treatment,
        mapping=ggplot2::aes(
            x=log10(Lactoferrin_Concentration+.8),
            y=infectivity_score_treatment_mean))+
#    ggplot2::geom_point(
#        data=infectivity_score_treatment_no_outliers,
#        mapping=ggplot2::aes(
#            x=log10(Lactoferrin_Concentration+.8),
#            y=infectivity_score_treatment_mean),
#        color="grey",
#        size=.4)+
#    ggplot2::geom_point(
#        data=infectivity_score_treatment_auc,
#        mapping=ggplot2::aes(
#            x=log10(Lactoferrin_Concentration+.8),
#            y=infectivity_score_treatment_AUC),
#        color="purple",
#        size=.6)+
#    ggplot2::geom_line(
#        data=sigmoid_fits_auc,
#        mapping=ggplot2::aes(
#            x=log10(Lactoferrin_Concentration+.8),
#            y=pred_value,
#            group=Remdesivir_Concentration),
#        size=.9,
#        color="purple") +
    ggplot2::scale_x_continuous(
        "Lactoferrin Dose (µg/mL)",
        breaks=log10(c(0, 0.39, 0.78, 1.56, 3.12, 6.25, 12.5, 25)+.8),
        labels=c("0", "0.39", "0.78", "1.56", "3.12", "6.25", "12.5", "25")) +
    ggplot2::facet_wrap(~rem_label, ncol=4) +
    ggplot2::scale_y_continuous(
        "Infectivity Score")   
#    ggplot2::scale_color_continuous(
#        "Remdesivir Concentration (nM)")

ggplot2::ggsave(
    filename="product/figures/dose_response/plate_1999B/infectivity_score_sigmoid_200517.pdf",
    width=8,
    height=4)

