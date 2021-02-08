library(plyr)
library(tidyverse)
library(ggplot2)
library(readxl)
library(MPStats)

plate_map <- readxl::read_excel("raw_data/FDA_Platemap_FULL_DrugAnnotation_Version3.xlsx")

###################################
### Viral intensity well scores ###
###################################
load("intermediate_data/image_scores_CX5_100X.Rdata")
meta_well_scores <- image_scores_CX5_100X %>%
    dplyr::distinct(
        master_plate_id,
        Compound,
        dose_nM,
        row,
        column)

viral_intensity_well_scores <- image_scores_CX5_100X %>%
    dplyr::select(
        master_plate_id,
        Compound,
        dose_nM,
        Image_Count_Cells,
        Image_Classify_Positive_PctObjectsPerBin) %>%
    dplyr::group_by(master_plate_id, dose_nM) %>%
        dplyr::mutate(
            cell_count_baseline = mean(Image_Count_Cells),
            prob_pose_baseline = mean(Image_Classify_Positive_PctObjectsPerBin[which(Compound == "Negative Control")])) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(master_plate_id, dose_nM, Compound) %>%
        dplyr::mutate(
            normed_cell_count = Image_Count_Cells / cell_count_baseline,
            normed_prob_pos = Image_Classify_Positive_PctObjectsPerBin / prob_pose_baseline) %>%
        dplyr::summarize(
            mean_normed_prob_pos = mean(normed_prob_pos),
            sem_normed_prob_pos = sd(normed_prob_pos)/sqrt(dplyr::n()),
            mean_normed_cell_count = mean(normed_cell_count),
            sem_normed_cell_count = sd(normed_cell_count)/sqrt(dplyr::n())) %>%
    dplyr::ungroup()

#viral_intensity_sigmoid_fits <- viral_intensity_scores %>%
#    dplyr::filter(keep_field) %>%
#    dplyr::filter(Compound != "PC", Compound != "NC") %>%
#    plyr::ddply(c("Compound"), function(compound_scores) {
#    compound_scores <- compound_scores %>%
#        dplyr::mutate(log_dose = log10(dose_nM) - 9)
#    tryCatch({
#        fit <- drc::drm(
#            formula = infectivity_score_field ~ log_dose,
#            data = compound_scores,
#            fct = drc::L.4(fixed = c(NA, NA, NA, NA)))
#        log_dose <- seq(log10(50) - 9, log10(2000) - 9, length.out = 200)
#        pred_value <- predict(fit, expand.grid(log_dose, 1))
#        fit_data <- data.frame(log_dose, pred_value) %>%
#          dplyr::mutate(
#            slope = fit$coefficients[1],
#            bottom = 0,
#            top = fit$coefficients[2],
#            ic50 = fit$coefficients[3],
#            dose_nM = 10^(log_dose + 9))
#    }, error = function(e) {
#        cat("    Failed to fit curve for compound '", compound_scores$Compound[1], "'\n", sep = "")
#        data.frame()
#    })
#})



#######################
## infectivity score ##
#######################
infectivity_score_treatment <- arrow::read_parquet(
    file = "product/infectivity_score_treatment_10XX_200519.parquet") %>%
    dplyr::mutate(
        master_plate_id = plate_id %>%
            stringr::str_extract("^....") %>%
            as.numeric()) %>%
    dplyr::filter(!(plate_id %>% stringr::str_detect("1005....C")))

infectivity_score_field <- arrow::read_parquet(
    file = "product/infectivity_score_field_10XX_200519.parquet") %>%
    dplyr::mutate(
        master_plate_id = plate_id %>% stringr::str_extract("^....") %>% as.numeric())



infectivity_score_sigmoid_fits <- infectivity_score_field %>%
    dplyr::filter(keep_field) %>%
    dplyr::filter(Compound != "PC", Compound != "NC") %>%
    plyr::ddply(c("Compound"), function(compound_scores) {
    compound_scores <- compound_scores %>%
        dplyr::mutate(log_dose = log10(dose_nM) - 9)
    tryCatch({
        fit <- drc::drm(
            formula = infectivity_score_field ~ log_dose,
            data = compound_scores,
            fct = drc::L.4(fixed = c(NA, NA, NA, NA)))
        log_dose <- seq(log10(50) - 9, log10(2000) - 9, length.out = 200)
        pred_value <- predict(fit, expand.grid(log_dose, 1))
        fit_data <- data.frame(log_dose, pred_value) %>%
          dplyr::mutate(
            slope = fit$coefficients[1],
            bottom = 0,
            top = fit$coefficients[2],
            ic50 = fit$coefficients[3],
            dose_nM = 10^(log_dose + 9))
    }, error = function(e) {
        cat("    Failed to fit curve for compound '", compound_scores$Compound[1], "'\n", sep = "")
        data.frame()
    })
})



################################
### UMAP cluster well scores ###
################################

umap_well_scores <- dplyr::bind_cols(
    dplyr::bind_rows(
        arrow::read_parquet(
            file = "product/SARS_10030050A_Cell_MasterDataTable.parquet",
            col_select = c("Image_Metadata_WellID", "Compound")) %>%
            dplyr::mutate(dose_nM = 50),
        arrow::read_parquet(
            file="product/SARS_10030250A_Cell_MasterDataTable.parquet",
            col_select = c("Image_Metadata_WellID", "Compound")) %>%
            dplyr::mutate(dose_nM = 250),
        arrow::read_parquet(
            file="product/SARS_10030500A_Cell_MasterDataTable.parquet",
            col_select = c("Image_Metadata_WellID", "Compound")) %>%
            dplyr::mutate(dose_nM = 500),
        arrow::read_parquet(
            file="product/SARS_10031000A_Cell_MasterDataTable.parquet",
            col_select = c("Image_Metadata_WellID", "Compound")) %>%
            dplyr::mutate(dose_nM = 1000),
        arrow::read_parquet(
            file="product/SARS_10032000A_Cell_MasterDataTable.parquet",
            col_select = c("Image_Metadata_WellID", "Compound")) %>%
            dplyr::mutate(dose_nM = 2000)),
    arrow::read_parquet(
               file = "~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/SARS_1003_Cell_umap2_into2M_15_0.0/umap_embedding.parquet")) %>%
    dplyr::mutate(master_plate_id = 1003) %>%
    dplyr::mutate(is_infected = UMAP_2 < -3.1) %>%
    dplyr::group_by(master_plate_id, dose_nM, Compound) %>%
    dplyr::summarize(
        cell_count = dplyr::n(),
        infected_count = sum(is_infected)) %>%
    dplyr::mutate(
        prob_pos = infected_count / cell_count)

########################################
### Stratominer distance well scores ###
########################################
stratominer_well_scores <- readxl::read_excel("raw_data/Stratominer_distance_hits_200425.xlsx") %>%
    dplyr::mutate(dose_nM = Barcode %>%
       stringr::str_extract("[0-9][0-9][0-9][0-9]A$") %>%
       stringr::str_replace("A", "") %>%
       as.numeric()) %>%
    dplyr::select(
        master_plate_id=FDA384Barcode,
        dose_nM,
        Compound,
        distance=DistanceScore) %>%
    dplyr::filter(Compound != "PC", Compound != "NC")    

stratominer_compound_scores <- stratominer_well_scores %>%
    dplyr::mutate(distance = ifelse(distance < .1, .1, distance)) %>%
    tidyr::pivot_wider(
        names_from=dose_nM,
        names_prefix="dose_",
        values_from=distance) %>%
    dplyr::mutate(
        dose_50 = ifelse(is.na(dose_50), 0.1, dose_50),
        dose_250 = ifelse(is.na(dose_250), dose_50, dose_250),
        dose_500 = ifelse(is.na(dose_500), dose_250, dose_500),
        dose_1000 = ifelse(is.na(dose_1000), dose_500, dose_1000),
        dose_2000 = ifelse(is.na(dose_2000), dose_1000, dose_2000)) %>%
    dplyr::mutate(
        score50_250 = log10(((dose_250-dose_50)/dose_50 + 1)),
        score250_500 = log10(((dose_500-dose_250)/dose_250 + 1)),
        score500_1000 = log10(((dose_1000-dose_500)/dose_500 + 1)),
        score1000_2000 = log10(((dose_2000-dose_1000)/dose_1000 + 1))) %>%
    dplyr::mutate(
        score = score50_250 + score250_500 + score500_1000 + score1000_2000) %>%
    dplyr::arrange(desc(score)) %>%
    dplyr::mutate(
        rank_compound = factor(
            x=paste0(dplyr::row_number(), "_", Compound),
            levels=paste0(dplyr::row_number(), "_", Compound))) %>%
    dplyr::select(
        master_plate_id, rank_compound, Compound,
        dose_50, dose_250, dose_500, dose_1000, dose_2000,
        score50_250, score250_500, score500_1000, score1000_2000,
        score)

stratominer_well_scores <- stratominer_well_scores %>%
    dplyr::left_join(stratominer_compound_scores, by=c("master_plate_id", "Compound")) %>%
    dplyr::select(
        master_plate_id,
        dose_nM,
        Compound,
        stratominer_distance=distance,
        rank_compound,
        stratominer_quality=score)


############################################
### Stratominer distance May 4th version ###
############################################
stratominer_well_scores_200504 <- readxl::read_excel(
    path="raw_data/S25_stratominer_field_level_200504.xlsx",
    sheet='withcontrols') %>%
    dplyr::transmute(
        master_plate_id=FDA384Barcode,
        dose_nM=CONC2,
        Compound,
        stratominer_NPC=NPC/100) %>%
    dplyr::filter(Compound != "PC", Compound != "NC")

##################
## Manual Score ##
##################
source("scripts/classifiers.R")
manual_score <- cell_features %>%
    dplyr::mutate(manual_score = manual_score(cell_features)) %>%
    dplyr::group_by(
         Condition,
         Remdesivir_Concentration,
         Lactoferrin_Concentration,
         row,
         column) %>%
    dplyr::summarize(
         manual_score_well_mean = mean(manual_score)) %>%
    dplyr::ungroup() %>%
# group by plate?
    dplyr::mutate(
          manual_score_plate_mean = mean(manual_score_well_mean)) %>%
    dplyr::group_by(
         master_plate_id,
         dose_nM,
         Compound,
         manual_score_plate_mean) %>%
    dplyr::summarize(
         manual_score_condition_mean = mean(manual_score_well_mean/manual_score_plate_mean),
         manual_score_condition_mean = sd(manual_score_well_mean/manual_score_plate_mean)/sqrt(dplyr::n())) %>%
    dplyr::ungroup()

qHTS_hit_score <- manual_score %>%
    dplyr::left_join(
        manual_score %>%
            tidyr::pivot_wider(
                 names_from=dose_nM,
                 names_prefix="dose_",
                 values_from=manual_score_condition_mean) %>%
            dplyr::mutate(
                 has_trend = qHTS_has_trend(.),
                 significant_inhibition = qHTS_significant_inhibition(threshold=0.5))



###################
### RF Scores V1 ##
###################

rf_scores_v1 <- readxl::read_excel("raw_data/RF_Hits_Plate5.xlsx") %>%
    dplyr::select(
        plate_id = Barcode,
        well_id = FDA384WellID,
        Compound,
        prob_pos = probPOSITIVE) %>%
    dplyr::mutate(
        master_plate_id = plate_id %>%
            stringr::str_extract("SARS_....") %>%
            stringr::str_replace("SARS_", "") %>%
            as.numeric(),
        dose_nM = plate_id %>%
            stringr::str_extract("....B") %>%
            stringr::str_replace("B", "") %>%
            as.numeric(),
        row = well_id %>%
            stringr::str_extract("^[A-Z]") %>%
            purrr::map_int(~which(LETTERS==., arr.ind=T)),
        column = well_id %>%
            stringr::str_extract("[0-9]+$") %>%
            as.integer())

################
### RF Scores ##
################

load("intermediate_data/rf_scores_field_10XX.Rdata")
rf_scores_points <- rf_scores_field_10XX %>%
    dplyr::mutate(rf_score = 1-infectivity_probpos_field)

sigmoid_fits <- rf_scores_points %>%
    dplyr::filter(Compound != "PC", Compound != "NC") %>%
    plyr::ddply(c("Compound"), function(compound_scores){
    #cat("Fitting dose response for compound '", compound_scores$Compound[1], "'...\n", sep="")
    compound_scores <- compound_scores %>%
        dplyr::mutate(log_dose = log10(dose_nM) - 9)
    tryCatch({
        fit <- drc::drm(
            formula=rf_score ~ log_dose,
            data=compound_scores,
            fct=drc::L.4(fixed=c(NA, 0, 1, NA)))
        log_dose <- seq(log10(4)-9, log10(2000)-9, length.out=200)
        pred_value <- predict(fit, expand.grid(log_dose, 1))
        fit_data <- data.frame(log_dose, pred_value) %>%
          dplyr::mutate(
            slope=fit$coefficients[1],
            bottom=0,
            top=fit$coefficients[2],
            ic50=fit$coefficients[3],
            dose_nM = 10^(log_dose + 9))
    }, error=function(e){
        cat("    Failed to fit curve for compound '", compound_scores$Compound[1], "'\n", sep="")
        data.frame()
    })
})



rf_scores <- rf_scores_well_10XX %>%
    dplyr::filter(!(Plate_Name %>% stringr::str_detect("^SARS_1005....B"))) %>%
    dplyr::transmute(
         master_plate_id,
         dose_nM,
         Compound,
         rf_score_median = infectivity_probpos_well_NC_median,
         rf_score_sd = infectivity_probpos_well_PC_sd,
         rf_score_low = pmax(0, rf_score_median-rf_score_sd),
         rf_score_high = pmin(rf_score_median+rf_score_sd, 1))


###################
## Selected Hits ##
###################

load("intermediate_data/hits_10XX.Rdata")



######################
### combine scores ###
######################
well_scores <- meta_well_scores %>%
    dplyr::left_join(
        viral_intensity_well_scores %>%
        dplyr::select(
            master_plate_id,
            Compound,
            dose_nM,
            viral_intensity_prob_pos = mean_normed_prob_pos,
            viral_intensity_prob_pos_sem = sem_normed_prob_pos,
            cell_count = mean_normed_cell_count,
            cell_count_sem = sem_normed_cell_count),
        by=c("master_plate_id", "dose_nM", "Compound")) %>%
    dplyr::left_join(
        infectivity_score_treatment %>%
        dplyr::select(
            master_plate_id,
            Compound,
            dose_nM,
            infectivity_score_treatment_mean,
            infectivity_score_treatment_sem),
        by=c("master_plate_id", "dose_nM", "Compound"))%>%
#    dplyr::left_join(
#        stratominer_well_scores %>%
#            dplyr::select(
#                master_plate_id, dose_nM, Compound,
#                stratominer_distance),
#        by=c("master_plate_id", "dose_nM", "Compound")) %>%
#    dplyr::left_join(
#        stratominer_well_scores_200504 %>%
#            dplyr::select(
#                master_plate_id, dose_nM, Compound,
#                stratominer_NPC),
#        by=c("master_plate_id", "dose_nM", "Compound")) %>%
#    dplyr::left_join(
#        umap_well_scores %>%
#            dplyr::select(master_plate_id, dose_nM, Compound, umap_prob_pos=prob_pos),
#        by=c("master_plate_id", "dose_nM", "Compound")) %>%
#    dplyr::left_join(
#        rf_scores %>%
#            dplyr::select(
#                master_plate_id,
#                dose_nM,
#                Compound,
#                rf_score_median,
#                rf_score_low,
#                rf_score_high),
#        by=c("master_plate_id", "dose_nM", "Compound")) %>%
    dplyr::left_join(
        hits_10XX %>%
            dplyr::select(master_plate_id, Compound) %>%
            dplyr::mutate(hit="Hit"),
        by=c("master_plate_id", "Compound"))

make_dose_response_plot <- function(
    well_scores, subtitle, score_types=c(), fits=NULL) {
    plot <- ggplot2::ggplot(data=well_scores) +
        ggplot2::theme_bw() +
        theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust=.5)) +
        ggplot2::geom_hline(yintercept=1, size=.3)
    if("cell_count" %in% names(score_types)){
        cat("Adding cell count layer with color ", score_types['cell_count'], "\n", sep="")
        plot <- plot +
            ggplot2::geom_errorbar(
               mapping=ggplot2::aes(
                 x=dose_nM,
                 ymin=cell_count-cell_count_sem,
                 ymax=cell_count+cell_count_sem),
               color=score_types['cell_count'],
               width=.1) +
            ggplot2::geom_line(
               mapping=ggplot2::aes(x=dose_nM, y=cell_count),
               color=score_types['cell_count'],
               size=1)
    }
    if("viral_intensity" %in% names(score_types)){
        cat("Adding viral intensity layer with color ", score_types['viral_intensity'], "\n", sep="")
        plot <- plot +
            ggplot2::geom_errorbar(
               mapping=ggplot2::aes(
                 x=dose_nM,
                 ymin=viral_intensity_prob_pos-viral_intensity_prob_pos_sem,
                 ymax=viral_intensity_prob_pos+viral_intensity_prob_pos_sem),
               color=score_types['viral_intensity'],
               width=.1) +
            ggplot2::geom_line(
               mapping=ggplot2::aes(x=dose_nM, y=viral_intensity_prob_pos),
               color=score_types['viral_intensity'],
               size=.3)
    }
    if("infectivity_score" %in% names(score_types)){
        cat("Adding infectivity score layer with color ", score_types['infectivity_score'], "\n", sep="")
        plot <- plot +
            ggplot2::geom_errorbar(
               mapping=ggplot2::aes(
                 x=dose_nM,
                 ymin=infectivity_score_treatment_mean-infectivity_score_treatment_sem,
                 ymax=infectivity_score_treatment_mean+infectivity_score_treatment_sem),
               color=score_types['infectivity_score'],
               width=.1) +
            ggplot2::geom_line(
               mapping=ggplot2::aes(x=dose_nM, y=infectivity_score_treatment_mean),
               color='red', #score_types['infectivity_score'],
               size=.9)
    }
    if("stratominer_distance" %in% names(score_types)){
        cat("Adding stratominer distance layer with color ", score_types['stratominer_distance'], "\n", sep="")        
        plot <- plot +
            ggplot2::geom_line(
               mapping=ggplot2::aes(x=dose_nM, y=stratominer_distance),
               color=score_types['stratominer_distance'],
               size=1)
    }
    if("stratominer_NPC" %in% names(score_types)){
        cat("Adding stratominer NCP layer with color ", score_types['stratominer_NPC'], "\n", sep="")
        plot <- plot +
            ggplot2::geom_line(
               mapping=ggplot2::aes(x=dose_nM, y=stratominer_NPC),
               color=score_types['stratominer_NPC'],
               size=.3)
    }
    if("umap_cluster" %in% names(score_types)){
        cat("Adding umap_cluster layer with color ", score_types['umap_cluster'], "\n", sep="")
        plot <- plot +
            ggplot2::geom_line(
               mapping=ggplot2::aes(x=dose_nM, y=umap_prob_pos),
               color=score_types['umap_cluster'],
               size=1)
    }    
    if("rf_score" %in% names(score_types)){
        cat("Adding RF score layer with color ", score_types['rf_score'], "\n", sep="")        
        plot <- plot +
#            ggplot2::geom_boxplot(
#                data=rf_scores_points %>%
#                    dplyr::filter(!(Compound %in% c("Positive Control", "Negative Control"))),
#                mapping=ggplot2::aes(
#                    x=dose_nM,
#                    y=rf_score,
#                    group=dose_nM),
#                color=score_types['rf_score'],
#                width=.05)
            ggplot2::geom_jitter(
                 data=rf_scores_points %>%
                     dplyr::filter(!(Compound %in% c("Positive Control", "Negative Control"))),
                mapping=ggplot2::aes(
                    x=dose_nM,
                    y=rf_score),
                color=score_types['rf_score'],
                width=0.05, height=0)
    }
    if("hits" %in% names(score_types)){
        cat("Adding hits indicator layer with color ", score_types['hits'], "\n", sep="")        
        plot <- plot +
            MPStats::geom_indicator(
               mapping=ggplot2::aes(indicator=hit),
               color=score_types['hits'],
               xpos="left", ypos="top", group=2, size=3)
    }        
    plot <- plot +
        ggplot2::facet_grid(row~column) +
        MPStats::geom_indicator(
            mapping=ggplot2::aes(indicator=Compound),
            xpos="left", ypos="top", group=1, size=3) +
        ggplot2::ggtitle(
            label="Dose Response",
            subtitle=subtitle) +
        ggplot2::scale_x_log10(
            "Dose nM",
            breaks=c(50, 250, 500, 1000, 2000),
            labels=c("50nM", "250nM", "0.5uM", "1uM", "2uM")) +
        ggplot2::scale_y_continuous(
            paste0(score_types, ": ", names(score_types), collapse="   "))
}

master_plate_ids <- c(1001, 1002, 1003, 1004, 1005)
p <- master_plate_ids %>%
    purrr::map(function(master_plate_id){
        cat("Making dose response plot for master plate '", master_plate_id, "'\n", sep="")
        make_dose_response_plot(
            well_scores = well_scores %>%
                dplyr::filter(master_plate_id==!!master_plate_id) %>%
                dplyr::filter(!(Compound %in% c("Positive Control", "Negative Control"))),
            score_type=c(
                cell_count="black",
                infectivity_score="red",
                #viral_intensity="red",
                #stratominer_NPC="red",                
                #rf_score="blue",
                hits="red"),
            subtitle=paste0("Plate ", master_plate_id))
        ggplot2::ggsave(
            paste0("product/figures/dose_response/", master_plate_id, "_200519.png"),
            width=30, height=25)
})               


# with out fits
p <- make_dose_response_plot(
    well_scores = well_scores %>%
        dplyr::semi_join(hits_10XX, by=c("Compound")),
    score_type=c(
        cell_count="black",
        infectivity_score="red",
        #viral_intensity="red",
        #stratominer_NPC="red",                
        #rf_score="blue",
        hits="red"),
    subtitle=paste0("qHTS Hits")) +
    ggplot2::facet_wrap(~Compound) +
    ggplot2::theme(
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
    ggplot2::coord_cartesian(ylim=c(-.1, 5))    
ggplot2::ggsave(
    paste0("product/figures/dose_response/qHTS_hits_200525.pdf"),
    width=22, height=20)


# with fits
p <- make_dose_response_plot(
    well_scores = well_scores %>%
        dplyr::semi_join(hits_10XX, by=c("Compound")),
    score_type=c(
        cell_count="black",
        infectivity_score="red",
        #viral_intensity="red",
        #stratominer_NPC="red",                
        #rf_score="blue",
        hits="red"),
    subtitle=paste0("qHTS Hits"),
    fits=sigmoid_fits %>%
        dplyr::semi_join(hits_10XX, by=c("Compound")))
p <- p + 
    ggplot2::facet_wrap(~Compound) +
    ggplot2::theme(
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
    ggplot2::geom_line(
        data=sigmoid_fits %>%
            dplyr::semi_join(hits_10XX, by=c("Compound")),
        mapping=ggplot2::aes(x=dose_nM, y=pred_value),
        color='red') +
    ggplot2::coord_cartesian(ylim=c(-.1, 3))
ggplot2::ggsave(
    paste0("product/figures/dose_response/qHTS_hits_with_fits_200525.pdf"),
    width=22, height=20)
