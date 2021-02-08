library(plyr)
library(tidyverse)
library(fuzzyjoin)
library(ggplot2)
library(readxl)
library(MPStats)


image_scores_20XX <- arrow::read_parquet("product/image_scores_20XX.parquet")
load("intermediate_data/rf_scores_field_20XX.Rdata")

meta_well_scores <- rf_scores_field_20XX %>%
    dplyr::distinct(
        master_plate_id,
        Compound,
        dose_nM) %>%
    dplyr::filter(Compound != "Metformin")

###################################
### Viral intensity well scores ###
###################################

viral_intensity_well_scores <- image_scores_20XX %>%
    dplyr::select(
        master_plate_id,
        Compound,
        dose_nM,       
        Image_Count_Cells,
        Image_Count_Nucleoli,        
        Image_Classify_Positive_PctObjectsPerBin) %>%
    dplyr::group_by(master_plate_id, dose_nM) %>%
        dplyr::mutate(
            cell_count_baseline = mean(Image_Count_Cells),
            nucleoli_count_baseline = mean(Image_Count_Nucleoli),
            prob_pose_baseline = mean(Image_Classify_Positive_PctObjectsPerBin)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(master_plate_id, dose_nM, Compound) %>%
        dplyr::mutate(
            normed_cell_count = Image_Count_Cells / cell_count_baseline,
            normed_nucleoli_count = Image_Count_Nucleoli / nucleoli_count_baseline,
            normed_prob_pos = Image_Classify_Positive_PctObjectsPerBin / prob_pose_baseline) %>%
        dplyr::summarize(
            mean_normed_prob_pos = mean(normed_prob_pos),
            sem_normed_prob_pos = sd(normed_prob_pos)/sqrt(dplyr::n()),
            mean_normed_cell_count = mean(normed_cell_count),
            sem_normed_cell_count = sd(normed_cell_count)/sqrt(dplyr::n()),
            mean_normed_nucleoli_count = mean(normed_nucleoli_count),
            sem_normed_nucleoli_count = sd(normed_nucleoli_count)/sqrt(dplyr::n())) %>%
    dplyr::ungroup()

###################
## UMAP Clusters ##
###################
umap_well_scores <- dplyr::bind_cols(
    arrow::read_parquet(
       file="product/top_hits_Cell_MasterDataTable.parquet",
       col_select=c("master_plate_id", "Compound", "dose_nM")),                
    arrow::read_parquet(
       file="~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/top_hits_umap2_15_0.0/hdbscan_clustering_min100.parquet")) %>%
    dplyr::group_by(master_plate_id, dose_nM, Compound) %>%
    dplyr::summarize(
        well_prob_cluster_9 = mean(cluster_label == 9),               
        well_prob_cluster_2 = mean(cluster_label == 2),
        well_prob_cluster_1 = mean(cluster_label == 1))




################
### RF Scores ##
################


rf_scores <- rf_scores_field_20XX %>%
    dplyr::filter(Compound != "BLANK") %>%
    dplyr::group_by(master_plate_id, Compound, dose_nM) %>%
    dplyr::summarize(
        rf_score_mean = mean(1-infectivity_probpos_field),
        rf_score_sem = sd(infectivity_probpos_field)/sqrt(dplyr::n()),
        rf_score_low = pmax(0, rf_score_mean-rf_score_sem),
        rf_score_high = pmin(rf_score_mean+rf_score_sem, 1)) %>%
    dplyr::ungroup()

rf_scores_points <- rf_scores_field_20XX %>%
    dplyr::filter(Compound != "BLANK", Compound != "PC", Compound != "NC") %>%
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
            fct=drc::L.4(fixed=c(NA, 0, NA, NA)))
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
        cat("    Failed to fit curve...\n")
        data.frame()
    })
})

rf_scores_cell_count <- rf_scores_field_20XX %>%
    dplyr::select(
        master_plate_id,
        Compound,
        dose_nM,       
        cell_count_field) %>%
    dplyr::group_by(master_plate_id, dose_nM) %>%
        dplyr::mutate(cell_count_plate = mean(cell_count_field)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(cell_count_field_normed = cell_count_field / cell_count_plate) %>%
    dplyr::group_by(master_plate_id, dose_nM, Compound) %>%
        dplyr::summarize(
            cell_count_well_mean = mean(cell_count_field_normed),
            cell_count_well_sem = sd(cell_count_field_normed)/sqrt(dplyr::n()),
            cell_count_well_low = cell_count_well_mean - cell_count_well_sem,
            cell_count_well_high = cell_count_well_mean + cell_count_well_sem) %>%
    dplyr::ungroup()




################################
## Stratominer 200514 version ##
################################
stratominer_well_scores <- arrow::read_parquet(
    file="intermediate_data/stratominer_well_scores_2016A-2019A_200514.parquet") %>%
    dplyr::mutate(master_plate_id = Plate_Name %>% stringr::str_extract("20..")) %>%
    dplyr::select(
        master_plate_id, Compound, dose_nM, stratominer_probpos_well=probPOSITIVE)    

####################
## Combine scores ##
####################

well_scores <- meta_well_scores %>%
    dplyr::left_join(
        viral_intensity_well_scores %>%
        dplyr::select(
            master_plate_id,
            Compound,
            dose_nM,
            viral_intensity_prob_pos = mean_normed_prob_pos,
            viral_intensity_prob_pos_sem = sem_normed_prob_pos,
            ##cell_count = mean_normed_cell_count,
            ##cell_count_sem = sem_normed_cell_count,
            nucleoli_count = mean_normed_nucleoli_count,
            nucleoli_count_sem = sem_normed_nucleoli_count),
        by=c("master_plate_id", "dose_nM", "Compound")) %>%
    fuzzyjoin::fuzzy_join(
        rf_scores %>%
            dplyr::select(
                master_plate_id,
                dose_nM,
                Compound,
                rf_score_mean),
        by=c("master_plate_id", "dose_nM", "Compound"),
        match_fun=list(`==`, function(v1, v2){abs(v1-v2)<.01},`==`),
        mode="left") %>%
    dplyr::select(-master_plate_id.y, -dose_nM.y, -Compound.y) %>%
    dplyr::rename(
        master_plate_id = master_plate_id.x,
        dose_nM = dose_nM.x,
        Compound = Compound.x) %>%
    dplyr::left_join(
        rf_scores_cell_count %>%
        dplyr::select(
            master_plate_id,
            Compound,
            dose_nM,
            cell_count_well_mean,
            cell_count_well_low,
            cell_count_well_high),
        by=c("master_plate_id", "dose_nM", "Compound")) %>%
    fuzzyjoin::fuzzy_join(
        stratominer_well_scores %>%
            dplyr::select(
                master_plate_id,
                dose_nM,
                Compound,
                stratominer_probpos_well),
        by=c("master_plate_id", "dose_nM", "Compound"),
        match_fun=list(`==`, function(v1, v2){abs(v1-v2)<.01},`==`),
        mode="left") %>%
    dplyr::select(-master_plate_id.y, -dose_nM.y, -Compound.y) %>%    
    dplyr::rename(
        master_plate_id = master_plate_id.x,
        dose_nM = dose_nM.x,
        Compound = Compound.x) %>%
    dplyr::left_join(
        umap_well_scores %>%
           dplyr::rename(
              umap_prob_cluster_9 = well_prob_cluster_9,                   
              umap_prob_cluster_1 = well_prob_cluster_1,
              umap_prob_cluster_2 = well_prob_cluster_2),      
        by=c("master_plate_id", "dose_nM", "Compound"))       





make_dose_response_plot <- function(
    well_scores,
    subtitle,
    score_types,
    rf_scores_points,
    fits) {
    plot <- ggplot2::ggplot(data=well_scores) +
        ggplot2::theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5)) +
        ggplot2::geom_hline(yintercept=1, size=.3)
    if("cell_count" %in% names(score_types)){
        cat("Adding cell count layer with color ", score_types['cell_count'], "\n", sep="")
        plot <- plot +
            ggplot2::geom_errorbar(
               mapping=ggplot2::aes(
                 x=dose_nM,
                 ymin=cell_count_well_low,
                 ymax=cell_count_well_high),
               color=score_types['cell_count'],
               width=.1) +
            ggplot2::geom_line(
               mapping=ggplot2::aes(x=dose_nM, y=cell_count_well_mean),
               color=score_types['cell_count'],
               size=1)
    }
    if("nucleoli_count" %in% names(score_types)){
        cat("Adding nucleoli count layer with color ", score_types['nucleoli_count'], "\n", sep="")
        plot <- plot +
            ggplot2::geom_errorbar(
               mapping=ggplot2::aes(
                 x=dose_nM,
                 ymin=nucleoli_count-nucleoli_count_sem,
                 ymax=nucleoli_count+nucleoli_count_sem),
               color=score_types['nucleoli_count'],
               width=.1) +
            ggplot2::geom_line(
               mapping=ggplot2::aes(x=dose_nM, y=nucleoli_count),
               color=score_types['nucleoli_count'],
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
               color="red",
               width=.1) +
            ggplot2::geom_line(
               mapping=ggplot2::aes(x=dose_nM, y=viral_intensity_prob_pos),
               color=score_types['viral_intensity'],
               size=1)
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
               size=1)
    }
    if("umap_cluster_1" %in% names(score_types)){
        cat("Adding umap_cluster_1 layer with color ", score_types['umap_cluster_1'], "\n", sep="")
        plot <- plot +
            ggplot2::geom_line(
               mapping=ggplot2::aes(x=dose_nM, y=umap_prob_cluster_1/umap_prob_cluster_9*10),
               color=score_types['umap_cluster_1'],
               size=1)
    }    
    if("umap_cluster_2" %in% names(score_types)){
        cat("Adding umap_cluster_2 layer with color ", score_types['umap_cluster_2'], "\n", sep="")
        plot <- plot +
            ggplot2::geom_line(
               mapping=ggplot2::aes(x=dose_nM, y=umap_prob_cluster_2/umap_prob_cluster_9*10),
               color=score_types['umap_cluster_2'],
               size=1)
    }    
    if("umap_cluster_9" %in% names(score_types)){
        cat("Adding umap_cluster_9 layer with color ", score_types['umap_cluster_9'], "\n", sep="")
        plot <- plot +
            ggplot2::geom_line(
               mapping=ggplot2::aes(x=dose_nM, y=umap_prob_cluster_9),
               color=score_types['umap_cluster_9'],
               size=1)
    }    
    if("rf_score" %in% names(score_types)){
        cat("Adding RF score layer with color ", score_types['rf_score'], "\n", sep="")
        plot <- plot +
#            ggplot2::geom_errorbar(
#               mapping=ggplot2::aes(
#                 x=dose_nM,
#                 ymin=tf_score_low,
#                 ymax=tf_score_high),
#               color=score_types['rf_score'],
#               width=.1) +
            ggplot2::geom_jitter(
               data=rf_scores_points %>% dplyr::filter(Compound != "PC", Compound != "NC"),  
               mapping=ggplot2::aes(x=dose_nM, y=rf_score),
               size=.3, color=score_types['rf_score'],
               alpha=1,
               width=0.05, height=0) +
        ggplot2::geom_line(
            data=fits,
            mapping=ggplot2::aes(x=dose_nM, y=pred_value),
            color=score_types['rf_score'])
    }
    plot <- plot +
        ggplot2::facet_wrap(~Compound) +
        ggplot2::ggtitle(
           label="Dose Response",
           subtitle=subtitle) +
        ggplot2::scale_x_log10(
           "Dose nM",
           breaks=c(4, 7.99, 16.0, 31.2, 62.3, 125, 256, 496, 991, 1998),
           labels=c("4", "8", "16", "31", "62", "125", "256", "496", "991", "1998")) +
        ggplot2::scale_y_continuous(
           paste0(score_types, ": ", names(score_types), collapse="   "))
}



p <- make_dose_response_plot(
    well_scores = well_scores %>%
        dplyr::filter(!Compound %in% c("NC", "PC", "BLANK")),
    subtitle="20XX Plates",
    score_types=c(
        cell_count="black",
        nucleoli_count="green",
        rf_score="red"),
    rf_scores_points=rf_scores_points %>% dplyr::filter(Compound != "Metformin"),
    fits=sigmoid_fits)


ggplot2::ggsave(
    "product/figures/plate_20XX_score_cell_nucleoli_rf_score_dose_response_200515.pdf",
    width=15, height=11)


##############
## Top Hits ##
##############

## RF Scores

top_hits <- readr::read_tsv("raw_data/top_hits_200513.tsv")

p <- make_dose_response_plot(
    well_scores = well_scores %>%
        dplyr::semi_join(top_hits, by=c("Compound")) %>%
        dplyr::mutate(Compound = Compound %>% stringr::str_extract("^[A-Za-z0-9-]+")),
    subtitle="Top Hits",
    score_types=c(
        cell_count="black",
        rf_score="magenta",
        umap_cluster_1="blue",
        umap_cluster_2="green",
        umap_cluster_9="grey"),
    rf_scores_points=rf_scores_points %>%
        dplyr::semi_join(top_hits, by=c("Compound")) %>%
        dplyr::mutate(Compound = Compound %>% stringr::str_extract("^[A-Za-z0-9-]+")),
    fits=sigmoid_fits %>%
        dplyr::semi_join(top_hits, by=c("Compound")) %>%
        dplyr::mutate(Compound = Compound %>% stringr::str_extract("^[A-Za-z0-9-]+")))        

p <- p + ggplot2::facet_wrap(~Compound, ncol=1, strip.position="left")
ggplot2::ggsave(
    "product/figures/dose_response/top_hits_score_cell_rf_score_dose_response_200515.pdf",
    width=3, height=15)


p <- p + ggplot2::facet_wrap(~Compound, ncol=3, strip.position="left")
ggplot2::ggsave(
    "product/figures/dose_response/top_hits_score_cell_rf_score_dose_response_3_col_200515.pdf",
    width=6, height=6)


##stratominer score

top_hits <- readr::read_tsv("raw_data/top_hits_200513.tsv")

p <- make_dose_response_plot(
    well_scores = well_scores %>%
        dplyr::semi_join(top_hits, by=c("Compound")) %>%
        dplyr::mutate(Compound = Compound %>% stringr::str_extract("^[A-Za-z0-9-]+")),
    subtitle="Top Hits",
    score_types=c(
        cell_count="black",
        stratominer_score="magenta"),
    rf_scores_points=rf_scores_points %>%
        dplyr::semi_join(top_hits, by=c("Compound")) %>%
        dplyr::mutate(Compound = Compound %>% stringr::str_extract("^[A-Za-z0-9-]+")),
    fits=sigmoid_fits %>%
        dplyr::semi_join(top_hits, by=c("Compound")) %>%
        dplyr::mutate(Compound = Compound %>% stringr::str_extract("^[A-Za-z0-9-]+")))        

p <- p + ggplot2::facet_wrap(~Compound, ncol=1, strip.position="left")
    
ggplot2::ggsave(
    "product/figures/dose_response/top_hits_score_cell_rf_score_dose_response_200514.pdf",
    width=3, height=15)
