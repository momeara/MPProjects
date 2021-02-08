library(plyr)
library(tidyverse)
library(fuzzyjoin)
library(ggplot2)
library(readxl)
library(MPStats)

Zprime <- function(positives, negatives, na.rm=TRUE) {
    1 - 3 * (sd(positives, na.rm=na.rm) + sd(negatives, na.rm=na.rm))/abs(mean(positives, na.rm=na.rm) - mean(negatives, na.rm=na.rm))
}



feature_columns <- readr::read_tsv("raw_data/cell_feature_columns.tsv")
image_scores_20XX <- arrow::read_parquet("product/image_scores_20XX.parquet")
meta_well_scores <- image_scores_20XX %>%
    dplyr::distinct(
        master_plate_id,
        Compound,
        dose_nM)

plate_id <- "2006A"
cell_features <- arrow::read_parquet(
    file=paste0("product/SARS_", plate_id, "_Cell_MasterDataTable.parquet"),
    col_select=c("plate_id", "Condition", "Image_Metadata_WellID", feature_columns$feature)) %>%
    dplyr::filter(Condition %in% c("PC", "NC")) %>%
    dplyr::mutate(
        infectivity_score =
            -5.064328 +
            Cells_Intensity_IntegratedIntensityEdge_Virus * 1.487025e-01 +
            Cells_Intensity_MeanIntensityEdge_Virus * -3.840196e+01 +
            Cells_Intensity_MaxIntensityEdge_Virus * 4.270269e+01 +
            Cells_Intensity_MaxIntensity_Virus * 4.254849e+01)

feature_columns <- feature_columns %>%
    dplyr::bind_rows(
        data.frame(feature="infectivity_score", transform="identity"))

Zprime <- function(positives, negatives) {
    1 - 3 * (sd(positives) + sd(negatives))/abs(mean(positives) - mean(negatives))
}

zprime_scores <- plyr::ldply(feature_columns$feature, function(feature_id){   
  positives <- cell_features %>%
    dplyr::filter(Condition == "PC") %>%
    magrittr::extract2(feature_id)
  negatives <- cell_features %>%
    dplyr::filter(Condition == "NC") %>%
    magrittr::extract2(feature_id)
  data.frame(
    feature_id=feature_id,
    Zprime=Zprime(positives, negatives))
})

                                              feature_id        Zprime
1                      Nuclei_Intensity_MinIntensity_CMO     -3.344305
2                  Nuclei_Intensity_MinIntensityEdge_CMO     -3.402473
3                   Nuclei_Intensity_MinIntensity_Lipids     -3.477715
4               Nuclei_Intensity_MinIntensityEdge_Lipids     -3.592840
5         Nuclei_Intensity_LowerQuartileIntensity_Lipids     -3.878905
6                 Nuclei_Intensity_MeanIntensityEdge_CMO     -4.005551
7            Nuclei_Intensity_LowerQuartileIntensity_CMO     -4.141382
8                Nuclei_Intensity_MedianIntensity_Lipids     -4.240294
9                  Nuclei_Intensity_MeanIntensity_Lipids     -4.310019
10                    Cells_Intensity_MinIntensity_Virus     -4.343327
11             Cytoplasm_Intensity_MeanIntensityEdge_CMO     -4.368224
12             Nuclei_Intensity_MeanIntensityEdge_Lipids     -4.432672
13                    Nuclei_Intensity_MeanIntensity_CMO     -4.456870
14                  Nuclei_Intensity_MedianIntensity_CMO     -4.574998
15                 Cytoplasm_Intensity_MeanIntensity_CMO     -4.669878
16            Cytoplasm_Intensity_MinIntensityEdge_Virus     -4.717068
17                Cytoplasm_Intensity_MinIntensity_Virus     -4.738681
...
                                       infectivity_score     -6.718468

## averaging to the well level first
zprime_scores_well_mean <- plyr::ldply(feature_columns$feature, function(feature_id){   
  positives <- cell_features %>%
    dplyr::filter(Condition == "PC") %>%
    dplyr::group_by(Image_Metadata_WellID) %>%
    dplyr::summarize(mean_feature_value=mean(!!sym(feature_id))) %>%  
    magrittr::extract2("mean_feature_value")
  negatives <- cell_features %>%
    dplyr::filter(Condition == "NC") %>%
    dplyr::group_by(Image_Metadata_WellID) %>%
    dplyr::summarize(mean_feature_value=mean(!!sym(feature_id))) %>%  
    magrittr::extract2("mean_feature_value")
  data.frame(
    feature_id=feature_id,
    Zprime=Zprime(positives, negatives))
})

#                                                feature_id        Zprime
#  1             Cells_Intensity_LowerQuartileIntensity_CMO    -0.1013488
#  2                    Cells_Intensity_MedianIntensity_CMO    -0.1127203
#  3                   Cytoplasm_Intensity_MinIntensity_CMO    -0.1537897
#  4                       Cells_Intensity_MinIntensity_CMO    -0.1538760
#  5                Cytoplasm_Intensity_MedianIntensity_CMO    -0.1548385
#  6               Cells_Intensity_MeanIntensityEdge_Lipids    -0.1548388
#  7                      Cells_Intensity_MeanIntensity_CMO    -0.1561189
#  8                  Cells_Intensity_MeanIntensityEdge_CMO    -0.1598360
#  9                  Cytoplasm_Intensity_MeanIntensity_CMO    -0.1628205
#  10              Cytoplasm_Intensity_MeanIntensity_Lipids    -0.1643084
# ...
#  142                                    infectivity_score    -0.9883708



##########################################
## Zprime for RF SCores for 20XX Series ##
##########################################


load("intermediate_data/rf_scores_field_10XX.Rdata")

# straight fraction with viral intensity > .01
zprime_plate_10XX <- rf_scores_field_10XX %>%
    plyr::ddply(c("Plate_Name"), function(rf_scores){
        Zprime(
            positives = rf_scores %>%
                dplyr::filter(Compound == "PC") %>%
                magrittr::extract2("Image_Classify_Positive_PctObjectsPerBin"),
            negatives = rf_scores %>%
                dplyr::filter(Compound == "NC") %>%
                magrittr::extract2("Image_Classify_Positive_PctObjectsPerBin")) %>%
            data.frame(Zprime=.)
    })
zprime_plate_10XX %>% dplyr::summarize(mean(Zprime, na.rm=TRUE))

# mean over fields
zprime_plate_10XX <- rf_scores_field_10XX %>%
    dplyr::filter(Compound %in% c("PC", "NC")) %>%
    dplyr::group_by(Plate_Name, Compound, Well_ID) %>%
    dplyr::summarize(infectivity_probpos_well = mean(Image_Classify_Positive_PctObjectsPerBin)) %>%
    dplyr::ungroup() %>%
    plyr::ddply(c("Plate_Name"), function(rf_scores){
        MPStats::Zprime(
            positives = rf_scores %>%
                dplyr::filter(Compound == "PC") %>%
                magrittr::extract2("infectivity_probpos_well"),
            negatives = rf_scores %>%
                dplyr::filter(Compound == "NC") %>%
                magrittr::extract2("infectivity_probpos_well")) %>%
            data.frame(Zprime=.)
    })
zprime_plate_10XX %>% dplyr::summarize(mean(Zprime, na.rm=TRUE))



# RF stabilized model over frames
zprime_plate_10XX <- rf_scores_field_10XX %>%
    plyr::ddply(c("Plate_Name"), function(rf_scores){
        MPStats::Zprime(
            positives = rf_scores %>%
                dplyr::filter(Compound == "PC") %>%
                magrittr::extract2("infectivity_probpos_field"),
            negatives = rf_scores %>%
                dplyr::filter(Compound == "NC") %>%
                magrittr::extract2("infectivity_probpos_field")) %>%
            data.frame(Zprime=.)
    })
zprime_plate_10XX %>% dplyr::summarize(mean(Zprime, na.rm=TRUE))

# mean over fields
zprime_plate_10XX <- rf_scores_field_10XX %>%
    dplyr::filter(Compound %in% c("PC", "NC")) %>%
    dplyr::group_by(Plate_Name, Compound, Well_ID) %>%
    dplyr::summarize(infectivity_probpos_well = mean(infectivity_probpos_field)) %>%
    dplyr::ungroup() %>%
    plyr::ddply(c("Plate_Name"), function(rf_scores){
        MPStats::Zprime(
            positives = rf_scores %>%
                dplyr::filter(Compound == "PC") %>%
                magrittr::extract2("infectivity_probpos_well"),
            negatives = rf_scores %>%
                dplyr::filter(Compound == "NC") %>%
                magrittr::extract2("infectivity_probpos_well")) %>%
            data.frame(Zprime=.)
    })
zprime_plate_10XX %>% dplyr::summarize(mean(Zprime, na.rm=TRUE))


# median over fields
zprime_plate_10XX <- rf_scores_field_10XX %>%
    dplyr::filter(Compound %in% c("PC", "NC")) %>%
    dplyr::group_by(Plate_Name, Compound, Well_ID) %>%
    dplyr::summarize(infectivity_probpos_well = median(infectivity_probpos_field)) %>%
    dplyr::ungroup() %>%
    plyr::ddply(c("Plate_Name"), function(rf_scores){
        MPStats::Zprime(
            positives = rf_scores %>%
                dplyr::filter(Compound == "PC") %>%
                magrittr::extract2("infectivity_probpos_well"),
            negatives = rf_scores %>%
                dplyr::filter(Compound == "NC") %>%
                magrittr::extract2("infectivity_probpos_well")) %>%
            data.frame(Zprime=.)
    })
zprime_plate_10XX %>% dplyr::summarize(mean(Zprime, na.rm=TRUE))


# max over fields
zprime_plate_10XX <- rf_scores_field_10XX %>%
    dplyr::filter(Compound %in% c("PC", "NC")) %>%
    dplyr::group_by(Plate_Name, Compound, Well_ID) %>%
    dplyr::summarize(infectivity_probpos_well = max(infectivity_probpos_field)) %>%
    dplyr::ungroup() %>%
    plyr::ddply(c("Plate_Name"), function(rf_scores){
        MPStats::Zprime(
            positives = rf_scores %>%
                dplyr::filter(Compound == "PC") %>%
                magrittr::extract2("infectivity_probpos_well"),
            negatives = rf_scores %>%
                dplyr::filter(Compound == "NC") %>%
                magrittr::extract2("infectivity_probpos_well")) %>%
            data.frame(Zprime=.)
    })
zprime_plate_10XX %>% dplyr::summarize(mean(Zprime, na.rm=TRUE))


##########################################
## Zprime for RF SCores for 20XX Series ##
##########################################



zprime_plate_20XX <- rf_scores_field_20XX %>%
    plyr::ddply(c("Plate_Name"), function(rf_scores){
        MPStats::Zprime(
            positives = rf_scores %>%
                dplyr::filter(Compound == "PC") %>%
                magrittr::extract2("infectivity_probpos_field"),
            negatives = rf_scores %>%
                dplyr::filter(Compound == "NC") %>%
                magrittr::extract2("infectivity_probpos_field")) %>%
            data.frame(Zprime=.)
    })
zprime_plate_20XX %>% dplyr::summarize(mean(Zprime, na.rm=TRUE))

# mean over fields
zprime_plate_20XX <- rf_scores_field_20XX %>%
    dplyr::filter(Compound %in% c("PC", "NC")) %>%
    dplyr::group_by(Plate_Name, Compound, Well_ID) %>%
    dplyr::summarize(infectivity_probpos_well = mean(infectivity_probpos_field)) %>%
    dplyr::ungroup() %>%
    plyr::ddply(c("Plate_Name"), function(rf_scores){
        MPStats::Zprime(
            positives = rf_scores %>%
                dplyr::filter(Compound == "PC") %>%
                magrittr::extract2("infectivity_probpos_well"),
            negatives = rf_scores %>%
                dplyr::filter(Compound == "NC") %>%
                magrittr::extract2("infectivity_probpos_well")) %>%
            data.frame(Zprime=.)
    })
zprime_plate_20XX %>% dplyr::summarize(mean(Zprime, na.rm=TRUE))


# median over fields
zprime_plate_20XX <- rf_scores_field_20XX %>%
    dplyr::filter(Compound %in% c("PC", "NC")) %>%
    dplyr::group_by(Plate_Name, Compound, Well_ID) %>%
    dplyr::summarize(infectivity_probpos_well = median(infectivity_probpos_field)) %>%
    dplyr::ungroup() %>%
    plyr::ddply(c("Plate_Name"), function(rf_scores){
        MPStats::Zprime(
            positives = rf_scores %>%
                dplyr::filter(Compound == "PC") %>%
                magrittr::extract2("infectivity_probpos_well"),
            negatives = rf_scores %>%
                dplyr::filter(Compound == "NC") %>%
                magrittr::extract2("infectivity_probpos_well")) %>%
            data.frame(Zprime=.)
    })
zprime_plate_20XX %>% dplyr::summarize(mean(Zprime, na.rm=TRUE))




