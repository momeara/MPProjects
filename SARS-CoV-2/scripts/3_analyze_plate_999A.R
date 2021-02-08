
library(plyr)
library(tidyverse)
library(readxl)

plate_map <- readxl::read_excel("raw_data/Plate_999A_Metadata.xlsx")

# field level averages and has the number of viral positive cells.
data <- readr::read_csv("raw_data/MyExpt_Image.csv")

meta_data <- data %>%
    dplyr::select(
               dplyr::starts_with("Metadata"),
               dplyr::starts_with("Classify"),
               dplyr::starts_with("Count")) %>%
    dplyr::mutate(
      row = Metadata_WellID %>%
          stringr::str_extract("^[A-Z]") %>%
          purrr::map_int(~which(LETTERS==., arr.ind=T)),
      column = Metadata_WellID %>%
          stringr::str_extract("[0-9]+$") %>%
          as.integer())
               

meta_data <- meta_data %>%
    dplyr::left_join(
               plate_map %>%
               dplyr::select(-Plate_ID),
               by=c("Metadata_WellID" = "Well_ID"))


field_counts <- meta_data %>%
    dplyr::select(
               Metadata_WellID,
               row, column,
               Metadata_Field,
               Compound,
               Concentration,
               dplyr::starts_with("Count")) %>%
    tidyr::pivot_longer(
               cols=dplyr::starts_with("Count"),
               names_to="object_type",
               names_prefix="Count_",
               values_to="count") %>%
    dplyr::filter(!Compound %in% c("PC", "NC"))

### smoothed dose-response summaries
p <- ggplot2::ggplot(field_counts) +
    ggplot2::theme_bw() +
    ggplot2::geom_smooth(
        mapping=ggplot2::aes(
                             x=Concentration,
                             y=count,
                             color=object_type)) +
    ggplot2::facet_grid(
                 object_type~Compound,
                 scales="free")

ggsave("product/Plate_999A_object_counts_dose_response_200418.pdf",
       width=10, height=10)

##############################################################

classify_summary <- meta_data %>%
    dplyr::select(
               Metadata_WellID,
               Metadata_Field,
               Compound,
               Concentration,
               dplyr::starts_with("Classify")) %>%
    tidyr::pivot_longer(
               cols=dplyr::starts_with("Classify"),
               names_to="object_type",
               names_prefix="Classify_",
               values_to="classify") %>%
    dplyr::filter(!Compound %in% c("PC", "NC"))


### smoothed dose-response summaries
p <- ggplot2::ggplot(classify_summary) +
    ggplot2::theme_bw() +
    ggplot2::geom_smooth(
        mapping=ggplot2::aes(
                             x=Concentration,
                             y=classify,
                             color=classify)) +
    ggplot2::facet_grid(
                 object_type~Compound,
                 scales="free")

ggsave("product/Plate_999A_classify_dose_response_200418.pdf",
       width=10, height=10)

################################################################

classify_summary <- meta_data %>%
    dplyr::select(
               Metadata_WellID,
               Metadata_Field,
               Compound,
               Concentration,
               dplyr::starts_with("Classify")) %>%
    tidyr::pivot_longer(
               cols=dplyr::starts_with("Classify"),
               names_to="object_type",
               names_prefix="Classify_",
               values_to="classify") %>%
    dplyr::filter(!Compound %in% c("PC", "NC"))


### smoothed dose-response summaries
p <- ggplot2::ggplot(classify_summary) +
    ggplot2::theme_bw() +
    ggplot2::geom_smooth(
        mapping=ggplot2::aes(
                             x=Concentration,
                             y=classify,
                             color=classify)) +
    ggplot2::facet_grid(
                 object_type~Compound,
                 scales="free")

ggsave("product/Plate_999A_classify_dose_response_200418.pdf",
       width=10, height=10)

##################################################################3

model_data <- data %>%
    dplyr::left_join(
               plate_map %>%
               dplyr::select(-Plate_ID),
               by=c("Metadata_WellID" = "Well_ID")) %>%
    dplyr::filter(Compound %in% c("PC", "NC")) %>%
    dplyr::mutate(Compound = factor(Compound)) %>%
    dplyr::select(
               Metadata_WellID,
               Metadata_Field,
               Compound,
               dplyr::starts_with("Mean")) %>%
    tidyr::pivot_longer(
               cols=dplyr::starts_with("Mean"),
               names_to="feature",
               names_prefix="Mean_",
               values_to="value") %>%
    na.omit()


feature_summaries <- model_data %>%
    dplyr::group_by(feature) %>%
    dplyr::summarize(
        n_PC = sum(Compound == "PC"),
        n_NC = sum(Compound == "NC"),
        std_dev = sd(value)) %>%
    dplyr::filter(
               !(feature %>% stringr::str_detect("Parent")),
               !(feature %>% stringr::str_detect("Center_")),
               n_PC > 64,
               n_NC > 64,
               std_dev > 0)
model_data <- model_data %>%
    dplyr::semi_join(
               feature_summaries,
               by="feature")
    

s <- model_data %>%
    dplyr::group_by(feature) %>%
    dplyr::mutate(value = scale(value)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Compound, feature) %>%
    dplyr::summarize(
               mean = mean(value),
               sd = sd(value),
               min = min(value),
               max = max(value)) %>%
    tidyr::pivot_wider(
        id_cols=feature,
        names_from=Compound,       
        values_from=c(mean, sd, min, max)) %>%
    dplyr::mutate(zprime = 1 - 3 * (sd_PC + sd_NC/(abs(mean_PC - mean_NC))))

top_features <- s %>%
    dplyr::arrange(desc(zprime)) %>%
    head(20)



###########

s10 <- model_data %>%
    dplyr::filter(as.numeric(Metadata_Field) < 10) %>%
    dplyr::group_by(feature) %>%
    dplyr::mutate(value = scale(value)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Compound, feature) %>%
    dplyr::summarize(
               mean = mean(value),
               sd = sd(value),
               min = min(value),
               max = max(value)) %>%
    tidyr::pivot_wider(
        id_cols=feature,
        names_from=Compound,       
        values_from=c(mean, sd, min, max)) %>%
    dplyr::mutate(zprime = 1 - 3 * (sd_PC + sd_NC/(abs(mean_PC - mean_NC))))

s <- s %>%
    dplyr::mutate(zprime16 = zprime) %>%
    dplyr::left_join(
               s10 %>% dplyr::select(feature, zprime10=zprime),
               by="feature")
               



top_model_data <- model_data %>%
    dplyr::semi_join(top_features, by="feature") %>%
    tidyr::pivot_wider(
               id_cols=c("Metadata_WellID", "Metadata_Field", "Compound"),
               names_from=feature,
               values_from=value)


model_data <- model_data %>%
    tidyr::pivot_wider(
               id_cols=c("Metadata_WellID", "Metadata_Field", "Compound"),
               names_from=feature,
               values_from=value)

z <- model_data %>%
    dplyr::select(
         Compound,
         Mean_Cells_Intensity_MinIntensityEdge_Virus,
         Mean_Cells_Intensity_MedianIntensity_Lipids) %>%
    dplyr::filter(!is.na(Mean_Cells_Intensity_MedianIntensity_Lipids))

model <- glm(Compound ~ Mean_Cells_Intensity_MinIntensityEdge_Virus + Mean_Cells_Intensity_MedianIntensity_Lipids, data=z, family = "binomial")
model <- glm(Compound ~ ., data=z, family = "binomial")
model <- glm(Compound ~ ., data=model_data, family="binomial")

###########################

model <- glm(
    formula=Compound ~ .,
    data=top_model_data %>% dplyr::select(-Metadata_WellID, -Metadata_Field),
    family="binomial") %>%
    broom::


library(glmnet)
lasso_best <- glmnet(
    x=top_model_data %>%
        dplyr::select(-Metadata_WellID, -Metadata_Field, -Compound) %>%
        as.matrix(),
    y=top_model_data %>%
        dplyr::select(Compound) %>%
        as.matrix(),
    family="binomial",
    alpha = .2)


pdf("product/Plate_999A_top_features_by_zprime_200418.pdf")
plotmo::plot_glmnet(lasso_best)
dev.off()
