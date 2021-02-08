
library(plyr)
library(tidyverse)
library(arrow)
library(caret)
library(MPstats)

source("scripts/make_scatter_boards.R")


single_data_path <- "intermediate_data/infected_patch_999A_20201112"
single_viral_features <- arrow::read_parquet(
    file = paste0(data_path, "/viral_features.parquet"))
single_viral_feature_columns <- readr::read_tsv(
    file = paste0(data_path, "/viral_feature_columns.tsv"))
single_viral_metadata_columns <- readr::read_tsv(
    file = paste0(data_path, "/viral_metadata_columns.tsv"))


combo_data_path <- "intermediate_data/infected_patch_1999B_2020A_2021A_20201017"
combo_viral_features <- arrow::read_parquet(
    file = paste0(data_path, "/viral_features.parquet"))
combo_viral_feature_columns <- readr::read_tsv(
    file = paste0(data_path, "/viral_feature_columns.tsv"))
combo_viral_metadata_columns <- readr::read_tsv(
    file = paste0(data_path, "/viral_metadata_columns.tsv"))


viral_embed_features <- viral_features %>%
    dplyr::filter(condition != "BLANK") %>%
    dplyr::mutate(
        dplyr::across(
            .cols = viral_feature_columns$feature,
            .fns = function(feature) {
                transform <- viral_feature_columns %>%
                    dplyr::filter(feature == dplyr::cur_column()) %>%
                    dplyr::pull(transform)
                if (transform == "log") {
                    feature <- log(feature)
                } else if (transform == "log1p") {
                    feature <- log(feature + 1)
                } else if (transform == "identity") {
                    feature <- feature
                } else {
                    stop(paste0("Unrecongized transform '", transform, "'", sep = ""))
                }
            }),
        .keep = "unused") %>%
    dplyr::group_by(plate_id) %>%
    dplyr::mutate_at(viral_feature_columns$feature, ~ scale(.)[, 1]) %>%
    dplyr::ungroup()

viral_embed_features %>%
    arrow::write_parquet(
        sink = paste0(
            data_path,
            "/viral_plate_scaled_MasterDataTable.parquet"))


#######################
# Embed Viral Patches #
#######################
system(paste0("
        cd ~/opt/MPLearn/vignettes/SARS-CoV-2/S25 &&
        /home/ubuntu/anaconda3/envs/sextonlab/bin/python \\
            ~/anaconda3/envs/sextonlab/bin/embed_umap \\
            --dataset ~/opt/MPStats/vignettes/SARS-CoV-2/", data_path, "/viral_plate_scaled_MasterDataTable.parquet \\
            --tag UMAP_viral_plate_scaled_n_neighbors=5_20201019 \\
            --feature_columns ~/opt/MPStats/vignettes/SARS-CoV-2/", data_path, "/viral_feature_no_transform_columns.tsv \\
            --no_standardize_features \\
            --umap_n_neighbors 5 \\
	    --umap_low_memory \\
	    --verbose
"))

system(paste0("
        cd ~/opt/MPLearn/vignettes/SARS-CoV-2/S25 &&
        /home/ubuntu/anaconda3/envs/sextonlab/bin/python \\
            ~/anaconda3/envs/sextonlab/bin/embed_umap \\
            --dataset ~/opt/MPStats/vignettes/SARS-CoV-2/", data_path, "/viral_plate_scaled_MasterDataTable.parquet \\
            --tag UMAP_viral_plate_scaled_n_neighbors=30_umap_negative_sample_rate=50_20201019 \\
            --feature_columns ~/opt/MPStats/vignettes/SARS-CoV-2/", data_path, "/viral_feature_no_transform_columns.tsv \\
            --no_standardize_features \\
            --umap_n_neighbors 30 \\
            --umap_negative_sample_rate 50 \\
	    --umap_low_memory \\
	    --verbose
"))


system(paste0("
        cd ~/opt/MPLearn/vignettes/SARS-CoV-2/S25 &&
        /home/ubuntu/anaconda3/envs/sextonlab/bin/python \\
            ~/anaconda3/envs/sextonlab/bin/embed_umap \\
            --dataset ~/opt/MPStats/vignettes/SARS-CoV-2/", data_path, "/viral_plate_scaled_MasterDataTable.parquet \\
            --tag UMAP_viral_plate_scaled_n_neighbors=30_umap_negative_sample_rate=50_a=2_b=.2_20201019 \\
            --feature_columns ~/opt/MPStats/vignettes/SARS-CoV-2/", data_path, "/viral_feature_no_transform_columns.tsv \\
            --no_standardize_features \\
            --umap_n_neighbors 30 \\
            --umap_negative_sample_rate 50 \\
            --umap_a 2 \\
            --umap_b .2 \\
	    --umap_low_memory \\
	    --verbose
"))

system(paste0("
        cd ~/opt/MPLearn/vignettes/SARS-CoV-2/S25 &&
        /home/ubuntu/anaconda3/envs/sextonlab/bin/python \\
            ~/anaconda3/envs/sextonlab/bin/embed_umap \\
            --dataset ~/opt/MPStats/vignettes/SARS-CoV-2/", data_path, "/viral_plate_scaled_MasterDataTable.parquet \\
            --tag UMAP_viral_plate_scaled_n_neighbors=30_umap_negative_sample_rate=50_a=2_b=.6_20201019 \\
            --feature_columns ~/opt/MPStats/vignettes/SARS-CoV-2/", data_path, "/viral_feature_no_transform_columns.tsv \\
            --no_standardize_features \\
            --umap_n_neighbors 30 \\
            --umap_negative_sample_rate 50 \\
            --umap_a 2 \\
            --umap_b .6 \\
	    --umap_low_memory \\
	    --verbose
"))



system(paste0("
        cd ~/opt/MPLearn/vignettes/SARS-CoV-2/S25 &&
        /home/ubuntu/anaconda3/envs/sextonlab/bin/python \\
            ~/anaconda3/envs/sextonlab/bin/embed_umap \\
            --dataset ~/opt/MPStats/vignettes/SARS-CoV-2/", data_path, "/viral_plate_scaled_MasterDataTable.parquet \\
            --tag UMAP_viral_plate_scaled_n_neighbors=30_umap_negative_sample_rate=50_a=20_b=.6_20201019 \\
            --feature_columns ~/opt/MPStats/vignettes/SARS-CoV-2/", data_path, "/viral_feature_no_transform_columns.tsv \\
            --no_standardize_features \\
            --umap_n_neighbors 30 \\
            --umap_negative_sample_rate 50 \\
            --umap_a 20 \\
            --umap_b .6 \\
	    --umap_low_memory \\
	    --verbose
"))

system(paste0("
        cd ~/opt/MPLearn/vignettes/SARS-CoV-2/S25 &&
        /home/ubuntu/anaconda3/envs/sextonlab/bin/python \\
            ~/anaconda3/envs/sextonlab/bin/embed_umap \\
            --dataset ~/opt/MPStats/vignettes/SARS-CoV-2/", data_path, "/viral_plate_scaled_MasterDataTable.parquet \\
            --tag UMAP_viral_plate_scaled_n_neighbors=30_umap_negative_sample_rate=50_a=50_b=.6_20201019 \\
            --feature_columns ~/opt/MPStats/vignettes/SARS-CoV-2/", data_path, "/viral_feature_no_transform_columns.tsv \\
            --no_standardize_features \\
            --umap_n_neighbors 30 \\
            --umap_negative_sample_rate 50 \\
            --umap_a 50 \\
            --umap_b .6 \\
	    --umap_low_memory \\
	    --verbose
"))

system(paste0("
        cd ~/opt/MPLearn/vignettes/SARS-CoV-2/S25 &&
        /home/ubuntu/anaconda3/envs/sextonlab/bin/python \\
            ~/anaconda3/envs/sextonlab/bin/embed_umap \\
            --dataset ~/opt/MPStats/vignettes/SARS-CoV-2/", data_path, "/viral_plate_scaled_MasterDataTable.parquet \\
            --tag UMAP_viral_no_zp0_plate_scaled_n_neighbors=30_umap_negative_sample_rate=50_a=50_b=.2_20201019 \\
            --feature_columns ~/opt/MPStats/vignettes/SARS-CoV-2/", data_path, "/viral_feature_no_transform_columns.tsv \\
            --no_standardize_features \\
            --umap_n_neighbors 30 \\
            --umap_negative_sample_rate 50 \\
            --umap_a 50 \\
            --umap_b .2 \\
	    --umap_low_memory \\
	    --verbose
"))



system(paste0("
        cd ~/opt/MPLearn/vignettes/SARS-CoV-2/S25 &&
        /home/ubuntu/anaconda3/envs/sextonlab/bin/python \\
            ~/anaconda3/envs/sextonlab/bin/embed_umap \\
            --dataset ~/opt/MPStats/vignettes/SARS-CoV-2/", data_path, "/viral_plate_scaled_MasterDataTable.parquet \\
            --tag UMAP_viral_plate_scaled_n_neighbors=5_negative_sample_rate=50_20201019 \\
            --feature_columns ~/opt/MPStats/vignettes/SARS-CoV-2/", data_path, "/viral_feature_no_transform_columns.tsv \\
            --no_standardize_features \\
            --umap_n_neighbors 5 \\
            --umap_negative_sample_rate 50 \\
	    --umap_low_memory \\
	    --verbose
"))


system(paste0("
        cd ~/opt/MPLearn/vignettes/SARS-CoV-2/S25 &&
        /home/ubuntu/anaconda3/envs/sextonlab/bin/python \\
            ~/anaconda3/envs/sextonlab/bin/embed_umap \\
            --dataset ~/opt/MPStats/vignettes/SARS-CoV-2/", data_path, "/viral_plate_scaled_MasterDataTable.parquet \\
            --tag UMAP_viral_plate_scaled_n_neighbors=15_20201019 \\
            --feature_columns ~/opt/MPStats/vignettes/SARS-CoV-2/", data_path, "/viral_feature_no_transform_columns.tsv \\
            --no_standardize_features \\
            --umap_n_neighbors 15 \\
	    --umap_low_memory \\
	    --verbose
"))

system(paste0("
        cd ~/opt/MPLearn/vignettes/SARS-CoV-2/S25 &&
        /home/ubuntu/anaconda3/envs/sextonlab/bin/python \\
            ~/anaconda3/envs/sextonlab/bin/embed_umap \\
            --dataset ~/opt/MPStats/vignettes/SARS-CoV-2/", data_path, "/viral_plate_scaled_MasterDataTable.parquet \\
            --tag UMAP_viral_plate_scaled_n_nieghbors=30_20201019 \\
            --feature_columns ~/opt/MPStats/vignettes/SARS-CoV-2/", data_path, "/viral_feature_no_transform_columns.tsv \\
            --no_standardize_features \\
            --umap_n_neighbors 30 \\
	    --umap_low_memory \\
	    --verbose
"))

system(paste0("
        cd ~/opt/MPLearn/vignettes/SARS-CoV-2/S25 &&
        /home/ubuntu/anaconda3/envs/sextonlab/bin/python \\
            ~/anaconda3/envs/sextonlab/bin/embed_umap \\
            --dataset ~/opt/MPStats/vignettes/SARS-CoV-2/", data_path, "/viral_plate_scaled_MasterDataTable.parquet \\
            --tag UMAP_viral_plate_scaled_umap_n_neighbors=300_20201019 \\
            --feature_columns ~/opt/MPStats/vignettes/SARS-CoV-2/", data_path, "/viral_feature_no_transform_columns.tsv \\
            --no_standardize_features \\
            --umap_n_neighbors 300 \\
	    --umap_low_memory \\
	    --verbose
"))



# too smooth
system(paste0("
        cd ~/opt/MPLearn/vignettes/SARS-CoV-2/S25 &&
        /home/ubuntu/anaconda3/envs/sextonlab/bin/python \\
            ~/anaconda3/envs/sextonlab/bin/embed_umap \\
            --dataset ~/opt/MPStats/vignettes/SARS-CoV-2/", data_path, "/viral_plate_scaled_MasterDataTable.parquet \\
            --tag UMAP_viral_plate_scaled_neighbors=15_neg_sampling_rate=20_epochs=2000_20201019 \\
            --feature_columns ~/opt/MPStats/vignettes/SARS-CoV-2/", data_path, "/viral_feature_no_transform_columns.tsv \\
            --no_standardize_features \\
	    --umap_low_memory \\
            --umap_negative_sample_rate 20 \\
            --umap_n_epochs 2000 \\
	    --verbose
"))

# too many small clusters?
system(paste0("
        cd ~/opt/MPLearn/vignettes/SARS-CoV-2/S25 &&
        /home/ubuntu/anaconda3/envs/sextonlab/bin/python \\
            ~/anaconda3/envs/sextonlab/bin/embed_umap \\
            --dataset ~/opt/MPStats/vignettes/SARS-CoV-2/", data_path, "/viral_plate_scaled_MasterDataTable.parquet \\
            --tag UMAP_viral_plate_scaled_neighbors=5_neg_sampling_rate=200_epochs=2000_b=.2_20201019 \\
            --feature_columns ~/opt/MPStats/vignettes/SARS-CoV-2/", data_path, "/viral_feature_no_transform_columns.tsv \\
            --no_standardize_features \\
	    --umap_low_memory \\
            --umap_b .2 \\
            --umap_n_neighbors 5 \\
            --umap_negative_sample_rate 200 \\
            --umap_n_epochs 2000 \\
	    --verbose
"))


#################################################
# Use Monocle3 to cluster and analyze embedding #
#################################################

embedding_path <- "/home/ubuntu/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/UMAP_viral_no_zp0_plate_scaled_n_neighbors=30_umap_negative_sample_rate=50_a=50_b=.2_20201019"

viral_embed_features <- dplyr::bind_cols(
    arrow::read_parquet(paste0(data_path, "/viral_plate_scaled_MasterDataTable.parquet")),
    arrow::read_parquet(paste0(embedding_path, "/umap_embedding.parquet")))

infected_cds <- MPStats::populate_cds(
    cell_features = viral_embed_features,
    cell_feature_columns = viral_feature_columns,
    cell_metadata_columns = viral_metadata_columns,
    embedding_type = c("UMAP"),
    embedding = viral_embed_features %>% dplyr::select(UMAP_1, UMAP_2),
    verbose = TRUE)

# as resolution gets bigger --> more clusters
infected_cds <- infected_cds %>%
    monocle3::cluster_cells(
        reduction_method = "UMAP",
        k = 200,
        resolution = .00001,
        num_iter = 10,
        verbose = TRUE)
infected_cds %>% MPStats::serialize_clusters(
    output_fname = paste0(embedding_path, "/clusters_leiden_res=5e-5.parquet"))

# as resolution gets bigger --> more clusters
infected_cds <- infected_cds %>%
    monocle3::cluster_cells(
        reduction_method = "UMAP",
        k = 200,
        resolution = .0001,
        num_iter = 10,
        verbose = TRUE)
infected_cds %>% MPStats::serialize_clusters(
    output_fname = paste0(embedding_path, "/clusters_leiden_res=5e-4.parquet"))

system(paste0("cd ", embedding_path, " && ln -s clusters_leiden_res=5e-4.parquet clusters.parquet"))

# as resolution gets bigger --> more clusters
infected_cds <- infected_cds %>%
    monocle3::cluster_cells(
        reduction_method = "UMAP",
        k = 200,
        resolution = .001,
        num_iter = 10,
        verbose = TRUE)
infected_cds %>% MPStats::serialize_clusters(
    output_fname = paste0(embedding_path, "/clusters_leiden_res=5e-3.parquet"))

system(paste0("cd ", embedding_path, " && ln -s clusters_leiden_res=5e-3.parquet clusters.parquet"))



###############################################
# plot clusters distribution by drug and dose #
###############################################

viral_embed_features <- dplyr::bind_cols(
    arrow::read_parquet(paste0(data_path, "/viral_plate_scaled_MasterDataTable.parquet")),
    arrow::read_parquet(paste0(embedding_path, "/umap_embedding.parquet")),
    arrow::read_parquet(paste0(embedding_path, "/clusters.parquet")))
    



features <- viral_embed_features %>%
    dplyr::filter(
        UMAP_2 < 5, UMAP_2 > -8,
        UMAP_1 < 5)


make_scatter_boards(
    features = features,
    feature_x = "UMAP_1",
    feature_y = "UMAP_2",
    feature_color = "cluster_label",
    scales = NULL,
    output_dir = paste0(embedding_path, "/viral_features"))


make_scatter_boards(
    features = features,
    feature_x = "UMAP_1",
    feature_y = "UMAP_2",
    feature_color = "AreaShape_Eccentricity",
    scales = NULL,
    output_dir = paste0(embedding_path, "/viral_features"))



make_scatter_boards(
    features = features,
    feature_x = "UMAP_1",
    feature_y = "UMAP_2",
    feature_color = "AreaShape_Area",
    scales = NULL,
    output_dir = paste0(embedding_path, "/viral_features"))


make_scatter_boards(
    features = features,
    feature_x = "UMAP_1",
    feature_y = "UMAP_2",
    feature_color = "AreaShape_FormFactor",
    scales = NULL,
    output_dir = paste0(embedding_path, "/viral_features"))



make_scatter_boards(
    features = features,
    feature_x = "UMAP_1",
    feature_y = "UMAP_2",
    feature_color = "Children_Nuclei_Count",
    scales = NULL,
    output_dir = "product/infected_patch_1999B_2020A_2021A_20201017/viral_features")


make_scatter_boards(
    features = features,
    feature_x = "UMAP_1",
    feature_y = "UMAP_2",
    feature_color = "AreaShape_MinorAxisLength",
    scales = NULL,
    output_dir = "product/infected_patch_1999B_2020A_2021A_20201017/viral_features")


######################################################
# Discriminative features between different clusters #
######################################################



feature_importance <- viral_embed_features %>%
    dplyr::distinct(cluster_label) %>%
    plyr::adply(1, function(df) {
        cat("Computing feature importance for cluster '", df$cluster_label[1], "'\n", sep = "")
        data_matrix <- viral_embed_features  %>%
            dplyr::mutate(label = dplyr::case_when(
                as.character(cluster_label) != df$cluster_label[1] ~ "other",
                as.character(cluster_label) == df$cluster_label[1] ~ "cluster",
                TRUE ~ NA_character_) %>%
                    factor(
                        levels = c("cluster", "other"),
                        labels = c("cluster", "other"))) %>%
            dplyr::select(
                label,
                tidyselect::one_of(viral_feature_columns$feature)) %>%
            tidyr::drop_na()
        in_training <- caret::createDataPartition(
            data_matrix$label,
            p = .80,
            list = FALSE)
        data_train <- data_matrix %>% dplyr::slice(in_training)
        data_test  <- data_matrix %>% dplyr::slice(-in_training)
        cat("    N train: ", nrow(data_train), "\n", sep = "")
        cat("    N test:  ", nrow(data_test), "\n", sep = "")
        fit_control <- caret::trainControl(## 10-fold CV
            method = "repeatedcv",
            number = 2,
            classProbs = TRUE)
        pls_fit <- caret::train(
            label ~ .,
            data = data_train,
            method = "pls",
            trControl = fit_control,
            preProc = c("center", "scale"),
            tuneLength = 3,
            verbose = TRUE)
        pred_fit <- predict(pls_fit, newdata = data_test)
        feature_importance <- pls_fit %>% caret::varImp()
        feature_importance$importance %>%
            dplyr::bind_cols(
                viral_feature_columns %>%
                dplyr::select(feature)) %>%
            dplyr::rename(feature_score = Overall) %>%
            dplyr::arrange(desc(feature_score)) %>%
            head(20) %>%
            dplyr::mutate(cluster_label = df$cluster_label[1])
    })


feature_importance %>%
    dplyr::group_by(cluster_label) %>%
    dplyr::arrange(desc(feature_score)) %>%
    dplyr::slice(1:3) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(feature) %>%
    dplyr::rowwise() %>%
    dplyr::do({
        cat("Making plot for feature '", .$feature[1], "'\n", sep = "")
        make_scatter_boards(
            features = features,
            feature_x = "UMAP_1",
            feature_y = "UMAP_2",
            feature_color = .$feature[1],
            scales = NULL,
            output_dir = paste0(embedding_path, "/viral_features"))
        data.frame()
    })




viral_features <- dplyr::bind_cols(
    arrow::read_parquet(file = paste0(data_path, "/viral_features.parquet")) %>%
    dplyr::filter(condition != "BLANK"),
    arrow::read_parquet(paste0(embedding_path, "/umap_embedding.parquet")),
    arrow::read_parquet(paste0(embedding_path, "/clusters.parquet")))    

viral_features %>%
    dplyr::group_by(cluster_label) %>%
    dplyr::summarize(
        dplyr::across(
            c("AreaShape_Area", "AreaShape_FormFactor", "Texture_AngularSecondMoment_NP_8_00"),
            mean)) %>%
    dplyr::mutate(
        cluster_label = cluster_label %>%
            factor(            
                levels = c(5, 3, 2, 4, 1),
                labels = c("Cluster 5", "Cluster 3", "Cluster 2", "Cluster 4", "Cluster 1")))
        

data <- viral_features %>%
    dplyr::count(
        plate_id,
        well_id,
        condition,
        cluster_label,
        drug_1,
        drug_1_concentration,
        drug_1_units,
        drug_1_label,
        drug_2,
        drug_2_concentration,
        drug_2_units,
        drug_2_label,
        name = "viral_count_per_well")

data <- data %>%
    dplyr::filter(condition == "Treatment") %>%
    dplyr::left_join(
        data %>%
        dplyr::filter(condition == "NC") %>%
        dplyr::group_by(plate_id, cluster_label) %>%
        dplyr::summarize(
            mean_nc_count_per_well = mean(viral_count_per_well),
            .groups = "drop"),
        by = c("plate_id", "cluster_label")) %>%    
    dplyr::group_by(
        plate_id,
        condition,
        cluster_label,
        drug_1,
        drug_1_concentration,
        drug_1_units,
        drug_1_label,
        drug_2,
        drug_2_concentration,
        drug_2_units,
        drug_2_label) %>%
    dplyr::summarize(
        normalized_viral_count_per_well = mean(viral_count_per_well) / mean(mean_nc_count_per_well),
        .groups = "drop")

data %>%
    dplyr::group_by(plate_id) %>%
    dplyr::do({
        data_per_plate <- .
        plate_id <- data_per_plate$plate_id[1]
        data_per_plate <- data_per_plate %>%
            dplyr::mutate(
                drug_1_label = reorder(drug_1_label, drug_1_concentration),
                drug_2_label = reorder(drug_2_label, drug_2_concentration),
                cluster_label = paste0("Cluster ", cluster_label))
	ggplot2::ggplot() +
	    ggplot2::theme_bw() +
            ggplot2::theme(legend.position = c(.8, .18)) +
	    ggplot2::geom_line(
	        data = data_per_plate,
	        mapping = ggplot2::aes(
	            x = drug_2_concentration,
	            y = normalized_viral_count_per_well,
	            group = drug_1_label,
	            color = drug_1_label),
                size = 1.5) +
            ggplot2::scale_y_continuous("Cluster count normalized per-plate to NC") +
            ggplot2::scale_x_log10(paste0("Dose ", data_per_plate$drug_2[1], " (", data_per_plate$drug_1_units[1], ")")) +
            ggplot2::scale_color_discrete(data_per_plate$drug_1[1]) +
            ggplot2::ggtitle(
                label = "Viral objects -> UMAP clusters by dose",
                subtitle = paste0("Plate: ", data_per_plate$plate_id[1])) +
	    ggplot2::facet_wrap(
	        facets = dplyr::vars(cluster_label))
	ggsave(paste0(embedding_path, "/viral_features/counts_per_cluster_drug2_drug1_", plate_id, ".pdf"),
	    width = 7,
	    height = 5)
        data.frame()
    }) %>%
    dplyr::ungroup()


data %>%
    dplyr::group_by(plate_id) %>%
    dplyr::do({
        data_per_plate <- .
        plate_id <- data_per_plate$plate_id[1]
        data_per_plate <- data_per_plate %>%
            dplyr::mutate(
                drug_1_label = reorder(drug_1_label, drug_1_concentration),
                drug_2_label = reorder(drug_2_label, drug_2_concentration),
                cluster_label = paste0("Cluster ", cluster_label))
	ggplot2::ggplot() +
	    ggplot2::theme_bw() +
            ggplot2::theme(legend.position = c(.8, .18)) +
	    ggplot2::geom_line(
	        data = data_per_plate,
	        mapping = ggplot2::aes(
	            x = drug_1_concentration,
	            y = normalized_viral_count_per_well,
	            group = drug_2_label,
	            color = drug_2_label),
                size = 1.5) +
            ggplot2::scale_y_continuous("Cluster count normalized per-plate to NC") +
            ggplot2::scale_x_log10(paste0("Dose ", data_per_plate$drug_1[1], " (", data_per_plate$drug_1_units[1], ")")) +
            ggplot2::scale_color_discrete(data_per_plate$drug_2[1]) +
            ggplot2::ggtitle(
                label = "Viral objects -> UMAP clusters by dose",
                subtitle = paste0("Plate: ", data_per_plate$plate_id[1])) +
	    ggplot2::facet_wrap(
	        facets = dplyr::vars(cluster_label))
	ggsave(paste0(embedding_path, "/viral_features/counts_per_cluster_drug1_drug2_", plate_id, ".pdf"),
	    width = 7,
	    height = 5)
        data.frame()
    }) %>%
    dplyr::ungroup()
