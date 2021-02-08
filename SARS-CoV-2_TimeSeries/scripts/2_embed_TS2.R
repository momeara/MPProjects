


library(plyr)
library(tidyverse)
library(fuzzyjoin)
library(ggplot2)
library(readxl)
library(MPStats)
library(arrow)
library(rdist)

cell_feature_columns <- readr::read_tsv(
    "product/cell_feature_columns_TS_202008.tsv")

cell_metadata_columns <- readr::read_tsv(
    "product/cell_metadata_columns_TS_202008.tsv")



#################
# Load features #
#################

system("cp -r ~/bucket/UMAP_embeddings/top_hits_plate_scaled_200522a_umap2_2M_15_0.0 intermediate_data/")


# 6 time-points + unifected control replicated across 3 plates
# 2M / 21 total conditions = 95238 cells per condition
cell_features_2M <- tibble::tibble(plate_id = c("TS2PL1", "TS2PL2", "TS2PL3")) %>%
    dplyr::rowwise() %>%
    dplyr::do({
        plate_id <- .$plate_id[1]
        cat("Gathering cell features from plate id '", plate_id, "' ...\n", sep = "")
        arrow::read_parquet(
            paste0("product/", plate_id, "_Cell_MasterDataTable.parquet")) %>%
            dplyr::mutate_at(
                cell_feature_columns$feature, ~ scale(.) %>% as.vector()) %>%
            dplyr::group_by(Condition, time_point) %>%
            dplyr::sample_n(floor(2000000 / 21)) %>%
            dplyr::ungroup()
    }) %>%
    dplyr::ungroup()

cell_features_2M %>%
    arrow::write_parquet(
        sink = "product/TS2_2M_Cell_MasterDataTable.parquet")


cell_features <- NULL
cell_features_2M <- NULL
gc()

system(paste0("
        cd ", parameters$base_dir, " &&
        ", parameters$python_env, "/bin/python \\
            ", parameters$python_env, "/bin/embed_umap \\
            --dataset ", parameters$base_dir, "/product/covid19cq1_SARS_TS2_2M_Cell_MasterDataTable.parquet \\
            --tag UMAP_embedding_TS2_2M_epochs=2000_20200901 \\
            --feature_columns ", parameters$base_dir, "/product/cell_feature_columns_TS_202008.tsv \\
            --no_standardize_features \\
            --umap_low_memory \\
            --umap_n_epochs 2000 \\
	    --verbose
")

######################
# Re-embed all cells #
######################

cell_features <- tibble::tibble(plate_id = c("TS2PL1", "TS2PL2", "TS2PL3")) %>%
    dplyr::rowwise() %>%
    dplyr::do({
        plate_id <- .$plate_id[1]
        cat("Gathering cell features from plate id '", plate_id, "' ...\n", sep = "")
        arrow::read_parquet(
            paste0("product/", plate_id, "_Cell_MasterDataTable.parquet")) %>%
            dplyr::mutate_at(
                cell_feature_columns$feature, ~ scale(.) %>% as.vector())
    }) %>%
    dplyr::ungroup()

cell_features %>%
    arrow::write_parquet(
        sink = "product/TS2_plate_scaled_Cell_MasterDataTable.parquet")


system("
        cd ~/opt/MPLearn/vignettes/SARS-CoV-2/S25 &&
        /home/ubuntu/anaconda3/envs/sextonlab/bin/python \\
            ~/anaconda3/envs/sextonlab/bin/embed_umap \\
            --ref_embed_dir intermediate_data/UMAP_embedding_TS2_2M_epochs=2000_20200901 \\
            --dataset /home/ubuntu/opt/MPStats/vignettes/SARS-CoV-2/product/covid19cq1_SARS_TS2_plate_scaled_Cell_MasterDataTable.parquet \\
            --tag UMAP_embedding_into_TS2_2M_epochs=2000_20200901 \\
            --feature_columns ~/opt/MPStats/vignettes/SARS-CoV-2/raw_data/cell_feature_columns_TS_202008.tsv \\
            --no_standardize_features \\
            --re_embed_batch_size 700000 \\
	    --verbose
")


#####################
# cluster all cells #
#####################

cell_features <- dplyr::bind_cols(
    arrow::read_parquet(
        file = "product/TS2_plate_scaled_Cell_MasterDataTable.parquet"),
    arrow::read_parquet(
        file = "intermediate_data/UMAP_embedding_into_TS2_2M_epochs=2000_20200901/umap_embedding.parquet"))

###n_cells	k	runtime min	n_clusters	Total Memory(Gb)
###500000	15	3		39		19
###500000	200	18		5		29.5
###2000000	15	10		147		49


all_cds <- MPStats::populate_cds(
    cell_features = cell_features,
    cell_feature_columns = cell_feature_columns,
    cell_metadata_columns = cell_metadata_columns,
    embedding_type = c("UMAP"),
    embedding = cell_features %>% dplyr::select(UMAP_1, UMAP_2),
    verbose = TRUE)

# as resolution gets bigger --> more clusters
all_cds <- all_cds %>%
    monocle3::cluster_cells(
        k = 200,
        resolution = .00001,
        num_iter = 10,
        verbose = TRUE)

all_cds %>% MPStats::serialize_clusters(
    output_fname = "intermediate_data/UMAP_embedding_TS2_2M_epochs=2000_20200901/clusters_leiden_k=15_res=1e-5.parquet")


#################################
## Re-embed only infected cells #
#################################

# defined in jupyter notebook notebooks/pseudo_time_202008.pynb
infected_cell_membership <- arrow::read_parquet(
    file = "~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/UMAP_embedding_into_TS2_2M_epochs=2000_20200901/infected_ROI_membreship.parquet")

infected_cell_features <- cell_features %>%
    dplyr::bind_cols(infected_cell_membership) %>%
    dplyr::filter(roi_0)

infected_cell_features %>% arrow::write_parquet(
   sink = "product/covid19cq1_SARS_TS2_plate_scaled_infected_Cell_MasterDataTable.parquet")


system("
        cd ~/opt/MPLearn/vignettes/SARS-CoV-2/S25 &&
        /home/ubuntu/anaconda3/envs/sextonlab/bin/python \\
            ~/anaconda3/envs/sextonlab/bin/embed_umap \\
            --dataset /home/ubuntu/opt/MPStats/vignettes/SARS-CoV-2/product/covid19cq1_SARS_TS2_plate_scaled_infected_Cell_MasterDataTable.parquet \\
            --tag UMAP_embedding_TS2_infected_epochs=2000 \\
            --feature_columns ~/opt/MPStats/vignettes/SARS-CoV-2/raw_data/cell_feature_columns_TS_202008.tsv \\
            --no_standardize_features \\
            --umap_n_epochs 2000 \\
	    --verbose
")


####################################################
# Use Monocle3 to cluster and estimate pseudo-time #
####################################################
infected_cell_features <- dplyr::bind_cols(
    arrow::read_parquet(
        "product/covid19cq1_SARS_TS2_plate_scaled_infected_Cell_MasterDataTable.parquet"),
    arrow::read_parquet(
        "~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/UMAP_embedding_TS2_infected_epochs=2000/umap_embedding.parquet"))


infected_cds <- MPStats::populate_cds(
    cell_features = infected_cell_features,
    cell_feature_columns = cell_feature_columns,
    cell_metadata_columns = cell_metadata_columns,
    embedding_type = c("UMAP"),
    embedding = infected_cell_features %>% dplyr::select(UMAP_1, UMAP_2),
    verbose = TRUE)

# as resolution gets bigger --> more clusters
infected_cds <- infected_cds %>%
    monocle3::cluster_cells(
        k = 200,
        resolution = .00001,
        num_iter = 10,
        verbose = TRUE)

infected_cds %>% MPStats::serialize_clusters(
    output_fname = "~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/UMAP_embedding_TS2_infected_epochs=2000/clusters_leiden_k=200_res=1e-5.parquet")

system("ln -s clusters_leiden_k=200_res=1e-5.parquet ~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/UMAP_embedding_TS2_infected_epochs=2000/clusters.parquet")

save(
    infected_cds,
    file = "~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/UMAP_embedding_TS2_infected_epochs=2000/infected_cds.Rdata")


infected_cds <- infected_cds %>%
    monocle3::learn_graph(
        close_loop = FALSE,
        verbose = TRUE)

infected_cds <- infected_cds %>%
    monocle3::order_cells(verbose = TRUE)


save(
    infected_cds,
    file = "~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/UMAP_embedding_TS2_infected_epochs=2000/infected_cds.Rdata")



colData(infected_cds)@listData$time_point <- factor(
    colData(infected_cds)@listData$time_point,
    levels = c("8 hours", "18 hours", "24 hours", "30 hours", "36 hours", "48 hours",  "Uninfected"),
    labels = c("8 hours", "18 hours", "24 hours", "30 hours", "36 hours", "48 hours",  "Uninfected"))

plot <- infected_cds %>%
    monocle3::plot_cells(
        show_trajectory_graph = FALSE,
        label_cell_groups = FALSE,
        cell_size = .2,
        alpha = 1)
ggplot2::ggsave(
    filename = "product/figures/TS2_infected/UMAP_20200903.png",
    plot = plot,
    width = 15,
    height = 15)


plot <- infected_cds %>%
    monocle3::plot_cells(
        show_trajectory_graph = FALSE,
        label_cell_groups = FALSE,
        cell_size = .2,
        alpha = 1) +
    ggplot2::facet_wrap(
        facets = vars(time_point))

ggplot2::ggsave(
    filename = "product/figures/TS2_infected/UMAP_by_time_point_20200903.png",
    plot = plot,
    width = 15,
    height = 15)

###############################################
# plot clusters distribution by drug and dose #
###############################################

figures_output_path <- "product/figures/TS2_infected"

infected_cell_features <- dplyr::bind_cols(
    arrow::read_parquet(
        "product/covid19cq1_SARS_TS2_plate_scaled_infected_Cell_MasterDataTable.parquet"),
    arrow::read_parquet(
        "~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/UMAP_embedding_TS2_infected_epochs=2000/umap_embedding.parquet"),
    arrow::read_parquet(
        "~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/UMAP_embedding_TS2_infected_epochs=2000/clusters.parquet"))

infected_cell_features <- infected_cell_features %>%
    dplyr::mutate(
        time_point = factor(
            time_point,
            levels = c("Uninfected", "8 hours", "18 hours", "24 hours", "30 hours", "36 hours", "48 hours"),
            labels = c("Uninfected", "8 hours", "18 hours", "24 hours", "30 hours", "36 hours", "48 hours")))

data <- infected_cell_features %>% 
    dplyr::group_by(cluster_label) %>%
        dplyr::mutate(cluster_count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(time_point, cluster_label) %>%
        dplyr::mutate(group_count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(
        plate_id,
        time_point,
        cluster_label,
        cluster_count,
        group_count) %>%
    dplyr::mutate(
        is_control = as.character(time_point) == "Uninfected")


data_per_plate <- infected_cell_features %>% 
    dplyr::group_by(plate_id, cluster_label) %>%
        dplyr::mutate(cluster_count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(plate_id, time_point, cluster_label) %>%
        dplyr::mutate(group_count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(
        plate_id,
        time_point,
        cluster_label,
        cluster_count,
        group_count) %>%
    dplyr::mutate(
        is_control = as.character(time_point) == "Uninfected")


plot <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = c(.1, .75)) +
    ggplot2::geom_jitter(
        data = data_per_plate,
        mapping = ggplot2::aes(
            x = time_point,
            y = group_count / cluster_count,
            color = cluster_label),
        size = 1,
        width = .05,
        height = 0) +
    ggplot2::geom_line(
        data = data,
        mapping = ggplot2::aes(
            x = time_point,
            y = group_count / cluster_count,
            group = paste0(cluster_label, is_control),
            color = cluster_label),
        size = .8) +
    ggplot2::geom_point(
        data = data,
        mapping = ggplot2::aes(
            x = time_point,
            y = group_count / cluster_count,
            color = cluster_label),
        size = 1.5) +
    ggplot2::scale_color_discrete("Cluster") +
    ggplot2::scale_x_discrete("Time Point") +
    ggplot2::scale_y_continuous("Fraction of Cells in Cluster")


ggplot2::ggsave(
    filename = paste0(
        figures_output_path,
        "/cluster_by_time_point_",
        MPStats::date_code(),
        ".pdf"),
    plot = plot,
    width = 6,
    height = 4)




data2 <- data %>%
    dplyr::filter(!(plate_id %in% c("1999B", "2020A"))) %>%
    dplyr::filter(!(Compound %in% c("Metformin", "NC", "PC", "BLANK")))

plot <- ggplot2::ggplot(data = data) +
    ggplot2::theme_bw() +
    ggplot2::geom_line(
        mapping = ggplot2::aes(
            x = dose_nM,
            y = group_count / cluster_count,
            color = cluster_label),
        size = .8) +
    ggplot2::geom_point(
        mapping = ggplot2::aes(
            x = dose_nM,
            y = group_count / cluster_count,
            color = cluster_label),
        size = 1) +
    ggplot2::facet_wrap(facets = dplyr::vars(Compound)) +
    list(plate_20XX_scale_x_dose_nM) +
    ggplot2::scale_y_continuous("Fraction of Cells in Cluster")


ggplot2::ggsave(
    filename = paste0(
        figures_output_path,
        "/clusters_by_drug_dose_",
        MPStats::date_code(),
        ".pdf"),
    plot = plot,
    width = 12,
    height = 12)


data2 %>% plyr::ddply(c("cluster_label"), function(df){
    cluster_label <- as.character(df$cluster_label[1])
    cat("Making plot for cluster '", cluster_label, "'\n", sep = "")
    plot <- ggplot2::ggplot(data = df) +
        ggplot2::theme_bw() +
        ggplot2::geom_line(
            mapping = ggplot2::aes(
                x = dose_nM,
                y = group_count / cluster_count),
            size = .8) +
        ggplot2::geom_point(
            mapping = ggplot2::aes(
                x = dose_nM,
                y = group_count / cluster_count),
            size = 1) +
        ggplot2::ggtitle(paste0("Infected cluster '", cluster_label, "'", sep = "")) +
        ggplot2::facet_wrap(facets = dplyr::vars(Compound)) +
        list(plate_20XX_scale_x_dose_nM) +
        ggplot2::scale_y_continuous("Fraction of Cells in Cluster")
    ggplot2::ggsave(
        filename = paste0("product/figures/20XX_infected/cluster=", cluster_label, "_by_drug_dose_20200802.pdf"),
        plot = plot,
        width = 12,
        height = 12)
})


data2 %>% plyr::ddply(c("Compound"), function(df){
    Compound <- as.character(df$Compound[1])
    cat("Making plot for compound '", Compound, "'\n", sep = "")
    plot <- ggplot2::ggplot(data = df) +
        ggplot2::theme_bw() +
        ggplot2::geom_line(
            mapping = ggplot2::aes(
                x = dose_nM,
                y = group_count / cluster_count),
            size = .8) +
        ggplot2::geom_point(
            mapping = ggplot2::aes(
                x = dose_nM,
                y = group_count / cluster_count),
            size = 1) +
        ggplot2::ggtitle(paste0("Infected clusters for compound '", Compound, "'", sep = "")) +
        ggplot2::facet_wrap(facets = dplyr::vars(cluster_label)) +
        list(plate_20XX_scale_x_dose_nM) +
        ggplot2::scale_y_continuous("Fraction of Cells in Cluster")
    ggplot2::ggsave(
        filename = paste0(
            figures_output_path,
            "/Compound=", Compound, "_by_drug_dose_",
            MPStats::date_code(),
            ".pdf"),
        plot = plot,
        width = 12,
        height = 12)
})

################################
# plot specific clusters/drugs #
################################

data <- infected_cell_features %>%
    dplyr::filter(!(plate_id %in% c("1999B", "2020A"))) %>%
    dplyr::filter(!(Compound %in% c("Metformin", "PC", "BLANK"))) %>%
    dplyr::mutate(
        cluster_label = dplyr::case_when(
            cluster_label %in% c(1, 3, 4, 6) ~ "cluster_1,3,4,6",
            cluster_label %in% c(2, 5, 8, 9) ~ "cluster_2,5,8,9",
            TRUE ~ "other")) %>%
    dplyr::group_by(cluster_label) %>%
        dplyr::mutate(cluster_count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Compound, dose_nM, cluster_label) %>%
        dplyr::mutate(group_count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(
        plate_id,
        Compound, dose_nM,
        cluster_label,
        cluster_count,
        group_count)

plot <- ggplot2::ggplot(data = data %>% dplyr::filter(Compound != "NC")) +
    ggplot2::theme_bw() +
    ggplot2::geom_line(
        mapping = ggplot2::aes(
            x = dose_nM,
            y = group_count / cluster_count,
            group = cluster_label,
            color = cluster_label),
        size = .8) +
    ggplot2::geom_point(
        mapping = ggplot2::aes(
            x = dose_nM,
            y = group_count / cluster_count,
            color = cluster_label),
        size = 1) +
    ggplot2::ggtitle("Infected coarse clusters") +
    ggplot2::facet_wrap(facets = dplyr::vars(Compound)) +
    list(plate_20XX_scale_x_dose_nM) +
    ggplot2::scale_y_continuous("Fraction of Cells in Cluster")
ggplot2::ggsave(
    filename = paste0("product/figures/20XX_infected/coarse_clusters_by_drug_dose_20200803.pdf"),
    plot = plot,
    width = 20,
    height = 15)

####

data_20XX <- infected_cell_features %>%
    dplyr::filter(!(plate_id %in% c("1999B", "2020A"))) %>%
    dplyr::filter(!(Compound %in% c("Metformin", "PC", "BLANK"))) %>%
    dplyr::mutate(
        cluster_label = dplyr::case_when(
            cluster_label %in% c(11, 6, 8, 9) ~ as.character(cluster_label),
            TRUE ~ "other")) %>%
    dplyr::group_by(cluster_label) %>%
        dplyr::mutate(cluster_count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Compound, dose_nM, cluster_label) %>%
        dplyr::mutate(group_count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(
        plate_id,
        Compound, dose_nM,
        cluster_label,
        cluster_count,
        group_count)

data_lf_rem <- infected_cell_features %>%
    dplyr::filter(
        plate_id %in% c("1999B", "2020A")) %>%
    dplyr::mutate(
        Remdesivir_Concentration = signif(Remdesivir_Concentration, 3),
        Lactoferin_Concentration = signif(Lactoferrin_Concentration, 3))

data_lf_rem <- dplyr::bind_rows(
    data_lf_rem %>%
    dplyr::filter(
        !is.na(Lactoferrin_Concentration),
        Compound == "Lactoferrin",
        is.na(Remdesivir_Concentration) | Remdesivir_Concentration == 0) %>%
    dplyr::mutate(
        dose_nM = Lactoferrin_Concentration),
    data_lf_rem %>%
    dplyr::filter(
        !is.na(Remdesivir_Concentration),
        Compound == "Remdesivir",
        is.na(Lactoferrin_Concentration) | Lactoferrin_Concentration == 0) %>%
    dplyr::mutate(
        dose_nM = Remdesivir_Concentration * 10^3)) %>%
    dplyr::mutate(
        cluster_label = dplyr::case_when(
            cluster_label %in% c(11, 6, 8, 9) ~ as.character(cluster_label),
            TRUE ~ "other")) %>%
    dplyr::filter(cluster_label != "other") %>%
    dplyr::group_by(cluster_label) %>%
        dplyr::mutate(cluster_count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Compound, dose_nM, cluster_label) %>%
        dplyr::mutate(group_count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(
        plate_id,
        Compound, dose_nM,
        cluster_label,
        cluster_count,
        group_count)

data <- dplyr::bind_rows(
    data_20XX,
    data_lf_rem)


plot <- ggplot2::ggplot(data = data %>% dplyr::filter(Compound != "NC")) +
    ggplot2::theme_bw() +
    ggplot2::geom_line(
        mapping = ggplot2::aes(
            x = dose_nM,
            y = group_count / cluster_count,
            group = cluster_label,
            color = cluster_label),
        size = .8) +
    ggplot2::geom_point(
        mapping = ggplot2::aes(
            x = dose_nM,
            y = group_count / cluster_count,
            color = cluster_label),
        size = 1) +
    ggplot2::ggtitle("Infected coarse clusters") +
    ggplot2::facet_wrap(facets = dplyr::vars(Compound)) +
    list(plate_20XX_scale_x_dose_nM) +
    ggplot2::scale_y_continuous("Fraction of Cells in Cluster")
ggplot2::ggsave(
    filename = paste0("product/figures/20XX_infected/infected_clusters_actually_by_drug_dose_20200803.pdf"),
    plot = plot,
    width = 20,
    height = 15)


######

plot <- ggplot2::ggplot(
    data = data %>%
        dplyr::filter(
            cluster_label %in% c("11", "6", "8", "9"),
            Compound %in% c("S1RA", "Thioguanine", "Z-FA-FMK", "Gilteritinib", "Lactoferrin", "Remdesivir")) %>%
        dplyr::mutate(Compound = factor(
            Compound,
            levels = c("Thioguanine", "S1RA", "Lactoferrin", "Gilteritinib", "Z-FA-FMK", "Remdesivir"),
            labels = c("Thioguanine", "S1RA", "Lactoferrin", "Gilteritinib", "Z-FA-FMK", "Remdesivir"))) %>%
        dplyr::filter(plate_id != "1999B")) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::geom_line(
        mapping = ggplot2::aes(
            x = dose_nM,
            y = group_count,
            group = cluster_label,
            color = cluster_label),
        size = .8) +
    ggplot2::geom_point(
        mapping = ggplot2::aes(
            x = dose_nM,
            y = group_count,
            color = cluster_label),
        size = 1) +
    ggplot2::ggtitle("Infected coarse clusters") +
    ggplot2::facet_wrap(
        facets = dplyr::vars(Compound),
        scales = "free") +
    ggplot2::scale_x_log10("Dose nM") +
    ggplot2::scale_y_continuous("Cells in Cluster")
ggplot2::ggsave(
    filename = paste0("product/figures/20XX_infected/infected_clusters_8,9_by_drug_dose_20200826.pdf"),
    plot = plot,
    width = 5,
    height = 4)



plot <- ggplot2::ggplot(
    data = dplyr::bind_row(
        data %>%
        dplyr::filter(
            cluster_label %in% c("8", "9"),
            Compound %in% c("S1RA", "Thioguanine", "Z-FA-FMK", "Gilteritinib")) %>%
        dplyr::select(
            Compound,
            dose_nM,
            group_count,
            cluster_label)
        data %>%
        dplyr::filter(
                



    +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::geom_line(
        mapping = ggplot2::aes(
            x = dose_nM,
            y = group_count,
            group = cluster_label,
            color = cluster_label),
        size = .8) +
    ggplot2::geom_point(
        mapping = ggplot2::aes(
            x = dose_nM,
            y = group_count,
            color = cluster_label),
        size = 1) +
    ggplot2::ggtitle("Infected coarse clusters") +
    ggplot2::facet_wrap(
        facets = dplyr::vars(Compound),
        scales = "free") +
    list(plate_20XX_scale_x_dose_nM) +
    ggplot2::scale_y_continuous("Cells in Cluster")
ggplot2::ggsave(
    filename = paste0("product/figures/20XX_infected/infected_clusters_8,9_by_drug_dose_20200803.pdf"),
    plot = plot,
    width = 5,
    height = 4)


############################################
# clusters by lactoferrin remdesivir plate #
############################################

# 1999B
# Lf_conc    R0 R0.0032 R0.0062 R0.0124 R0.025 R0.048  R0.1  R0.2
#   0.000 0.616   0.681   0.754   0.742  0.749  0.720 0.746 0.762
#   0.390 0.577   0.649   0.661   0.733  0.681  0.677 0.717 0.738
#   0.780 0.564   0.560   0.618   0.671  0.693  0.694 0.712 0.701
#   1.562 0.608   0.634   0.678   0.712  0.601  0.751 0.661 0.734
#   3.125 0.605   0.612   0.615   0.633  0.640  0.727 0.732 0.725
#   6.250 0.608   0.747   0.732   0.721  0.703  0.738 0.687 0.692
#  12.500 0.675   0.709   0.751   0.735  0.705  0.762 0.747 0.769
#  25.000 0.709   0.803   0.770   0.803  0.730  0.784 0.786 0.857

# 2020A
#   Lf_Conc    R0 R2e.04 R4e.04 R8e.04 R0.0016 R0.0032 R0.0062 R0.0124 R0.025
# 1       0 0.994  0.636  0.186  0.397      NA   0.657   0.739   0.820  0.746
# 2      25 0.747  0.662     NA  0.780      NA      NA      NA      NA     NA
# 3      50 0.636     NA  0.701     NA      NA      NA      NA      NA     NA
# 4      75 0.770     NA  0.793  0.704   0.783      NA      NA      NA     NA
# 5     100 0.705     NA     NA     NA   0.627   0.789      NA      NA     NA
# 6     125 0.771     NA     NA     NA      NA   0.722   0.588      NA     NA
# 7     150 0.765     NA     NA     NA      NA      NA   0.731   0.705     NA
# 8     175 0.610     NA     NA     NA      NA      NA      NA   0.759  0.791
# 9     200 0.816     NA     NA     NA      NA      NA      NA      NA  0.768

lactoferrin_scale_x_continuous <- list(
    ggplot2::theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust=.5)),
    ggplot2::scale_x_log10(
        "Lactoferrin Dose (µg/mL)",
        breaks=c(0, 0.39, 0.78, 1.56, 3.12, 6.25, 12.5, 25, 50, 75, 100, 125, 150, 175, 200)*10^-6,
        labels=c("0", "0.39", "0.78", "1.56", "3.12", "6.25", "12.5", "25", "50", "75", "100", "125", "150", "175", "200")))

# single agent 
data <- infected_cell_features %>%
    dplyr::filter(
        plate_id %in% c("1999B", "2020A")) %>%
    dplyr::mutate(
        Remdesivir_Concentration = signif(Remdesivir_Concentration, 3),
        Lactoferin_Concentration = signif(Lactoferrin_Concentration, 3))

data <- dplyr::bind_rows(
    data %>%
    dplyr::filter(
        !is.na(Lactoferrin_Concentration),
        Compound == "Lactoferrin",
        is.na(Remdesivir_Concentration) | Remdesivir_Concentration == 0) %>%
    dplyr::mutate(
        dose_nM = Lactoferrin_Concentration * 10^-6),
    data %>%
    dplyr::filter(
        !is.na(Remdesivir_Concentration),
        Compound == "Remdesivir",
        is.na(Lactoferrin_Concentration) | Lactoferrin_Concentration == 0) %>%
    dplyr::mutate(
        dose_nM = Remdesivir_Concentration)) %>%
    dplyr::mutate(
        cluster_label = dplyr::case_when(
            cluster_label %in% c(2, 5, 8, 9) ~ as.character(cluster_label),
            TRUE ~ "other")) %>%
    dplyr::filter(cluster_label != "other") %>%
    dplyr::group_by(cluster_label) %>%
        dplyr::mutate(cluster_count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Compound, dose_nM, cluster_label) %>%
        dplyr::mutate(group_count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(
        plate_id,
        Compound, dose_nM,
        cluster_label,
        cluster_count,
        group_count)
    



    
plot <- ggplot2::ggplot(data = data) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::geom_line(
        mapping = ggplot2::aes(
            x = dose_nM,
            y = group_count,
            group = cluster_label,
            color = cluster_label),
        size = .8) +
    ggplot2::geom_point(
        mapping = ggplot2::aes(
            x = dose_nM,
            y = group_count,
            color = cluster_label),
        size = 1) +
    ggplot2::ggtitle("Infected coarse clusters") +
    ggplot2::facet_grid(
        rows = dplyr::vars(plate_id),
        cols = dplyr::vars(Compound),
        scales = "free") +
#    ggplot2::facet_wrap(
#        facets = dplyr::vars(Compound),
#        scales = "free") +
    ggplot2::scale_x_log10("Drug dose") +
#    list(plate_20XX_scale_x_dose_nM) +
    ggplot2::scale_y_continuous("Cells in Cluster")
ggplot2::ggsave(
    filename = paste0("product/figures/20XX_infected/plates_lf_rem_single_agent_infected_clusters_8,9_by_drug_dose_20200803.pdf"),
    plot = plot,
    width = 5,
    height = 6)


    

##########################################
# PLS enriched features for each cluster #
##########################################

feature_importance <- infected_cell_features %>%
    dplyr::distinct(cluster_label) %>%
    plyr::adply(1, function(df) {
        cat("Computing feature importance for cluster '", df$cluster_label[1], "'\n", sep = "")
        data_matrix <- infected_cell_features  %>%
            dplyr::mutate(label = dplyr::case_when(
                as.character(cluster_label) != df$cluster_label[1] ~ "other",
                as.character(cluster_label) == df$cluster_label[1] ~ "cluster",
                TRUE ~ NA_character_) %>%
                    factor(
                        levels = c("cluster", "other"),
                        labels = c("cluster", "other"))) %>%
            dplyr::select(
                label,
                tidyselect::one_of(cell_feature_columns$feature)) %>%
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
                cell_feature_columns %>%
                dplyr::select(feature)) %>%
            dplyr::rename(feature_score = Overall) %>%
            dplyr::arrange(desc(feature_score)) %>%
            head(20) %>%
            dplyr::mutate(cluster_label = df$cluster_label[1])
    })

p <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_point(
        data = feature_importance %>%
            dplyr::mutate(
                cluster_index = cluster_label %>%
                    stringr::str_extract("[0-9]+$") %>%
                    as.numeric(),
                cluster_label = factor(cluster_label) %>%
                    reorder(cluster_index)),
        mapping = ggplot2::aes(
            x = feature_score,
            y = feature)) +
    ggplot2::scale_x_continuous("PLS feature importance") +
    ggplot2::scale_y_discrete(
        "Feature",
        drop = FALSE) +
    ggplot2::facet_wrap(
        ~cluster_label,
        ncol = 3,
        scales = "free_y")

ggplot2::ggsave(
        filename = paste0(
            figures_output_path,
            "/clusters_feature_importance_",
            MPStats::date_code(),
            ".pdf"),
    plot = p,
    width = 15,
    height = 25)

feature_importance %>%
    readr::write_tsv(
        figures_output_path,
        "/clusters_feature_importance_",
        MPStats::date_code(),
        ".tsv")





##########################################
# Monocle3 based representative features #
##########################################

marker_test_res <- monocle3::top_markers(
    infected_cds,
    cores = 30,
    verbose = TRUE)

top_specific_markers <- marker_test_res %>%
    dplyr::group_by(cell_group) %>%
    dplyr::top_n(1, pseudo_R2)

top_specific_marker_ids <- top_specific_markers %>%
    dplyr::pull(gene_id) %>%
    unique()

plot <- monocle3::plot_genes_by_group(
    infected_cds,
    top_specific_marker_ids,
    ordering_type = "maximal_on_diag",
    max.size = 3)

ggplot2::ggsave(
    filename = "product/figures/pseudo_time_202006/UMAP_embedding_TS_scaled_infected_a=1_b=0.5_n_neighbors=500_negative_sample_rate=20_n_epochs=2000_top_markers_200722.pdf",
    plot = plot,
    width = 12,
    height = 12)

##################################
# Get instances for each cluster #
##################################

source("scripts/database.R")
con <- get_primary_database_connection(schema = "covid19cq1")
con %>% DBI::dbListTables()
con %>% dplyr::tbl('SARS_2009A_Per_Nuclei') %>% dplyr::glimpse()


#####################
# cluster_neighbors #
#####################



########################
cell_features_20XX <- c(
    '2006A', '2007A', '2008A', '2009A',
    '2010A', '2010A',          '2012A',
    '2013A', '2014A', '2015A', '2016A',
    '2017A',          '2019A') %>%
    plyr::ldply(function(plate_id) {
        cat("Reading features from plate ", plate_id, " ...\n", sep="")
        arrow::read_parquet(
            paste0("product/covid19cq1_SARS_", plate_id, "_plate_scaled_Cell_MasterDataTable.parquet"))
    })
cell_features_rem_lf <- arrow::read_parquet(
    "product/lf_rem_plate_scaled_1999B_2020A_Cell_MasterDataTable.parquet")
infected_ROI_membership <- arrow::read_parquet(
    "~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/top_hits_plate_scaled_200522a_umap2_2M_15_0.0/infected_ROI_membreship.parquet")
cell_features <- dplyr::bind_rows(
    cell_features_20XX,
    cell_features_rem_lf) %>%
    dplyr::bind_cols(
        infected_ROI_membership)



infected_cell_features <- dplyr::bind_cols(
    arrow::read_parquet(
        "product/top_hits_infected_plate_scaled_200730_Cell_MasterDataTable.parquet"),
    arrow::read_parquet(
        "~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/UMAP_embedding_top_hits_infected_plate_scaled_epochs=2000_re_embed_epochs=2000_200730/umap_embedding.parquet"),
    arrow::read_parquet(
        "~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/UMAP_embedding_top_hits_infected_plate_scaled_epochs=2000_re_embed_epochs=2000_200730/clusters.parquet"))

cell_features$cluster_label <- -1
cell_features[cell_features$roi_0, ]$cluster_label <- infected_cell_features$cluster_label


source("scripts/microscope_geometry.R")
neighbor_distance_threshold_mm <- 30

infected_cell_neighbors <- cell_features %>%
    dplyr::mutate(
        field_index = as.numeric(Image_Metadata_FieldID),
        well_x = cq1_field_pixel_to_well_mm_x(
            coordinate_x = Nuclei_Location_Center_X,
            field_index = field_index),
        well_y = cq1_field_pixel_to_well_mm_y(
            coordinate_y = Nuclei_Location_Center_Y,
            field_index = field_index)) %>%
    dplyr::select(
        Image_Metadata_PlateID,
        Image_Metadata_WellID,
        well_x,
        well_y,
        cluster_label) %>%
    dplyr::group_by(
        Image_Metadata_PlateID,
        Image_Metadata_WellID) %>%
    do({
        data <- .
        z <- dplyr::select(data, well_x, well_y)
        neighbors <- rdist::cdist(X = z, Y = z, metric = "euclidean") %>%
            as.data.frame.table(responseName = "distance_mm") %>%
            dplyr::mutate_if(is.factor, as.integer) %>%
            dplyr::rename(
                cell_index_1 = Var1,
                cell_index_2 = Var2) %>%
            dplyr::filter(
                distance_mm <= neighbor_distance_threshold_mm,
                cell_index_1 != cell_index_2) %>%
            dplyr::left_join(
                dplyr::transmute(data,
                    cell_index_1 = dplyr::row_number(),
                    cluster_label_1 = cluster_label),
                by = c("cell_index_1")) %>%
            dplyr::left_join(
                 dplyr::transmute(data,
                    cell_index_2 = dplyr::row_number(),
                    cluster_label_2 = cluster_label),
                by = c("cell_index_2"))
    }) %>%
    dplyr::ungroup()



# x1 = 3
# x2 = 10
# ground truth model:
#    count = x1 * x2
toy_data <-
    data.frame(
        x1 = c("1", "1", "2"),
        x2 = c("1", "2", "2"),
        count = c(9, 30, 100))

toy_model <- glm(
    formula = count ~ x1 + x2 + 0,
    data = toy_data,
    family = stats::poisson())

toy_model %>% predict %>% exp



cluster_neighbor_counts <- infected_cell_neighbors %>%
    dplyr::count(
        cluster_label_1,
        cluster_label_2,
        name = "cluster_neighbor_count") %>%
    dplyr::filter(as.numeric(cluster_label_1) <= as.numeric(cluster_label_2)) %>%
    dplyr::mutate(
        cluster_neighbor_count =
            dplyr::case_when(
                cluster_label_1 == cluster_label_2 ~ as.double(cluster_neighbor_count) / 2,
                TRUE ~ as.double(cluster_neighbor_count)))
        
cluster_neighbor_counts <- cluster_neighbor_counts %>%
    dplyr::filter(!(cluster_label_1 == "-1" & cluster_label_2 == "-1"))

model <- glm(
    formula = cluster_neighbor_count ~ cluster_label_1 + cluster_label_2,
    data = cluster_neighbor_counts,
    family = stats::poisson(link = "log"))

cluster_neighbor_counts <- cluster_neighbor_counts %>%
    dplyr::mutate(
        baseline_counts = model %>%
            predict %>%
            exp %>%
            as.numeric) %>%
    dplyr::mutate(
        cluster_neighbor_risk = cluster_neighbor_count / baseline_counts)

cluster_neighbor_counts %>%
    readr::write_tsv("product/figures/20XX_infected/cluster_neighbors/cluster_neighbor_counts_20200811.tsv")

cluster_neighbor_counts %>%
    dplyr::select(
        cluster_label_1,
        cluster_label_2,
        cluster_neighbor_risk) %>%
    dplyr::mutate(
        cluster_neighbor_risk = signif(cluster_neighbor_risk * 100, 3)) %>%
    tidyr::pivot_wider(
        id_cols = cluster_label_1,
        names_from = cluster_label_2,
        values_from = cluster_neighbor_risk) %>%
    data.frame    
    

    
# https://kateto.net/network-visualization
cluster_neighbor_network <- igraph::graph_from_data_frame(
  d = cluster_neighbor_counts,
  directed = FALSE) %>%
  igraph::simplify(
    remove.multiple = FALSE,
    remove.loops = TRUE)

igraph::E(cluster_neighbor_network)$width <- igraph::E(cluster_neighbor_network)$cluster_neighbor_risk * 5

igraph::E(cluster_neighbor_network)$weight <- igraph::E(cluster_neighbor_network)$cluster_neighbor_risk * 1000

pdf("product/20XX_infected/cluster_neighbors/cluster_neighbor_graph_20200811.pdf")
cluster_neighbor_network %>% plot(
  vertex.label.font = 1,
  layout = cluster_neighbor_network %>%
    igraph::layout_with_drl(
      options = list(
        simmer.attraction = 1.5,
        init.temperature = 10,
        init.iterations = 10000)))
def.off()

#################################################################
#
#################################################################
