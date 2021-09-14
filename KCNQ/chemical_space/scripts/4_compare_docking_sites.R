


library(tidyverse)
library(ggrepel)

data <- dplyr::bind_rows(
    readr::read_tsv(
        "intermediate_data/KCNQ2_7CR2_retigabine_AB_20201028,project_20210323,,20210402/pose_features.tsv") %>%
    dplyr::mutate(
        site = "retigabine",
        database = "project"),
    readr::read_tsv(
        "intermediate_data/KCNQ2_7CR2_retigabine_AB_20201028,project_decoys_20210323,,20210402/pose_features.tsv") %>%
    dplyr::mutate(
        site = "retigabine",
        database = "decoys"),
    readr::read_tsv(
        "intermediate_data/KCNQ2_7CR1_ztz240_A20210419,project_20210323,,20210519/pose_features.tsv") %>%
    dplyr::mutate(
        site = "ztz240",
        database = "project"),
    readr::read_tsv(
        "intermediate_data/KCNQ2_7CR1_ztz240_A20210419,project_decoys_20210323,,20210519/pose_features.tsv") %>%
    dplyr::mutate(
        site = "ztz240",
        database = "decoys"))



plot_data <- data %>%
    dplyr::select(
        site,
        database,
        name,
        pose_id,
        total_energy) %>%
    tidyr::pivot_wider(
        id_cols = c(database, name, pose_id),
        names_from = site,
        values_from = total_energy,
        values_fill = 0)

plot <- ggplot2::ggplot(data) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(
        label = "Docking site binding energy",
        subtitle = "Retigabine vs. ztz240 sites") +
    ggplot2::scale_x_continuous(
        "Retigabine Site (Dock Energy)") +
    ggplot2::scale_y_continuous(
        "Ztz240  Site (Dock Energy)") +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "grey30", size = 1.4) +
    ggplot2::geom_point(
        data = plot_data %>% dplyr::filter(database == "decoys"),
        mapping = ggplot2::aes(
            x = retigabine,
            y = ztz240,
            color = database),
        alpha = .7,
        color = "grey70") +
    ggplot2::geom_point(
        data = plot_data %>% dplyr::filter(database == "project"),
        mapping = ggplot2::aes(
            x = retigabine,
            y = ztz240,
            color = database),
        alpha = .8,
        color = "#000080") +
    ggrepel::geom_label_repel(
        data = plot_data %>% dplyr::filter(
           database == "project",
           (retigabine < -33) | (ztz240 < -33)),
        mapping = ggplot2::aes(
            x = retigabine,
            y = ztz240,
            label = name),
        size = 1.5,
        min.segment.length = 0,
        force = 10) +
    ggplot2::coord_equal()

ggplot2::ggsave(
    filename = "product/figures/compare_docking_sites/retigabine_ztz240_project_decoys.pdf",
    width = 8,
    height = 7,
    useDingbats = FALSE)
