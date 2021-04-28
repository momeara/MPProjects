
library(plyr)
library(dplyr)
library(ggrepel)

source("parameters.R")


project_activities <- readr::read_tsv(
    file = "raw_data/activities_20210323.tsv")

project_substances <- readr::read_tsv(
    file = "raw_data/substances_20210323.tsv") %>%
    dplyr::select(
        substance_source,
        substance_name,
        substance_dock_id) %>%
    dplyr::mutate(
        substance_source = ifelse(substance_name == "SCR-2682", "", substance_source)) %>%
    dplyr::mutate(substance_dock_id = substance_dock_id %>% stringr::str_sub(-14, -1))

project_substances_docked_features <- readr::read_tsv(
    file = "intermediate_data/KCNQ2_7CR2_retigabine_AB_20201028,project_20210323,,20210402/docking_features.tsv") %>%
    dplyr::mutate(substance_dock_id = name %>% stringr::str_replace("[.][0-9]+$", "")) %>%
    dplyr::mutate(score_tranche = dplyr::case_when(
        total_energy < -28 ~ "Top",
        total_energy < -10 ~ "Mid",
        TRUE ~ "Bad"))

project_sar <- project_activities %>%
    dplyr::left_join(
        project_substances %>% dplyr::select(substance_name, substance_dock_id),
        by = "substance_name") %>%
    dplyr::left_join(
        project_substances_docked_features,
        by = "substance_dock_id") 


project_sar_plot_data <- project_sar %>%
    dplyr::filter(receptor %in% c("KCNQ2", "KCNQ2/3")) %>%
    dplyr::select(
        reference,
        substance_name,
        total_energy,
        value)    
dock_energy_offset <- 50
save_plot <- TRUE
project_sar %>%
    dplyr::group_by(reference) %>%
    dplyr::do({
        data <- .
        cat("Plotting SAR for reference: '", data$reference[1], "' with ", nrow(data), " substances at\n", sep = "")
        plot <- ggplot2::ggplot() +
            ggplot2::theme_bw() +
            ggplot2::geom_point(
                data = data,
                mapping = ggplot2::aes(
                    x = log10(total_energy + dock_energy_offset),
                    y = as.numeric(value))) +
            ggrepel::geom_text_repel(
                data = data,
                mapping = ggplot2::aes(
                    x = log10(total_energy + dock_energy_offset),
                    y = as.numeric(value),
                    label = substance_name),
                force = 10) +
            ggplot2::ggtitle(
                label = "Structure activity relationship at the Retigabine binding site from 7CR2",
                subtitle = paste0("Study: ", data$reference[1]))
        if (save_plot) {
            fname <- paste0(
                "product/figures/retigabine_site_SAR/",
                "SAR_", gsub("[^a-zA-Z0-9,.]", "", data$reference[1]), ".pdf")
            cat("saving plot '", fname, "' ...\n", sep = "")
            ggplot2::ggsave(
                filename = fname,
                width = 8,
                height = 8,
                useDingbats = FALSE)
            fname <- paste0(
                "product/figures/retigabine_site_SAR/",
                "SAR_", gsub("[^a-zA-Z0-9,.]", "", data$reference[1]), ".png")
            cat("saving plot '", fname, "' ...\n", sep = "")
            ggplot2::ggsave(
                filename = fname,
                width = 8,
                height = 8)
        }
        data.frame()
    })        


dock_energy_offset <- 50
project_sar_plot_data <- project_sar %>%
    dplyr::filter(reference %in% c(
        "(Bock, 2018, 10.1039/C8OB02530D)",
        "(Surur, 2019, 10.1002/open.201800244)",
        "(Surur, 2019, 10.1039/C9OB00511K)")) %>%
    dplyr::select(
        reference,
        substance_name,
        total_energy,
        value) %>%
    dplyr::filter(!is.na(total_energy)) %>%
    dplyr::mutate(
        experimental_EC50 = dplyr::case_when(
            value == ">10" ~ 100,
            TRUE ~ as.numeric(value)),
        dock_energy_transformed = dplyr::case_when(
            is.na(total_energy) ~ log10(0 + dock_energy_offset),
            TRUE ~ log10(total_energy + dock_energy_offset)))

save_plot <- TRUE
plot <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_point(
        data = project_sar_plot_data,        
        mapping = ggplot2::aes(
            x = dock_energy_transformed,
            y = log10(experimental_EC50),
            color = reference)) +
    ggrepel::geom_text_repel(
        data = project_sar_plot_data %>%
            dplyr::filter((value == ">10") |  (!is.na(total_energy))),
        mapping = ggplot2::aes(
            x = dock_energy_transformed,
            y = log10(experimental_EC50),
            label = substance_name),
        size = 3,
        force = 10) +
    ggplot2::ggtitle(
        label = "Structure activity relationship at the Retigabine binding site from 7CR2",
        subtitle = paste0("Link Lab KCNQ2/3")) +
    ggplot2::scale_y_continuous(
        "Thallium Flux EC50 (uM)",
        limits = log10(c(.003, 100)),
        breaks = log10(c(.003, .01, .03, .1, .3, 1, 3, 10, 30, 100)),
        labels = c(.003, .01, .03, .1, .3, 1, 3, 10, 30, 100)) +
    ggplot2::scale_x_continuous(
        "Dock Score",
        limits = log10(c(-36, -14) + dock_energy_offset),
        breaks = log10(c(-35, -30, -25, -20, -15) + dock_energy_offset),
        labels = c(-35, -30, -25, -20, -15)) +
    ggplot2::theme(legend.position = "bottom")

if (save_plot) {
    fname <- paste0("product/figures/retigabine_site_SAR/SAR_Link_lab.pdf")
    cat("saving plot '", fname, "' ...\n", sep = "")
    ggplot2::ggsave(
        filename = fname,
        width = 8,
        height = 8,
        useDingbats = FALSE)
    fname <- paste0("product/figures/retigabine_site_SAR/SAR_Link_lab.png")
    cat("saving plot '", fname, "' ...\n", sep = "")
    ggplot2::ggsave(
        filename = fname,
        width = 8,
        height = 8)
}
