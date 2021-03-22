

library(plyr)
library(dplyr)
library(monocle3)
library(MPStats)


# off chip control
load("intermediate_data/cds_CZ.Rdata")
plot <- cds_CZ %>% monocle3::plot_cells()
ggplot2::ggsave(filename = "product/figures/CZ/umap_20210208.pdf")


# on chip control
load("intermediate_data/cds_CZ_2.Rdata")
plot <- cds_CZ_2 %>% monocle3::plot_cells()
ggplot2::ggsave(filename = "product/figures/CZ_2/umap_20210208.pdf")

# compare off and on chip
load("intermediate_data/cds_CZ_CZ_1.Rdata")

plot <- cds_CZ_CZ_1 %>% monocle3::plot_cells()
ggplot2::ggsave(
    filename = "product/figures/CZ_CZ_1/umap_aligned_k=200_clusters_5e-4_20210208.pdf")

plot <- cds_CZ_CZ_1 %>% monocle3::plot_cells(color_cells_by = "sample")
ggplot2::ggsave(
    filename = "product/figures/CZ_CZ_1/umap_aligned_k=200_20210208.pdf")


# all on-chip samples
load("intermediate_data/cds_CZ_x.Rdata")

plot <- cds_CZ_x %>%
    monocle3::plot_cells(color_cells_by = "sample")
ggplot2::ggsave(
    filename = "product/figures/CZ_x/umap_aligned_k=200_1e-5_20210208.pdf")

plot <- cds_CZ_x %>%
    monocle3::plot_cells()
ggplot2::ggsave(
    filename = "product/figures/CZ_x/umap_aligned_k=200_clusters_1e-5_20210208.pdf")

# Ouchi2019
load("intermediate_data/cds_Ouchi2019.Rdata")
plot <- cds_Ouchi2019 %>%
    monocle3::plot_cells(color_cells_by = "cluster")
ggplot2::ggsave(
    filename = "product/figures/Ouchi2019/umap_a=10,b=1_1e-4_20210208.pdf")

# Ouchi2019 + on_chip control
load("intermediate_data/cds_Ouchi2019_CZ_2.Rdata")
plot <- cds_Ouchi2019_CZ_2 %>%
    monocle3::plot_cells(
        color_cells_by = "sample")
ggplot2::ggsave(
    filename = "product/figures/Ouchi2019_CZ_2/umap_aligned_k=200_1e-5_20210208.pdf")

plot <- cds_Ouchi2019_CZ_2 %>%
    monocle3::plot_cells()
ggplot2::ggsave(
    filename = "product/figures/Ouchi2019_CZ_2/umap_aligned_k=200_clusters_1e-5_20210208.pdf")

# Ouchi2019 + all on chip samples
plot <- cds_Ouchi2019_CZ_x %>%
    monocle3::plot_cells(
        color_cells_by = "sample")
ggplot2::ggsave(
    filename = "product/figures/Ouchi2019_CZ_x/umap_aligned_k=200_1e-5_20210208.pdf")

plot <- cds_Ouchi2019_CZ_x %>%
    monocle3::plot_cells()
ggplot2::ggsave(
    filename = "product/figures/Ouchi2019_CZ_x/umap_aligned_k=200_clusters_1e-5_20210208.pdf")



z <- cds_Ouchi2019_CZ_x[, cds_Ouchi2019_CZ_x$sample  == "Ouchi2019"]
plot <- z %>%
    monocle3::plot_cells()
ggplot2::ggsave(
    filename = "product/figures/Ouchi2019_CZ_x/umap_Ouchi2019_only_aligned_k=200_clusters_1e-5_20210208.pdf")

