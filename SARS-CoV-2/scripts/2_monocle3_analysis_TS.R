library(plyr)
library(tidyverse)
library(arrow)
library(monocle3)
library(BiocNeighbors)
library(MPStats)

# gathered in 3_pseudo_time_202006_UMAP_embedding
cell_features <- arrow::read_parquet(
    "raw_data/covid19cq1_TS_scaled_Cell_MasterDataTable.parquet")

embedding <- arrow::read_parquet(
    "~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/UMAP_embedding_TS_scaled_full/umap_embedding.parquet")

cds <- MPStats::populate_cds(
    cell_features = cell_features,
    cell_feature_columns = cell_feature_columns,
    cell_metadata_columns = cell_metadata_columns,
    embedding = embedding)

options(shiny.port = 4999)
infected_cells <- cds %>% monocle3::choose_cells()

save(
    infected_cells,
    file = "~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/UMAP_embedding_TS_scaled_full/infected_cell_csd.Rdata")

infected_cells <- infected_cells %>% monocle3::reduce_dimension(
    verbose=TRUE)


cell_features <- cell_features %>%
    dplyr::mutate(cell_index = row_number()) %>%
    dplyr::left_join(
        data.frame(
            cell_index = colData(infected_cells)@rownames %>% as.numeric(),
            is_infected = TRUE),
        by = "cell_index") %>%
    dplyr::mutate(
        is_infected = dplyr::case_when(
            is.na(is_infected) ~ FALSE,
            TRUE ~ TRUE))

cell_features %>%
    arrow::write_parquet(
       "raw_data/covid19cq1_TS_scaled_Cell_MasterDataTable.parquet")

    


library(spatstat)
library(spdep)

image_width <- 2560
image_height <- 2160

# #layout
# #
# #0   1   2
# #
# #3   4   5
# #
# #6   7   8
# ################

# exclude wells with no infecteted cells for now
point_processes <- cell_features %>%
    dplyr::group_by(plate_id, Compound, Image_Metadata_WellID) %>%
    do({
        cat(
            "Generating point process for ",
            "plate '", as.character(.$plate_id[1]), "', ",
            "compound '", as.character(.$Compound[1]), "', ",
            "well '", as.numeric(.$Image_Metadata_WellID[1]), "'\n",
            sep = "")
        well_data <- dplyr::mutate(., 
                well_x = (as.numeric(Image_Metadata_FieldID) - 1) %% 3 * image_width +
                    Nuclei_Location_Center_X,
                well_y = floor((as.numeric(Image_Metadata_FieldID) - 1) / 3) * image_height +
                    Nuclei_Location_Center_Y)
        if (sum(well_data$is_infected) == 0) {
            cat(" No infected cells!\n")
            data.frame()
        } else {
            point_process <- spatstat::ppp(
                x = well_data$well_x,
                y = well_data$well_y,
                window = spatstat::owin(
                    xrange = c(0, 3 * image_width),
                    yrange = c(0, 3 * image_height)),
                marks = well_data %>%
                    dplyr::transmute(
                        is_infected = factor(
                            is_infected,
                            levels = c(TRUE, FALSE),
                            labels = c("Infected", "Not Infected"))))
            tibble::tibble(point_process = list(point_process))
        }
    }) %>%
    dplyr::ungroup()

pcf_infected_not_infected <- point_processes %>%
    dplyr::group_by(plate_id, Compound) %>%
    dplyr::do({
        cat(
            "Computing Kcross for plate '", as.character(.$plate_id[1]), "' ",
            "Compound '", as.character(.$Compound[1]), "'\n", sep = "")
        treatment_point_processes <- .
        expand.grid(
            from = c("Infected", "Not Infected"),
            to = c("Infected", "Not Infected")) %>%
            plyr::adply(1, function(cell_pair_types) {
                cat(
                    "  Compounting Kcross ",
                    "'", as.character(cell_pair_types$from[1]), "' -> ",
                    "'", as.character(cell_pair_types$to[1]), "'\n",
                    sep = "")
                tryCatch({
                    data.frame(
                        Kcross = lapply(
                            treatment_point_processes$point_process,
                            function(point_process) {
                                point_process %>%
                                    spatstat::Kcross(
                                        from = cell_pair_types$from[1],
                                        to = cell_pair_types$to[1]) %>%
                                    spatstat::pcf(
                                        spar = 0.8,
                                        method = "b")}) %>%
                            spatstat::as.anylist() %>%
                            spatstat::pool() %>%
                            list())
                    },
                    error = function(e) {
                        cat("    ERROR:", e$message, "\n")
                        data.frame()
                    })
            })
    }) %>%
    dplyr::ungroup()

data <- pcf_infected_not_infected %>%
    dplyr::select(
        -Kcross.pooltheo,
        -Kcross.vartheo,
        -Kcross.hitheo,
        -Kcross.lotheo) %>%
    tidyr::pivot_wider(
        id_cols = c(plate_id, Compound, Kcross.r),
        names_from = c(from, to),
        values_from = c(
            Kcross.poolpcf,
            Kcross.varpcf,
            Kcross.hipcf,
            Kcross.lopcf))

plot <- ggplot2::ggplot(
    data = pcf_infected_not_infected %>%
        dplyr::filter(
            from == "Infected",
            to == "Infected")) +
    ggplot2::theme_bw() +
    ggplot2::geom_hline(
        yintercept = 1) +
    ggplot2::geom_line(
        mapping = ggplot2::aes(
            x = Kcross.r,
            y = Kcross.poolpcf,
            color = paste0(from, " -> ", to))) +
    ggplot2::scale_y_continuous(
        "Pair Correlation Function") +
    ggplot2::scale_x_continuous(
        "Radial distance Nuclei center to Nuclei center (pixels)") +
    ggplot2::scale_color_discrete(
        "From -> To") +
    ggplot2::facet_grid(
        rows = dplyr::vars(Compound),
        cols = dplyr::vars(plate_id)) +
    ggplot2::ggtitle(
        label = "Assortive clustering over time")

ggplot2::ggsave(
    filename = "product/figures/pseudo_time_202006/pcf_lattice_20200714.pdf",
    plot = plot,
    height = 10,
    width = 10)

            
        
    

field <- cell_features %>%
    dplyr::filter(
        plate_id == "48h",
        Image_Metadata_WellID == "0311",
        Image_Metadata_FieldID == "0004") %>%
    dplyr::filter(is_infected)

point_process <- spatstat::ppp(
    x = field$Nuclei_Location_Center_X,
    y = field$Nuclei_Location_Center_Y,
    window = spatstat::owin(
        xrange = c(0, 2560),
        yrange = c(0, 2160)),
    marks = field %>%
        dplyr::mutate(
            is_infected = factor(is_infected, levels=c(TRUE, FALSE), labels=c("Infected", "Not Infected"))) %>%
        dplyr::select(is_infected))


pdf("product/figures/pseudo_time_202006/example_Kest.pdf", width=10, height=10)
par(mfrow=c(2,2))
point_process %>% spatstat::Kcross(from = "Infected", to = "Infected") %>% plot()
point_process %>% spatstat::Kcross(from = "Infected", to = "Not Infected") %>% plot()
point_process %>% spatstat::Kcross(from = "Not Infected", to = "Infected") %>% plot()
point_process %>% spatstat::Kcross(from = "Not Infected", to = "Not Infected") %>% plot()
dev.off()

pdf("product/figures/pseudo_time_202006/example_Kcross_pfc.pdf", width=10, height=10)
par(mfrow=c(2,2))
point_process %>%
    spatstat::Kcross(from = "Infected", to = "Infected") %>%
    spatstat::pcf(spar = 0.8, method = "b") %>%
    plot(main = "Pair Correlation Function: Infected --> Infected")
point_process %>% spatstat::Kcross(from = "Infected", to = "Not Infected") %>%
    spatstat::pcf(spar = 0.8, method = "b") %>% 
    plot(main = "Pair Correlation Function: Infected --> Not Infected")
point_process %>% spatstat::Kcross(from = "Not Infected", to = "Infected") %>%
    spatstat::pcf(spar = 0.8, method = "b") %>%
    plot(main = "Pair Correlation Function: Not Infected --> Infected")
point_process %>% spatstat::Kcross(from = "Not Infected", to = "Not Infected") %>%
    spatstat::pcf(spar = 0.8, method = "b") %>%
    plot(main = "Pair Correlation Function: Not Infected --> Not Infected")
dev.off()


pdf("product/figures/pseudo_time_202006/example_miplot.pdf")
point_process %>%
    spatstat::miplot()
dev.off()

pdf("product/figures/pseudo_time_202006/example_fry.pdf")
point_process %>%
    spatstat::fryplot(
        from = spatstat::marks(point_process) == "Infected",
        to = spatstat::marks(point_process) == "Infected",
        main = "Fry plot from infected to infected",
        pch = 16, cex=.3)
point_process %>%
    spatstat::fryplot(
        from = spatstat::marks(point_process) == "Infected",
        to = spatstat::marks(point_process) == "Not Infected",
        main = "Fry plot from infected to not infected",
        pch = 16, cex=.3)
point_process %>%
    spatstat::fryplot(
        from = spatstat::marks(point_process) == "Not Infected",
        to = spatstat::marks(point_process) == "Infected",
        main = "Fry plot from not infected to infected",
        pch = 16, cex=.3)
point_process %>%
    spatstat::fryplot(
        from = spatstat::marks(point_process) == "Not Infected",
        to = spatstat::marks(point_process) == "Not Infected",
        main = "Fry plot from not infected to not infected",
        pch = 16, cex=.3)
dev.off()




spdep::can.be.simmed # depricated?

?spdep::nb2listw(
    neighbors = 

spdep::moran.mc(
    x = field::cell_index,
    





distance_graph <- cell_features %>%
    dplyr::select(
        plate_id,
        Image_Metadata_WellID,
        Image_Metadata_FieldID,
        cell_index,
        Nuclei_Location_Center_X,
        Nuclei_Location_Center_Y) %>%
    plyr::ddply(
        c("plate_id", "Image_Metadata_WellID", "Image_Metadata_FieldID"),
        function(field) {
            if (nrow(field) < 50) {
                return(data.frame())
            }
            cat("plate: ", as.character(field$plate_id[1]),
                " Well: ", field$Image_Metadata_WellID[1],
                " Field: ", field$Image_Metadata_FieldID[1], "\n",
                sep = "")
            knn_neighbors <- field %>%
                dplyr::select(
                    Nuclei_Location_Center_X,
                    Nuclei_Location_Center_Y) %>%            
                BiocNeighbors::findKNN(k = 8)
            index <- as(knn_neighbors$index, "dgTMatrix")
            distance <- as(knn_neighbors$distance, "dgTMatrix")
            data.frame(
                from = field$cell_index[index@i + 1],
                to = field$cell_index[index@x],
                value = distance@x)
            })

graph <- distance_graph %>%
    dplyr::select(from, to) %>%
    igraph::graph_from_data_frame(
        directed = TRUE,
        vertices = cell_features %>%
            dplyr::select(
                cell_index,
                is_infected)

distance_graph %>%
    plyr::ddply(
        c("plate_id", "Image_Metadata_WellID", "Image_Metadata_FieldID"),
        function(field) {
            graph <- field %>%
                dplyr::select(from, to) %>%
                igraph::graph_from_data_frame(directed = TRUE)
            centrality_scores <- graph %>%
                igraph::eigen_centrality(
                    directed = TRUE,
                    
                

distance_graph <- distance_graph %>%
