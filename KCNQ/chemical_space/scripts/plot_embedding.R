#' Plot substance embedding
#' Group label by substance source
plot_embedding <- function(
    substance_data,
    label,
    dataset_tag,
    background_data = NULL,
    plot_tag,
    save_plot = TRUE) {
    output_path <- paste0("product/figures/", dataset_tag)
    if (!dir.exists(output_path)) {
        cat("creating '", output_path, "'\n", sep = "")
        dir.create(output_path)
    }
    label_data <- substance_data %>%
        dplyr::group_by({{label}}) %>%
        dplyr::summarize(
            UMAP_1 = mean(UMAP_1),
            UMAP_2 = mean(UMAP_2),
            .groups = "drop")
    mapping <- ggplot2::aes(
                x = UMAP_1,
                y = UMAP_2,
                color = {{label}})
    plot <- ggplot2::ggplot()

    if (!is.null(background_data)) {
        cat("Adding background with ", nrow(background_data), " ...\n", sep = "")
        plot <- plot +
            ggplot2::geom_point(
                data = background_data,
                mapping = ggplot2::aes(
                    x = UMAP_1,
                    y = UMAP_2),
                color = "darkgrey",
                size = .8,
                alpha = .3,
                shape = 16)
    }
    
    plot <- plot +
        ggplot2::geom_point(
            data = substance_data,
            mapping = ggplot2::aes(
                x = UMAP_1,
                y = UMAP_2,
                color = {{label}}),
            size = 1.2) +
        ggrepel::geom_text_repel(
            data = label_data,
            mapping = ggplot2::aes(
                x = UMAP_1,
                y = UMAP_2,
                label = {{label}}),
            force = 10) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none") +
        ggplot2::scale_x_continuous("UMAP 1") +
        ggplot2::scale_y_continuous("UMAP 2") +
        ggplot2::coord_fixed() +
        ggplot2::ggtitle(
            label = "Emedding based on ECP4 fingerprints",
            subtitle = "Literature reported substances")
    if(save_plot){
        fname <- paste0(
            "product/figures/",
            dataset_tag, "/",
            "embedding_", plot_tag, ".pdf")
        cat("saving plot '", fname, "' ...\n", sep = "")
        ggplot2::ggsave(
            filename = fname,
            width = 12,
            height = 12,
            useDingbats = FALSE)
    
        fname <- paste0(
            "product/figures/",
            dataset_tag, "/",
            "embedding_", plot_tag, ".png")
        cat("saving plot '", fname, "' ...\n", sep = "")
        ggplot2::ggsave(
            filename = fname,
            width = 12,
            height = 12)
    }
    plot
}
