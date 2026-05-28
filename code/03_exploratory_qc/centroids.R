rm(list = ls())
setwd("~/Documents/Studies/Thesis/code/03_exploratory_qc")
load("pca_baseline.RData")

library(dplyr)
library(ggplot2)
library(ggrepel)

make_centroid_pca_plots <- function(pca_result, plot_title = "PCA centroids") {
    pca_df <- pca_result$pca_df
    pc_pct <- pca_result$pc_pct

    centroid_df <- pca_df %>%
        group_by(population, temperature, region) %>%
        summarise(
            across(starts_with("PC"), mean, na.rm = TRUE),
            n_samples = n(),
            .groups = "drop"
        ) %>%
        mutate(
            region_temp_label = paste0(region, " ", temperature, "\u00B0C"),
            centroid_label = paste0(population, " ", temperature, "\u00B0C")
        )

    color_map <- c(
        "South 20\u00B0C" = "red4",
        "South 15\u00B0C" = "violetred1",
        "North 20\u00B0C" = "blue",
        "North 15\u00B0C" = "deepskyblue3",
        "East 20\u00B0C"  = "#088000",
        "East 15\u00B0C"  = "#5DBD56"
    )

    all_pops <- c("C.Fin", "E", "Ka", "L", "NA", "NL", "Upp", "VF")
    pch_map <- setNames(seq_along(all_pops), all_pops)

    plots <- vector("list", 5)

    for (i in 1:5) {
        xpc <- paste0("PC", i)
        ypc <- paste0("PC", i + 1)

        plots[[i]] <- ggplot(
            centroid_df,
            aes(x = .data[[xpc]], y = .data[[ypc]])
        ) +
            geom_vline(
                xintercept = 0,
                linetype = "dashed",
                color = "grey70",
                linewidth = 0.6
            ) +
            geom_hline(
                yintercept = 0,
                linetype = "dashed",
                color = "grey70",
                linewidth = 0.6
            ) +
            geom_point(
                aes(
                    color = region_temp_label,
                    shape = population
                ),
                size = 4,
                stroke = 1.2
            ) +
            geom_text_repel(
                aes(label = centroid_label),
                size = 3.5,
                box.padding = 0.5,
                point.padding = 0.3,
                min.segment.length = 0,
                max.overlaps = Inf,
                show.legend = FALSE
            ) +
            scale_color_manual(
                values = color_map,
                name = "Region x temperature"
            ) +
            scale_shape_manual(
                values = pch_map,
                name = "Population",
                na.translate = FALSE
            ) +
            labs(
                title = plot_title,
                x = paste0(xpc, " (", round(pc_pct[i], 2), "%)"),
                y = paste0(ypc, " (", round(pc_pct[i + 1], 2), "%)")
            ) +
            theme_classic(base_size = 14) +
            theme(
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
                axis.line = element_line(colour = "black"),
                legend.position = "right"
            )
    }

    list(
        centroid_df = centroid_df,
        plots = plots
    )
}

pca_baseline_centroids <- make_centroid_pca_plots(
    pca_result = pca_baseline,
    plot_title = "PCA centroids (All filtered genes)"
)

ggsave(
    filename = "08_PCA_baseline_centroids.pdf",
    plot = gridExtra::marrangeGrob(
        grobs = pca_baseline_centroids$plots,
        nrow = 2,
        ncol = 1
    ),
    width = 10,
    height = 12
)
