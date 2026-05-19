rm(list = ls())
setwd(file.path(Sys.getenv("THESIS_DIR"), "code/04_de"))

load("resultsTemperatureEffects.RData")

library(pheatmap)
library(ggplot2)
library(gridExtra)
library(ggrepel)

vsd <- vst(dds_group, blind = FALSE)
vst_mat <- assay(vsd)

for (pop in populations) {
    genes_pop <- temp_genes_by_pop[[pop]]

    samples_pop <- rownames(colData(vsd))[
        colData(vsd)$population == pop
    ]

    mat_pop <- vst_mat[genes_pop, samples_pop, drop = FALSE]

    # Remove genes with no variation across selected samples
    mat_pop <- mat_pop[apply(mat_pop, 1, sd, na.rm = TRUE) > 0, , drop = FALSE]

    if (nrow(mat_pop) < 2) {
        message(pop, ": fewer than 2 variable genes after filtering; skipping heatmap.")
        next
    }

    # Row-scale genes to emphasize relative expression patterns across samples
    mat_pop_scaled <- t(scale(t(mat_pop)))

    annot_col <- data.frame(
        temperature = colData(vsd)[samples_pop, "temperature"]
    )
    rownames(annot_col) <- samples_pop

    annotation_colors <- list(
        temperature = temperature_colors
    )

    pdf(paste0("results/population-specific/", pop, "_heatmap_temp_sig.pdf"))
    pheatmap(
        mat_pop_scaled,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        annotation_col = annot_col,
        annotation_colors = annotation_colors,
        show_rownames = FALSE,
        main = paste0(pop, ": significant temperature-responsive genes")
    )
    dev.off()

    expr_mat <- mat_pop

    pca <- prcomp(t(expr_mat), center = TRUE, scale. = FALSE)
    pc_pct <- (pca$sdev^2) / sum(pca$sdev^2) * 100

    scores <- as.data.frame(pca$x)
    scores$sample <- rownames(scores)

    meta_pop <- meta[match(scores$sample, meta$sample), , drop = FALSE]

    stopifnot(!any(is.na(meta_pop$sample)))
    stopifnot(all(meta_pop$sample == scores$sample))

    pca_df <- cbind(scores, meta_pop[, setdiff(colnames(meta_pop), "sample"), drop = FALSE])

    pca_df$temperature <- factor(pca_df$temperature)

    all_pops <- c("C.Fin", "E", "Ka", "L", "NA", "NL", "Upp", "VF")
    pch_map <- setNames(seq_along(all_pops), all_pops)

    pca_df$population <- factor(
        pca_df$population,
        levels = all_pops
    )

    plots <- list()

    max_pc_pair <- min(5, ncol(pca$x) - 1)

    for (i in seq_len(max_pc_pair)) {
        xpc <- paste0("PC", i)
        ypc <- paste0("PC", i + 1)

        p <- ggplot(
            pca_df,
            aes(x = .data[[xpc]], y = .data[[ypc]], color = temperature)
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
                aes(color = temperature, shape = population),
                size = 4,
                stroke = 1.2
            ) +
            geom_text_repel(
                aes(label = sample),
                size = 3.5,
                box.padding = 0.5,
                point.padding = 0.3,
                min.segment.length = 0,
                max.overlaps = Inf,
                show.legend = FALSE
            ) +
            scale_color_manual(
                values = temperature_colors_bold,
                drop = FALSE
            ) +
            scale_shape_manual(
                values = pch_map,
                guide = "none",
                drop = FALSE
            ) +
            labs(
                title = paste0(pop, ": PCA of significant temperature-responsive genes"),
                x = paste0(xpc, " (", round(pc_pct[i], 2), "%)"),
                y = paste0(ypc, " (", round(pc_pct[i + 1], 2), "%)"),
                color = "Temperature"
            ) +
            theme_classic(base_size = 14) +
            theme(
                panel.border = element_rect(
                    colour = "black",
                    fill = NA,
                    linewidth = 0.8
                ),
                axis.line = element_line(colour = "black"),
                legend.position = "right"
            )

        plots[[i]] <- p
    }

    ggsave(
        filename = paste0(
            "results/population-specific/",
            pop,
            "_pca_temp_sig.pdf"
        ),
        plot = gridExtra::marrangeGrob(
            grobs = plots,
            nrow = 2,
            ncol = 1
        ),
        width = 10,
        height = 12
    )
}
