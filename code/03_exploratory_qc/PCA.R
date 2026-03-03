make_pca_plots <- function(
  vsd_mat,
  meta,
  gene_mode = c("baseline", "top_variable_100", "top_variable_500"),
  plot_title = ""
) {
    gene_mode <- match.arg(gene_mode)

    # ----------------------------
    # 1. Gene selection
    # ----------------------------
    if (gene_mode == "top_variable_100") {
        gene_var <- apply(vsd_mat, 1, var)
        top_n <- min(100, length(gene_var))
        top_idx <- order(gene_var, decreasing = TRUE)[1:top_n]
        expr_mat <- vsd_mat[top_idx, , drop = FALSE]
    } else if (gene_mode == "top_variable_500") {
        gene_var <- apply(vsd_mat, 1, var)
        top_n <- min(500, length(gene_var))
        top_idx <- order(gene_var, decreasing = TRUE)[1:top_n]
        expr_mat <- vsd_mat[top_idx, , drop = FALSE]
    } else {
        expr_mat <- vsd_mat
    }

    # ----------------------------
    # 2. PCA
    # ----------------------------
    pca <- prcomp(t(expr_mat), center = TRUE, scale. = FALSE)
    pc_pct <- (pca$sdev^2) / sum(pca$sdev^2) * 100

    scores <- as.data.frame(pca$x)
    scores$sample <- rownames(scores)
    stopifnot(!anyDuplicated(meta$sample))
    stopifnot(all(scores$sample %in% meta$sample))
    stopifnot(all(meta$sample %in% scores$sample))

    pca_df <- merge(
        scores,
        meta,
        by.x = "sample",
        by.y = "sample",
        sort = FALSE
    )

    stopifnot(nrow(pca_df) == nrow(scores))
    stopifnot(all(pca_df$sample == scores$sample))

    # ----------------------------
    # 3. Region x temperature colors
    # ----------------------------
    pca_df$region_temp_label <- paste0(
        as.character(pca_df$region),
        " ",
        as.character(pca_df$temperature),
        "\u00B0C"
    )

    color_map <- c(
        "South 20\u00B0C" = "red4",
        "South 15\u00B0C" = "violetred1",
        "North 20\u00B0C" = "blue",
        "North 15\u00B0C" = "deepskyblue3",
        "East 20\u00B0C"  = "#088000",
        "East 15\u00B0C"  = "#5DBD56"
    )

    pca_df$color <- unname(color_map[pca_df$region_temp_label])

    # ----------------------------
    # 4. Shape mapping & Labeling
    # ----------------------------
    populations <- unique(pca_df$population)
    pch_map <- setNames(seq_along(populations), populations)
    pca_df$shape <- pch_map[pca_df$population]

    # Specific samples to label
    labeled_samples <- c("P32262_306", "P32262_274", "P32262_272", "P32262_248")
    pca_df$label <- ifelse(pca_df$sample %in% labeled_samples, as.character(pca_df$sample), NA)

    # ----------------------------
    # 5. Build PC1vsPC2 ... PC5vsPC6
    # ----------------------------
    plots <- vector("list", 5)

    for (i in 1:5) {
        xpc <- paste0("PC", i)
        ypc <- paste0("PC", i + 1)

        p <- ggplot(
            pca_df,
            aes(x = .data[[xpc]], y = .data[[ypc]])
        ) +
            geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.6) +
            geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.6) +
            geom_point(
                aes(color = color, shape = population),
                size = 3
            ) +
            geom_text_repel(
                aes(label = label),
                size = 3.5,
                box.padding = 0.5,
                point.padding = 0.3,
                min.segment.length = 0,
                max.overlaps = Inf,
                show.legend = FALSE
            ) +
            scale_color_identity(
                guide = "legend",
                name = "Region x temperature",
                breaks = unname(color_map),
                labels = names(color_map)
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

        plots[[i]] <- p
    }

    # ----------------------------
    # 6. Variance Explained Summary
    # ----------------------------
    pca_variance <- data.frame(
        PC = paste0("PC", seq_along(pc_pct)),
        variance_explained_pct = pc_pct,
        cumulative_variance_pct = cumsum(pc_pct)
    )

    return(list(
        pca_df = pca_df,
        pc_pct = pc_pct,
        pca_variance = pca_variance,
        plots = plots
    ))
}
