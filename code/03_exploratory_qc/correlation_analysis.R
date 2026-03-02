run_correlation_analysis <- function(
  vsd_mat,
  meta,
  gene_mode = c("baseline", "top_variable"),
  cluster_method = "average",
  scale_min = 0.79,
  scale_max = 1,
  palette_mode = c("red_white", "purple_yellow"),
  plot_title = "Correlation Heatmap"
) {
    gene_mode <- match.arg(gene_mode)
    color_values <- if (palette_mode == "red_white") {
        colorRampPalette(c("#313695", "#ffffbf", "#a50026"))(100)
    } else {
        colorRampPalette(c("purple", "black", "yellow"))(100)
    }

    # ----------------------------
    # 1. Gene selection
    # ----------------------------
    if (gene_mode == "top_variable") {
        gene_var <- apply(vsd_mat, 1, var)
        top_n <- min(500, length(gene_var))
        top_idx <- order(gene_var, decreasing = TRUE)[1:top_n]
        expr_mat <- vsd_mat[top_idx, , drop = FALSE]
    } else {
        expr_mat <- vsd_mat
    }

    # ----------------------------
    # 2. Pearson correlation
    # ----------------------------
    cor_mat <- cor(expr_mat, method = "pearson")

    # ----------------------------
    # 3. Hierarchical clustering
    # ----------------------------
    dist_mat <- as.dist(1 - cor_mat)
    hc <- hclust(dist_mat, method = cluster_method)
    ord <- hc$order
    cor_ord <- cor_mat[ord, ord]

    # ----------------------------
    # 4. Pairwise sample correlations
    # ----------------------------
    upper_idx <- which(upper.tri(cor_mat), arr.ind = TRUE)
    pairwise_tbl <- data.frame(
        sample1 = rownames(cor_mat)[upper_idx[, 1]],
        sample2 = colnames(cor_mat)[upper_idx[, 2]],
        correlation = cor_mat[upper.tri(cor_mat)],
        stringsAsFactors = FALSE
    )

    pop_map <- setNames(as.character(meta$population), rownames(meta))
    pop1 <- pop_map[pairwise_tbl$sample1]
    pop2 <- pop_map[pairwise_tbl$sample2]
    pairwise_tbl$group_label <- ifelse(pop1 == pop2, "within_population", "between_population")

    within_vals <- pairwise_tbl$correlation[pairwise_tbl$group_label == "within_population"]
    between_vals <- pairwise_tbl$correlation[pairwise_tbl$group_label == "between_population"]

    # ----------------------------
    # 5. Summary statistics
    # ----------------------------
    summary_tbl <- data.frame(
        group_label = c("within_population", "between_population"),
        mean_correlation = c(mean(within_vals, na.rm = TRUE), mean(between_vals, na.rm = TRUE)),
        stringsAsFactors = FALSE
    )

    heatmap_plot <- function(file = NULL, width = 3000, height = 3000, res = 300) {
        if (!requireNamespace("pheatmap", quietly = TRUE)) stop("Install pheatmap.")

        annotation_df <- meta[, c("population", "region", "temperature"), drop = FALSE]
        annotation_df <- annotation_df[colnames(cor_ord), , drop = FALSE]

        # Ensure factors (important for correct legend display)
        annotation_df$population <- factor(annotation_df$population)
        annotation_df$region <- factor(annotation_df$region)
        annotation_df$temperature <- factor(annotation_df$temperature)

        annotation_colors <- list(
            population = c(
                "Ka" = "#FFDE52",
                "Upp" = "#FF7070",
                "NA" = "#6f9fe2ff",
                "NL" = "#59FFC0",
                "VF" = "#6a3d9a",
                "C.Fin" = "#FF33E4",
                "E" = "#B7FF5E",
                "L" = "#00CDFF"
            ),
            region = c(
                "North" = "#1f78b4",
                "South" = "#e31a1c",
                "East"  = "#33a02c"
            ),
            temperature = c(
                "15" = "#a6cee3",
                "20" = "#fb9a99"
            )
        )

        if (!is.null(file)) jpeg(file, width = width, height = height, res = res)

        summary_label <- paste0(
            "Mean sample Pearson correlation | Within: ", round(summary_tbl$mean_correlation[1], 4),
            " | Between: ", round(summary_tbl$mean_correlation[2], 4)
        )

        pheatmap_obj <- pheatmap::pheatmap(
            cor_ord,
            cluster_rows = hc,
            cluster_cols = hc,
            annotation_col = annotation_df,
            annotation_row = annotation_df,
            annotation_colors = annotation_colors,
            color = color_values,
            breaks = seq(scale_min, scale_max, length.out = 101),
            border_color = NA,
            legend = TRUE,
            show_rownames = TRUE,
            show_colnames = TRUE,
            fontsize_row = 9,
            fontsize_col = 6,
            angle_col = 90,
            treeheight_row = 80,
            treeheight_col = 80,
            main = plot_title,
            silent = TRUE
        )

        grid::grid.newpage()
        grid::grid.draw(pheatmap_obj$gtable)

        # Add summary text BELOW the heatmap
        grid::grid.text(
            summary_label,
            x = 0.5,
            y = 0.01,
            just = c("center", "bottom"),
            gp = grid::gpar(fontsize = 12)
        )

        if (!is.null(file)) dev.off()
    }

    list(
        plot_heatmap = heatmap_plot
    )
}
