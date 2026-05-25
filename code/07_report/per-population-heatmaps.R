rm(list = ls())
setwd(file.path(Sys.getenv("THESIS_DIR"), "code/04_de"))
load("resultsTemperatureEffects.RData")

library(DESeq2)
library(pheatmap)
library(gridExtra)
library(grid)
library(gtable)

vsd <- vst(dds_group, blind = FALSE)
vst_mat <- assay(vsd)

panel_order <- c(
    "VF", "L",
    "C.Fin", "NA",
    "NL", "E",
    "Ka", "Upp"
)

panel_labels <- LETTERS[seq_along(panel_order)]

for (i in seq_along(panel_order)) {
    pop <- panel_order[i]

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

    show_legend <- i == length(panel_order)

    png(
        filename = paste0("results/population-specific/", pop, "_heatmap_temp_sig_scale.png"),
        width = 2000,
        height = 2000,
        res = 300
    )

    pheatmap(
        mat_pop_scaled,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        treeheight_col = 15,
        annotation_col = annot_col,
        annotation_colors = annotation_colors,
        show_rownames = FALSE,
        show_colnames = FALSE,
        legend = show_legend,
        annotation_legend = show_legend
    )

    dev.off()
}
