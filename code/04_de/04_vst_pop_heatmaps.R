setwd(file.path(Sys.getenv("THESIS_DIR"), "code/04_de"))

load("resultsTemperatureEffects.RData")

library(pheatmap)

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

    pdf(paste0(pop, "_sign_temp_responsive_genes.pdf"))
    pheatmap(
        mat_pop_scaled,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        annotation_col = annot_col,
        show_rownames = FALSE,
        main = paste0(pop, ": significant temperature-responsive genes")
    )
    dev.off()
}
