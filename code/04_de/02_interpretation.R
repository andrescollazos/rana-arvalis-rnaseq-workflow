load("resultsInteraction.RData")
load("resultsTemperatureEffects.RData")
library(pheatmap)
library(ggplot2)
library(reshape2)

# Subset + scale
lfc_interaction <- lfc_mat[interaction_genes, ]
lfc_scaled <- t(scale(t(lfc_interaction)))
lfc_scaled <- lfc_scaled[complete.cases(lfc_scaled), ]

# Column annotations from meta (no redefinition)
pop_annot <- unique(meta[, c("population", "lat_group", "lineage")])
pop_annot <- pop_annot[match(colnames(lfc_scaled), pop_annot$population), ]
rownames(pop_annot) <- pop_annot$population
pop_annot$population <- NULL

# -----------------------------
# 1. With column clustering
# -----------------------------
pdf("1.differential_plasticity.pdf")
p <- pheatmap(
    lfc_scaled,
    show_rownames = FALSE,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = pop_annot,
    main = "Differential plasticity across populations in response to temperature"
)
print(p)
pop_order <- c("NA", "NL", "VF", "C.Fin", "E", "L", "Upp", "Ka")

lfc_scaled_ordered <- lfc_scaled[, pop_order]
pop_annot_ordered <- pop_annot[pop_order, ]

p2 <- pheatmap(
    lfc_scaled_ordered,
    show_rownames = FALSE,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_col = pop_annot_ordered,
    main = "Differential plasticity across populations in response to temperature (fixed columns)"
)
print(p2)
dev.off()
