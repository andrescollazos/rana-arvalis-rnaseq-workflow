rm(list = ls())
setwd(file.path(Sys.getenv("THESIS_DIR"), "code/04_de"))
load("02_interpretation.RData")

library(DESeq2)

lfc_mat_shrunk <- matrix(NA, nrow = nrow(dds_group), ncol = length(populations))
rownames(lfc_mat_shrunk) <- rownames(dds_group)
colnames(lfc_mat_shrunk) <- populations

for (pop in populations) {
    res <- lfcShrink(
        dds_group,
        contrast = c("group", paste0(pop, "_20"), paste0(pop, "_15")),
        type = "ashr"
    )
    lfc_mat_shrunk[, pop] <- res$log2FoldChange
    # padj_mat stays the same — shrinkage doesn't change significance calls with ashr
}

temp_summary_table_shrunk <- data.frame(
    population = populations,
    n_sig = sapply(populations, function(pop) {
        sum(!is.na(padj_mat[, pop]) & padj_mat[, pop] < 0.05)
    }),
    n_up = sapply(populations, function(pop) {
        sum(!is.na(padj_mat[, pop]) &
            padj_mat[, pop] < 0.05 &
            lfc_mat_shrunk[, pop] > 0)
    }),
    n_down = sapply(populations, function(pop) {
        sum(!is.na(padj_mat[, pop]) &
            padj_mat[, pop] < 0.05 &
            lfc_mat_shrunk[, pop] < 0)
    }),
    mean_abs_log2FC = sapply(populations, function(pop) {
        idx <- !is.na(padj_mat[, pop]) &
            padj_mat[, pop] < 0.05 &
            !is.na(lfc_mat_shrunk[, pop])
        if (!any(idx)) {
            return(NA_real_)
        }
        mean(abs(lfc_mat_shrunk[idx, pop]))
    }),
    median_abs_log2FC = sapply(populations, function(pop) {
        idx <- !is.na(padj_mat[, pop]) &
            padj_mat[, pop] < 0.05 &
            !is.na(lfc_mat_shrunk[, pop])
        if (!any(idx)) {
            return(NA_real_)
        }
        median(abs(lfc_mat_shrunk[idx, pop]))
    }),
    var_abs_log2FC = sapply(populations, function(pop) {
        idx <- !is.na(padj_mat[, pop]) &
            padj_mat[, pop] < 0.05 &
            !is.na(lfc_mat_shrunk[, pop])
        if (sum(idx) < 2) {
            return(NA_real_)
        }
        var(abs(lfc_mat_shrunk[idx, pop]))
    })
)

View(temp_summary_table_shrunk)
View(temp_summary_table)

library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)

# Build long-format data frame of significant genes only
sig_lfc_df <- do.call(rbind, lapply(populations, function(pop) {
    idx <- !is.na(padj_mat[, pop]) &
        padj_mat[, pop] < 0.05 &
        !is.na(lfc_mat_shrunk[, pop])
    if (!any(idx)) {
        return(NULL)
    }
    data.frame(
        population = pop,
        log2FC = lfc_mat_shrunk[idx, pop]
    )
}))

# Order populations by n_sig (or set your preferred order)
sig_lfc_df$population <- factor(
    sig_lfc_df$population,
    levels = c("VF", "C.Fin", "L", "NA", "Upp", "NL", "E", "Ka") # adjust as you prefer
)

# Violin plot
pdf("2.profile_lfc.pdf")
ggplot(sig_lfc_df, aes(x = population, y = log2FC, fill = population)) +
    geom_violin(scale = "width", trim = FALSE, alpha = 0.7, color = NA) +
    geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white", alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    labs(
        x = "Population",
        y = expression("Shrunken log"[2] * " fold change (20°C vs 15°C)"),
        title = "Effect-size distribution of temperature-responsive genes",
        subtitle = "Significant genes only (padj < 0.05)"
    ) +
    theme_classic(base_size = 12) +
    theme(
        legend.position = "none",
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 0)
    )
dev.off()

# Violin plot + points
pdf("2.profile_lfc_violin_points.pdf")
ggplot(sig_lfc_df, aes(x = population, y = log2FC, fill = population)) +
    geom_violin(scale = "width", trim = FALSE, alpha = 0.6, color = NA) +
    geom_jitter(width = 0.15, size = 0.4, alpha = 0.3) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    labs(
        x = "Population",
        y = expression("Shrunken log"[2] * " fold change (20°C vs 15°C)"),
        title = "Effect-size distribution of temperature-responsive genes",
        subtitle = "Significant genes only (padj < 0.05)"
    ) +
    theme_classic(base_size = 12) +
    theme(
        legend.position = "none",
        plot.title = element_text(face = "bold")
    )
dev.off()


save.image("03_lfc_shrunken.RData")


# -----------------------------
# 1. Differential plasticity across populations in response to temperature
# -----------------------------

lfc_interaction <- lfc_mat_shrunk[interaction_genes, ]
lfc_scaled <- t(scale(t(lfc_interaction)))
lfc_scaled <- lfc_scaled[complete.cases(lfc_scaled), ]

# Column annotations from meta (no redefinition)
pop_annot <- unique(meta[, c("population", "lat_group", "lineage")])
pop_annot <- pop_annot[match(colnames(lfc_scaled), pop_annot$population), ]
rownames(pop_annot) <- pop_annot$population
pop_annot$population <- NULL

pdf("1.differential_plasticity_shrunken.pdf")
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

cor_mat <- cor(
    lfc_scaled[interaction_genes, ],
    use = "pairwise.complete.obs",
    method = "pearson"
)
pdf("1.correlation_plasticity_shrunken.pdf")
pheatmap(
    cor_mat,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = pop_annot,
    main = "Correlation of differential plasticity (LRT genes)"
)
dev.off()
