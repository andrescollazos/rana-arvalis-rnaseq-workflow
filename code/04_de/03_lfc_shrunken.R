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

# Order populations by n_sig
sig_lfc_df$population <- factor(
    sig_lfc_df$population,
    levels = c("VF", "C.Fin", "L", "NA", "Upp", "NL", "E", "Ka")
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

# Column annotations from meta
pop_annot <- unique(meta[, c("population", "lat_group", "lineage")])
pop_annot <- pop_annot[match(colnames(lfc_scaled), pop_annot$population), ]
rownames(pop_annot) <- pop_annot$population
pop_annot$population <- NULL

colnames(pop_annot) <- c("Latitude", "Region")

pop_annot$Region <- factor(
    pop_annot$Region,
    levels = c("North", "South", "East"),
    labels = c("North Sweden", "South Sweden", "East")
)

pop_annot$Latitude <- factor(
    pop_annot$Latitude,
    levels = c("North", "South")
)

annotation_colors <- list(
    Region = region_colors,
    Latitude = latitude_colors
)

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
dev.off()

cor_mat <- cor(
    lfc_scaled[interaction_genes, ],
    use = "pairwise.complete.obs",
    method = "pearson"
)
png(
    "1.correlation_plasticity_shrunken_scaled.png",
    width = 2400,
    height = 2000,
    res = 300
)

pheatmap(
    cor_mat,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = pop_annot,
    annotation_colors = annotation_colors,
    main = "Correlation of differential plasticity profiles",
    display_numbers = TRUE,
    number_format = "%.2f"
)

dev.off()


# -----------------------------
# PARAMETERS
# -----------------------------
top_n <- 100

# -----------------------------
# HELPER FUNCTIONS
# -----------------------------

# row-wise z-score scaling
scale_rows <- function(mat) {
    t(apply(mat, 1, function(x) {
        if (all(is.na(x))) {
            return(rep(NA_real_, length(x)))
        }
        s <- sd(x, na.rm = TRUE)
        m <- mean(x, na.rm = TRUE)

        if (is.na(s) || s == 0) {
            return(rep(0, length(x)))
        } else {
            return((x - m) / s)
        }
    }))
}

# average pairwise absolute difference across populations
avg_pairwise_diff <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2) {
        return(NA_real_)
    }
    mean(dist(x, method = "manhattan"))
}

# -----------------------------
# Top 100 genes for all the interaction genes
# -----------------------------

# -----------------------------
# SCALAR 1: variance across populations
# -----------------------------
var_score <- apply(lfc_interaction, 1, function(x) {
    if (sum(!is.na(x)) < 2) {
        return(NA_real_)
    }
    var(x, na.rm = TRUE)
})

top_var_genes <- names(sort(var_score, decreasing = TRUE))[1:min(top_n, sum(!is.na(var_score)))]
lfc_top_var <- lfc_interaction[top_var_genes, , drop = FALSE]
lfc_top_var_scaled <- scale_rows(lfc_top_var)

pdf("3.top_100_DE_variance_shrunken.pdf")
p <- pheatmap(
    lfc_top_var_scaled,
    show_rownames = FALSE,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = pop_annot,
    main = paste0("Top ", nrow(lfc_top_var_scaled), " interaction genes by variance across populations")
)
print(p)
dev.off()

cor_mat_top_var <- cor(
    lfc_top_var_scaled,
    use = "pairwise.complete.obs",
    method = "pearson"
)

pdf("3.top_100_DE_variance_correlation_shrunken.pdf")
pheatmap(
    cor_mat_top_var,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = pop_annot,
    display_numbers = TRUE,
    number_format = "%.2f",
    main = paste0(
        "Correlation of differential plasticity\n(Top ",
        nrow(lfc_top_var_scaled),
        " interaction genes by variance, shrunken)"
    )
)
dev.off()

# -----------------------------
# SCALAR 2: mean absolute deviation from gene mean
# -----------------------------
mad_score <- apply(lfc_interaction, 1, function(x) {
    x2 <- x[!is.na(x)]
    if (length(x2) < 2) {
        return(NA_real_)
    }
    mean(abs(x2 - mean(x2)))
})

top_mad_genes <- names(sort(mad_score, decreasing = TRUE))[1:min(top_n, sum(!is.na(mad_score)))]
lfc_top_mad <- lfc_interaction[top_mad_genes, , drop = FALSE]
lfc_top_mad_scaled <- scale_rows(lfc_top_mad)

pdf("3.top_100_DE_mad_shrunken.pdf")
p <- pheatmap(
    lfc_top_mad_scaled,
    show_rownames = FALSE,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = pop_annot,
    main = paste0("Top ", nrow(lfc_top_mad_scaled), " interaction genes by mean absolute deviation, shrunken")
)
print(p)
dev.off()

cor_mat_top_mad <- cor(
    lfc_top_mad_scaled,
    use = "pairwise.complete.obs",
    method = "pearson"
)

pdf("3.top_100_DE_mad_correlation_shrunken.pdf")
pheatmap(
    cor_mat_top_mad,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = pop_annot,
    display_numbers = TRUE,
    number_format = "%.2f",
    main = paste0(
        "Correlation of differential plasticity\n(Top ",
        nrow(lfc_top_mad_scaled),
        " interaction genes by mean absolute deviation, shrunken)"
    )
)
dev.off()

# -----------------------------
# SCALAR 3: average pairwise absolute difference
# -----------------------------
pairdiff_score <- apply(lfc_interaction, 1, avg_pairwise_diff)

top_pairdiff_genes <- names(sort(pairdiff_score, decreasing = TRUE))[1:min(top_n, sum(!is.na(pairdiff_score)))]
lfc_top_pairdiff <- lfc_interaction[top_pairdiff_genes, , drop = FALSE]
lfc_top_pairdiff_scaled <- scale_rows(lfc_top_pairdiff)

pdf("3.top_100_DE_pairdiff.pdf")
p <- pheatmap(
    lfc_top_pairdiff_scaled,
    show_rownames = FALSE,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = pop_annot,
    main = paste0("Top ", nrow(lfc_top_pairdiff_scaled), " interaction genes by average pairwise difference")
)
print(p)
dev.off()

cor_mat_top_pairdiff <- cor(
    lfc_top_pairdiff_scaled,
    use = "pairwise.complete.obs",
    method = "pearson"
)

pdf("3.top_100_DE_pairdiff_correlation_shrunken.pdf")
pheatmap(
    cor_mat_top_pairdiff,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = pop_annot,
    display_numbers = TRUE,
    number_format = "%.2f",
    main = paste0(
        "Correlation of differential plasticity\n(Top ",
        nrow(lfc_top_pairdiff_scaled),
        " interaction genes by average pairwise difference, shrunken)"
    )
)
dev.off()
