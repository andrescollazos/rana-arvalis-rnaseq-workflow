setwd(file.path(Sys.getenv("THESIS_DIR"), "code/04_de"))

# Load data
load("02_interpretation.RData")

library(pheatmap)

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

pdf("3.top_100_DE_variance.pdf")
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

pdf("3.top_100_DE_mad.pdf")
p <- pheatmap(
    lfc_top_mad_scaled,
    show_rownames = FALSE,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = pop_annot,
    main = paste0("Top ", nrow(lfc_top_mad_scaled), " interaction genes by mean absolute deviation")
)
print(p)
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

# -----------------------------
# Concordant genes
# -----------------------------
concordant_ids <- names(class_vec)[class_vec == "concordant"]
concordant_ids <- intersect(concordant_ids, rownames(lfc_mat))

lfc_concordant <- lfc_mat[concordant_ids, , drop = FALSE]

# remove genes with too few values
keep_conc <- rowSums(!is.na(lfc_concordant)) >= 2
lfc_concordant <- lfc_concordant[keep_conc, , drop = FALSE]

lfc_concordant_scaled <- scale_rows(lfc_concordant)

pdf("3.top_100_DE_concordant.pdf")
p <- pheatmap(
    lfc_concordant_scaled,
    show_rownames = FALSE,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = pop_annot,
    main = paste0("Concordant interaction genes (n = ", nrow(lfc_concordant_scaled), ")")
)
print(p)
dev.off()

# -----------------------------
# Top 100 genes for the concordant genes
# -----------------------------

# Variance across populations
var_score_conc <- apply(lfc_concordant, 1, function(x) {
    if (sum(!is.na(x)) < 2) {
        return(NA_real_)
    }
    var(x, na.rm = TRUE)
})

top_var_conc_genes <- names(sort(var_score_conc, decreasing = TRUE))[1:min(top_n, sum(!is.na(var_score_conc)))]
lfc_top_var_conc <- lfc_concordant[top_var_conc_genes, , drop = FALSE]
lfc_top_var_conc_scaled <- scale_rows(lfc_top_var_conc)

pdf("3.top_100_DE_variance_conc.pdf")
p <- pheatmap(
    lfc_top_var_conc_scaled,
    show_rownames = FALSE,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = pop_annot,
    main = paste0("Top ", nrow(lfc_top_var_conc_scaled), " concordant genes by variance across populations")
)
print(p)
dev.off()

# Mean absolute deviation from gene mean
mad_score_conc <- apply(lfc_concordant, 1, function(x) {
    x2 <- x[!is.na(x)]
    if (length(x2) < 2) {
        return(NA_real_)
    }
    mean(abs(x2 - mean(x2)))
})

top_mad_conc_genes <- names(sort(mad_score_conc, decreasing = TRUE))[1:min(top_n, sum(!is.na(mad_score_conc)))]
lfc_top_mad_conc <- lfc_concordant[top_mad_conc_genes, , drop = FALSE]
lfc_top_mad_conc_scaled <- scale_rows(lfc_top_mad_conc)

pdf("3.top_100_DE_mad_conc.pdf")
p <- pheatmap(
    lfc_top_mad_conc_scaled,
    show_rownames = FALSE,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = pop_annot,
    main = paste0("Top ", nrow(lfc_top_mad_conc_scaled), " concordant genes by mean absolute deviation")
)
print(p)
dev.off()

# Average pairwise absolute difference across populations
pairdiff_score_conc <- apply(lfc_concordant, 1, avg_pairwise_diff)

top_pairdiff_conc_genes <- names(sort(pairdiff_score_conc, decreasing = TRUE))[1:min(top_n, sum(!is.na(pairdiff_score_conc)))]
lfc_top_pairdiff_conc <- lfc_concordant[top_pairdiff_conc_genes, , drop = FALSE]
lfc_top_pairdiff_conc_scaled <- scale_rows(lfc_top_pairdiff_conc)

pdf("3.top_100_DE_pairdiff_conc.pdf")
p <- pheatmap(
    lfc_top_pairdiff_conc_scaled,
    show_rownames = FALSE,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = pop_annot,
    main = paste0("Top ", nrow(lfc_top_pairdiff_conc_scaled), " concordant genes by average pairwise difference")
)
print(p)
dev.off()

save.image("03_top_100_DE.RData")
