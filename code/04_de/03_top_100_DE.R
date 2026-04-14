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
