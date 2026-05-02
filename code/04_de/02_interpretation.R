setwd(file.path(Sys.getenv("THESIS_DIR"), "code/04_de"))
load("resultsInteraction.RData")
load("resultsTemperatureEffects.RData")
library(pheatmap)
library(ggplot2)
library(reshape2)

# -----------------------------
# 1. Differential plasticity across populations in response to temperature
# -----------------------------

lfc_interaction <- lfc_mat[interaction_genes, ]
lfc_scaled <- t(scale(t(lfc_interaction)))
lfc_scaled <- lfc_scaled[complete.cases(lfc_scaled), ]

# Column annotations from meta (no redefinition)
pop_annot <- unique(meta[, c("population", "lat_group", "lineage")])
pop_annot <- pop_annot[match(colnames(lfc_scaled), pop_annot$population), ]
rownames(pop_annot) <- pop_annot$population
pop_annot$population <- NULL

# -----------------------------
# 1. Differential plasticity across populations in response to temperature
# -----------------------------
pdf("1.differential_plasticity3.pdf")
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
pdf("1.correlation_plasticity.pdf")
pheatmap(
    cor_mat,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = pop_annot,
    main = "Correlation of differential plasticity (LRT genes)"
)
dev.off()

# -----------------------------
# 2. Distribution of plasticity across populations (pattern prevalence)
# -----------------------------
library(ComplexUpset)
library(dplyr)

# Binary matrix: significant or not per population
sig_binary <- as.data.frame(
    sapply(populations, function(pop) {
        !is.na(padj_mat[, pop]) & padj_mat[, pop] < 0.05
    }),
    stringsAsFactors = FALSE
)

sig_binary$gene_id <- rownames(padj_mat)

# Keep only genes significant in at least one population
sig_binary <- sig_binary[rowSums(sig_binary[, populations, drop = FALSE]) > 0, ]

# Ensure logical columns
sig_binary[, populations] <- lapply(sig_binary[, populations, drop = FALSE], as.logical)

pdf("2.distribution_plasticity2.pdf")
# UpSet plot
upset(
    sig_binary,
    intersect = populations,
    name = "Population",
    base_annotations = list(
        "Intersection size" = intersection_size(text = list(size = 2))
    ),
    min_size = 20
)
dev.off()

# Shared plasticity across all populations
shared_plasticity <- sig_binary$gene_id[
    rowSums(sig_binary[, populations, drop = FALSE]) == length(populations)
]

# Population-specific plasticity
population_specific <- lapply(populations, function(pop) {
    others <- setdiff(populations, pop)

    sig_binary$gene_id[
        sig_binary[[pop]] &
            rowSums(sig_binary[, others, drop = FALSE]) == 0
    ]
})

names(population_specific) <- populations

lat_north <- unique(meta$population[meta$lat_group == "North"])
lat_south <- unique(meta$population[meta$lat_group == "South"])

lin_north <- unique(meta$population[meta$lineage == "North"])
lin_east <- unique(meta$population[meta$lineage == "East"])
lin_south <- unique(meta$population[meta$lineage == "South"])

# Latitude-restricted (at least one in group, none outside)
latitude_restricted <- list(
    North = sig_binary$gene_id[
        rowSums(sig_binary[, lat_north, drop = FALSE]) == length(lat_north) &
            rowSums(sig_binary[, lat_south, drop = FALSE]) == 0
    ],
    South = sig_binary$gene_id[
        rowSums(sig_binary[, lat_south, drop = FALSE]) == length(lat_south) &
            rowSums(sig_binary[, lat_north, drop = FALSE]) == 0
    ]
)

# Lineage-restricted (all in group, none outside)

lineage_restricted <- list(
    North = sig_binary$gene_id[
        rowSums(sig_binary[, lin_north, drop = FALSE]) == length(lin_north) &
            rowSums(sig_binary[, setdiff(populations, lin_north), drop = FALSE]) == 0
    ],
    East = sig_binary$gene_id[
        rowSums(sig_binary[, lin_east, drop = FALSE]) == length(lin_east) &
            rowSums(sig_binary[, setdiff(populations, lin_east), drop = FALSE]) == 0
    ],
    South = sig_binary$gene_id[
        rowSums(sig_binary[, lin_south, drop = FALSE]) == length(lin_south) &
            rowSums(sig_binary[, setdiff(populations, lin_south), drop = FALSE]) == 0
    ]
)

upset_gene_sets <- list(
    shared_plasticity = shared_plasticity,
    population_specific = population_specific,
    latitude_restricted = latitude_restricted,
    lineage_restricted = lineage_restricted
)

list_lengths <- lapply(upset_gene_sets, function(x) {
    if (is.list(x) && !is.null(names(x))) {
        sapply(x, length)
    } else {
        length(x)
    }
})

list_lengths

# -----------------------------
# 3. Per-population summaries
# -----------------------------

temp_summary_table <- data.frame(
    population = populations,
    n_sig = sapply(populations, function(pop) {
        sum(!is.na(padj_mat[, pop]) & padj_mat[, pop] < 0.05)
    }),
    n_up = sapply(populations, function(pop) {
        sum(!is.na(padj_mat[, pop]) &
            padj_mat[, pop] < 0.05 &
            lfc_mat[, pop] > 0)
    }),
    n_down = sapply(populations, function(pop) {
        sum(!is.na(padj_mat[, pop]) &
            padj_mat[, pop] < 0.05 &
            lfc_mat[, pop] < 0)
    }),
    mean_abs_log2FC = sapply(populations, function(pop) {
        idx <- !is.na(padj_mat[, pop]) &
            padj_mat[, pop] < 0.05 &
            !is.na(lfc_mat[, pop])

        if (!any(idx)) {
            return(NA_real_)
        }
        mean(abs(lfc_mat[idx, pop]))
    }),
    median_abs_log2FC = sapply(populations, function(pop) {
        idx <- !is.na(padj_mat[, pop]) &
            padj_mat[, pop] < 0.05 &
            !is.na(lfc_mat[, pop])

        if (!any(idx)) {
            return(NA_real_)
        }
        median(abs(lfc_mat[idx, pop]))
    }),
    var_abs_log2FC = sapply(populations, function(pop) {
        idx <- !is.na(padj_mat[, pop]) &
            padj_mat[, pop] < 0.05 &
            !is.na(lfc_mat[, pop])

        if (sum(idx) < 2) {
            return(NA_real_)
        }
        var(abs(lfc_mat[idx, pop]))
    })
)
View(temp_summary_table)

# -----------------------------
# 4. Concordant vs discordant plasticity
# -----------------------------

# restrict to interaction genes
lfc_int <- lfc_mat[interaction_genes, ]
padj_int <- padj_mat[interaction_genes, ]

# direction matrix
dir_mat <- matrix(NA, nrow = nrow(lfc_int), ncol = ncol(lfc_int))
rownames(dir_mat) <- rownames(lfc_int)
colnames(dir_mat) <- colnames(lfc_int)

for (i in seq_len(nrow(lfc_int))) {
    for (j in seq_len(ncol(lfc_int))) {
        if (!is.na(padj_int[i, j]) && padj_int[i, j] < 0.05) {
            dir_mat[i, j] <- sign(lfc_int[i, j])
        }
    }
}

# classify genes
class_vec <- rep(NA, nrow(dir_mat))
names(class_vec) <- rownames(dir_mat)

for (i in seq_len(nrow(dir_mat))) {
    vals <- dir_mat[i, ]
    vals <- vals[!is.na(vals)]

    if (length(vals) < 2) {
        next
    }

    if (all(vals == 1) || all(vals == -1)) {
        class_vec[i] <- "concordant"
    } else if (any(vals == 1) && any(vals == -1)) {
        class_vec[i] <- "discordant"
    }
}

# extract gene sets
concordant_genes <- sum(class_vec == "concordant", na.rm = TRUE)
discordant_genes <- sum(class_vec == "discordant", na.rm = TRUE)

# summary
cat("Concordant genes:", concordant_genes, "\n")
cat("Discordant genes:", discordant_genes, "\n")


save.image("02_interpretation.RData")
