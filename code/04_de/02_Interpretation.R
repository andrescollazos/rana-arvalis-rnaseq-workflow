setwd(file.path(Sys.getenv("THESIS_DIR"), "code/04_de"))

library(pheatmap)
library(ggplot2)
library(tidyr)
library(dplyr)

# Load workspace
load("resultsConcordant.RData")

# -----------------------------------------------------------------------------
# Heatmap for concordant shared plasticity

pheatmap(
    lfc_shared,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    clustering_method = "complete",
    show_rownames = FALSE,
    main = "Concordant shared plasticity (LFC)"
)

# -----------------------------------------------------------------------------
# Heatmap for differential plasticity
pdf("heatmap_diff.pdf")
pheatmap(
    lfc_diff,
    scale = "row",
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    clustering_method = "complete",
    show_rownames = FALSE,
    main = "Differential plasticity (LFC)"
)
dev.off()

# -----------------------------------------------------------------------------
# Prepare data for plotting (differential plasticity only)

df_diff <- lfc_diff %>%
    as.data.frame() %>%
    mutate(gene = rownames(.)) %>%
    pivot_longer(
        cols = -gene,
        names_to = "population",
        values_to = "LFC"
    )

# -----------------------------------------------------------------------------
# Define colors

population_colors <- c(
    "C.Fin" = "#FF33E4",
    "E" = "#B7FF5E",
    "Ka" = "#FFDE52",
    "L" = "#00CDFF",
    "NA" = "#6f9fe2ff",
    "NL" = "#59FFC0",
    "Upp" = "#FF7070",
    "VF" = "#6a3d9a"
)

# -----------------------------------------------------------------------------
# Violin plot
pdf("boxplot_diff_clipped.pdf")
ggplot(df_diff, aes(x = population, y = LFC, fill = population)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = population_colors) +
    coord_cartesian(ylim = c(-3, 3)) +
    theme_minimal() +
    labs(
        title = "Distribution of temperature effects (LFC) per population",
        x = "Population",
        y = "Log2 fold change (20 vs 15)"
    )
dev.off()

# -----------------------------------------------------------------------------
# Categorization based on differential plasticity genes

df_diff <- lfc_diff %>%
    as.data.frame() %>%
    mutate(gene = rownames(.))

# Function to classify each gene
classify_pattern <- function(x) {
    x <- na.omit(x)

    signs <- sign(x)
    unique_signs <- unique(signs)

    # Same direction across populations
    if (length(unique_signs) == 1) {
        if (unique_signs == 1) {
            return("Same direction (up)")
        } else if (unique_signs == -1) {
            return("Same direction (down)")
        }
    }

    # Sign change
    if (length(unique_signs) > 1) {
        return("Sign change")
    }

    return(NA)
}

# Apply classification
pattern_labels <- apply(df_diff[, -ncol(df_diff)], 1, classify_pattern)

# Add magnitude-based category: population-specific response
# Define "strong" as |LFC| > 1 (you can adjust threshold)
strong_threshold <- 1

population_specific <- apply(
    df_diff[, -ncol(df_diff)],
    1,
    function(x) sum(abs(na.omit(x)) > strong_threshold) <= 2
)

# Final combined classification
final_category <- ifelse(
    population_specific,
    "Population-specific response",
    pattern_labels
)

# Replace NA if any remain
final_category[is.na(final_category)] <- "Other"

# -----------------------------------------------------------------------------
# Count categories

category_counts <- table(final_category)

df_counts <- as.data.frame(category_counts)
colnames(df_counts) <- c("category", "count")

# -----------------------------------------------------------------------------
# Barplot (counts)

ggplot(df_counts, aes(x = category, y = count, fill = category)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(
        title = "Differential plasticity pattern categories",
        x = "Pattern category",
        y = "Number of genes"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
