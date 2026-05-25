setwd(file.path(Sys.getenv("THESIS_DIR"), "code/04_de"))
load("resultsInteraction.RData")
load("resultsTemperatureEffects.RData")

# Keep only interaction genes present in the contrast matrices
genes <- intersect(interaction_genes, rownames(lfc_mat))

lfc_int <- lfc_mat[genes, , drop = FALSE]
padj_int <- padj_mat[genes, , drop = FALSE]

# Significant population-specific temperature contrasts
sig_int <- !is.na(padj_int) & padj_int < 0.05

# Direction among significant contrasts only
up_int <- sig_int & lfc_int > 0
down_int <- sig_int & lfc_int < 0

# Summary per interaction gene
interaction_pattern_table <- data.frame(
    gene_id = genes,
    LRT_padj = res_LRT_interaction_df[genes, "padj"],
    n_sig_populations = rowSums(sig_int, na.rm = TRUE),
    n_up_populations = rowSums(up_int, na.rm = TRUE),
    n_down_populations = rowSums(down_int, na.rm = TRUE),
    mean_lfc = rowMeans(lfc_int, na.rm = TRUE),
    max_abs_lfc = apply(abs(lfc_int), 1, max, na.rm = TRUE)
)

# Which populations are significant for each gene
interaction_pattern_table$significant_populations <- apply(
    sig_int,
    1,
    function(x) paste(colnames(sig_int)[x], collapse = ";")
)

interaction_pattern_table$up_populations <- apply(
    up_int,
    1,
    function(x) paste(colnames(up_int)[x], collapse = ";")
)

interaction_pattern_table$down_populations <- apply(
    down_int,
    1,
    function(x) paste(colnames(down_int)[x], collapse = ";")
)

# Classify patterns
interaction_pattern_table$pattern <- "no_significant_population_contrast"

interaction_pattern_table$pattern[
    interaction_pattern_table$n_sig_populations == 1
] <- "single_population_response"

interaction_pattern_table$pattern[
    interaction_pattern_table$n_sig_populations >= 2 &
        interaction_pattern_table$n_sig_populations <= 3
] <- "small_subset_response"

interaction_pattern_table$pattern[
    interaction_pattern_table$n_sig_populations >= 4 &
        interaction_pattern_table$n_sig_populations <= 6
] <- "multi_population_response"

interaction_pattern_table$pattern[
    interaction_pattern_table$n_sig_populations >= 7
] <- "broad_response"

interaction_pattern_table$pattern[
    interaction_pattern_table$n_up_populations > 0 &
        interaction_pattern_table$n_down_populations > 0
] <- "opposite_direction_response"

# View main result
V(interaction_pattern_table)

# Count patterns
table(interaction_pattern_table$pattern)

# How many significant population-specific contrasts per interaction gene
table(interaction_pattern_table$n_sig_populations)

# Number of significant interaction-gene contrasts per population
population_contribution <- data.frame(
    population = colnames(sig_int),
    n_sig_genes = colSums(sig_int, na.rm = TRUE),
    n_up_genes = colSums(up_int, na.rm = TRUE),
    n_down_genes = colSums(down_int, na.rm = TRUE)
)

View(population_contribution)

# Para ordenar los genes donde la señal está más localizada:
localized_genes <- interaction_pattern_table[
    interaction_pattern_table$n_sig_populations <= 3,
]

localized_genes <- localized_genes[
    order(localized_genes$n_sig_populations, localized_genes$LRT_padj),
]

View(head(localized_genes))
