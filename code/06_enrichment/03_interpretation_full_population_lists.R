rm(list = ls())
setwd(file.path(Sys.getenv("THESIS_DIR"), "code/06_enrichment"))
load("02_ORA_enrichment_BP.RData")

library(rrvgo)
library(pheatmap)
library(dplyr)

full_up_names <- paste0(populations, "_full_up")
full_down_names <- paste0(populations, "_full_down")

# Build semantic recurrence matrix
build_semantic_recurrence <- function(
  ora_tables_sig,
  list_names,
  populations,
  ont = "BP",
  orgdb = "org.Hs.eg.db",
  threshold = 0.7
) {
    # -----------------------------
    # 1. Extract significant GO terms
    # -----------------------------
    term_df <- do.call(rbind, lapply(seq_along(list_names), function(i) {
        list_name <- list_names[i]
        pop <- populations[i]
        df <- ora_tables_sig[[list_name]]

        if (is.null(df) || nrow(df) == 0) {
            return(NULL)
        }

        data.frame(
            population = pop,
            list_name = list_name,
            GO = df$ID,
            term_name = df$Description,
            p.adjust = df$p.adjust,
            score = -log10(df$p.adjust),
            stringsAsFactors = FALSE
        )
    }))

    if (is.null(term_df) || nrow(term_df) == 0) {
        stop("No enriched terms found.")
    }

    # -----------------------------
    # 2. Semantic similarity
    # -----------------------------
    all_terms <- unique(term_df$GO)

    simMatrix <- calculateSimMatrix(
        all_terms,
        orgdb = orgdb,
        ont = ont,
        method = "Rel"
    )

    # Use best score per GO term across populations
    term_scores <- term_df %>%
        group_by(GO) %>%
        summarise(score = max(score, na.rm = TRUE), .groups = "drop")

    scores <- term_scores$score
    names(scores) <- term_scores$GO

    # -----------------------------
    # 3. Reduce semantic redundancy
    # -----------------------------
    reduced <- reduceSimMatrix(
        simMatrix = simMatrix,
        scores = scores,
        threshold = threshold,
        orgdb = orgdb
    )

    # reduced has:
    # term, parent, parentTerm, score, uniqueness, dispensability

    reduced_map <- data.frame(
        GO = rownames(reduced),
        representative_GO = reduced$parent,
        representative_name = reduced$parentTerm,
        stringsAsFactors = FALSE
    )

    term_df_reduced <- term_df %>%
        left_join(reduced_map, by = "GO")

    # -----------------------------
    # 4. Build population × semantic cluster matrix
    # -----------------------------
    cluster_scores <- term_df_reduced %>%
        group_by(representative_GO, representative_name, population) %>%
        summarise(
            score = max(score, na.rm = TRUE),
            n_terms_in_cluster = n_distinct(GO),
            .groups = "drop"
        )

    mat <- matrix(
        0,
        nrow = length(unique(cluster_scores$representative_GO)),
        ncol = length(populations)
    )

    rownames(mat) <- unique(cluster_scores$representative_GO)
    colnames(mat) <- populations

    for (i in seq_len(nrow(cluster_scores))) {
        mat[
            cluster_scores$representative_GO[i],
            cluster_scores$population[i]
        ] <- cluster_scores$score[i]
    }

    # Add readable row labels
    rep_names <- cluster_scores %>%
        distinct(representative_GO, representative_name)

    row_labels <- rep_names$representative_name
    names(row_labels) <- rep_names$representative_GO

    rownames(mat) <- row_labels[rownames(mat)]

    # -----------------------------
    # 5. Keep recurrent semantic clusters
    # -----------------------------
    keep <- rowSums(mat > 0) >= 2
    mat_recurrent <- mat[keep, , drop = FALSE]

    return(list(
        term_df = term_df,
        reduced = reduced,
        term_df_reduced = term_df_reduced,
        matrix_all = mat,
        matrix_recurrent = mat_recurrent
    ))
}

# Build UP semantic recurrence
bp_up_recurrence <- build_semantic_recurrence(
    ora_tables_sig = ora_tables_BP_sig,
    list_names = full_up_names,
    populations = populations,
    ont = "BP",
    orgdb = "org.Hs.eg.db",
    threshold = 0.7
)

# Build DOWN semantic recurrence
bp_down_recurrence <- build_semantic_recurrence(
    ora_tables_sig = ora_tables_BP_sig,
    list_names = full_down_names,
    populations = populations,
    ont = "BP",
    orgdb = "org.Hs.eg.db",
    threshold = 0.7
)

# Plot the semantic recurrence matrix
pdf("BP/recurrent_GO_bp_full_up.pdf", width = 12, height = 10)
p_up <- pheatmap(
    bp_up_recurrence$matrix_recurrent,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    main = "Recurrent BP modules among upregulated full lists",
    fontsize_row = 6,
    border_color = NA
)
print(p_up)
dev.off()

pdf("BP/recurrent_GO_bp_full_down.pdf", width = 12, height = 10)
p_down <- pheatmap(
    bp_down_recurrence$matrix_recurrent,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    main = "Recurrent BP modules among downregulated full lists",
    fontsize_row = 6,
    border_color = NA
)
print(p_down)
dev.off()

# inspect reduced tables
View(bp_up_recurrence$term_df_reduced)
View(bp_down_recurrence$term_df_reduced)
