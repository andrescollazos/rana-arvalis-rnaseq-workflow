rm(list = ls())
setwd(file.path(Sys.getenv("THESIS_DIR"), "code/04_de"))
load("02_interpretation.RData")
setwd(file.path(Sys.getenv("THESIS_DIR"), "code/05_GO_annotation"))
load("04_GO_terms_mapping.RData")
setwd(file.path(Sys.getenv("THESIS_DIR"), "code/06_enrichment"))

# 1. Create mapping table
mapping_table <- unique(
    loc_to_xtrop_gene[, c("gene_id", "protein_id", "xtrop_gene_id")]
)

# 2. Create lists
alpha <- 0.05

# -----------------------------
# Population-specific full significant lists
# -----------------------------

population_full_lists <- list()

for (pop in populations) {
    sig_up <- rownames(lfc_mat)[
        !is.na(padj_mat[, pop]) &
            padj_mat[, pop] < alpha &
            !is.na(lfc_mat[, pop]) &
            lfc_mat[, pop] > 0
    ]

    sig_down <- rownames(lfc_mat)[
        !is.na(padj_mat[, pop]) &
            padj_mat[, pop] < alpha &
            !is.na(lfc_mat[, pop]) &
            lfc_mat[, pop] < 0
    ]

    population_full_lists[[paste0(pop, "_full_up")]] <- sig_up
    population_full_lists[[paste0(pop, "_full_down")]] <- sig_down
}

# -----------------------------
# Population-specific unique significant lists
# -----------------------------

population_unique_lists <- list()

for (pop in populations) {
    genes <- population_specific[[pop]]

    up_genes <- genes[
        !is.na(lfc_mat[genes, pop]) &
            lfc_mat[genes, pop] > 0
    ]

    down_genes <- genes[
        !is.na(lfc_mat[genes, pop]) &
            lfc_mat[genes, pop] < 0
    ]

    population_unique_lists[[paste0(pop, "_unique_up")]] <- up_genes
    population_unique_lists[[paste0(pop, "_unique_down")]] <- down_genes
}

# -----------------------------
# Concordant and discordant sets
# -----------------------------

concordant_ids <- names(class_vec)[
    !is.na(class_vec) & class_vec == "concordant"
]

concordant_up <- concordant_ids[
    apply(dir_mat[concordant_ids, , drop = FALSE], 1, function(x) {
        vals <- x[!is.na(x)]
        all(vals == 1)
    })
]

concordant_down <- concordant_ids[
    apply(dir_mat[concordant_ids, , drop = FALSE], 1, function(x) {
        vals <- x[!is.na(x)]
        all(vals == -1)
    })
]

discordant <- names(class_vec)[
    !is.na(class_vec) & class_vec == "discordant"
]

# -----------------------------
# Final ORA input lists
# -----------------------------

ora_input_lists <- c(
    list(
        concordant_up   = concordant_up,
        concordant_down = concordant_down,
        discordant      = discordant
    ),
    population_unique_lists,
    population_full_lists
)

sapply(ora_input_lists, length)

# -----------------------------
# Export tables with metadata
# -----------------------------

final_tables <- list()

for (list_name in names(ora_input_lists)) {
    ids <- ora_input_lists[[list_name]]

    df <- data.frame(
        gene_id = ids,
        list_name = list_name,
        stringsAsFactors = FALSE
    )

    # Add LFC and padj when applicable (population-based lists)
    pop_match <- sub("_(full|unique)_(up|down)$", "", list_name)

    if (pop_match %in% populations) {
        df$log2FC <- lfc_mat[df$gene_id, pop_match]
        df$padj <- padj_mat[df$gene_id, pop_match]
    } else {
        df$log2FC <- NA
        df$padj <- NA
    }

    # Merge mapping (keep many-to-many)
    df <- merge(
        df,
        mapping_table,
        by = "gene_id",
        all.x = TRUE
    )

    final_tables[[list_name]] <- df
}

sapply(final_tables, nrow)

# Create output directory
dir.create("ORA_input_tables", showWarnings = FALSE)

# Export each table separately
for (name in names(final_tables)) {
    write.csv(
        final_tables[[name]],
        file = file.path("ORA_input_tables", paste0(name, ".csv")),
        row.names = FALSE
    )
}

# Export all tables combined
combined_table <- do.call(rbind, final_tables)

write.csv(
    combined_table,
    file = "ORA_input_tables/all_ORA_inputs_combined.csv",
    row.names = FALSE
)

save.image("00_ora_input_tables.RData")
