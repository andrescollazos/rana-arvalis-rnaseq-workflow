rm(list = ls())
setwd(file.path(Sys.getenv("THESIS_DIR"), "code/04_de"))
load("02_interpretation.RData")
setwd(file.path(Sys.getenv("THESIS_DIR"), "code/05_GO_annotation"))
load("04_GO_terms_mapping.RData")

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

# 3. Merge to avoid duplication
final_tables <- lapply(final_lists, function(ids) {
    df <- merge(
        data.frame(gene_id = ids),
        mapping_table,
        by = "gene_id",
        all.x = TRUE
    )
    df[!duplicated(df$gene_id), ]
})

View(final_tables$C.Fin_up)

for (name in names(final_tables)) {
    write.csv(final_tables[[name]], paste0(name, ".csv"), row.names = FALSE)
}
