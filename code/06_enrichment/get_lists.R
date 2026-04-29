rm(list = ls())
setwd(file.path(Sys.getenv("THESIS_DIR"), "code/04_de"))
load("02_interpretation.RData")
setwd(file.path(Sys.getenv("THESIS_DIR"), "code/05_GO_annotation"))
load("04_GO_terms_mapping.RData")

# 1. Construir tabla de mapping (puede haber múltiples por gen)
mapping_table <- unique(
    loc_to_xtrop_gene[, c("gene_id", "protein_id", "xtrop_gene_id")]
)

# 2. Crear listas (IDs crudos, como ya los defines)
alpha <- 0.05

population_lists <- list()

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

    population_lists[[paste0(pop, "_up")]] <- sig_up
    population_lists[[paste0(pop, "_down")]] <- sig_down
}

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

final_lists <- c(
    population_lists,
    list(
        concordant_up   = concordant_up,
        concordant_down = concordant_down,
        discordant      = discordant
    )
)

# -----------------------------
# 3. Merge SIN duplicar
# -----------------------------
final_tables <- lapply(final_lists, function(ids) {
    df <- merge(
        data.frame(gene_id = ids),
        mapping_table,
        by = "gene_id",
        all.x = TRUE
    )
    df[!duplicated(df$gene_id), ]
})

for (name in names(final_tables)) {
    write.csv(final_tables[[name]], paste0(name, ".csv"), row.names = FALSE)
}
