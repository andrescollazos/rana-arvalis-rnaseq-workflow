rm(list = ls())
setwd(file.path(Sys.getenv("THESIS_DIR"), "code/04_de"))
load("02_interpretation.RData")
setwd(file.path(Sys.getenv("THESIS_DIR"), "code/05_GO_annotation"))
load("04_GO_terms_mapping.RData")
setwd(file.path(Sys.getenv("THESIS_DIR"), "code/06_enrichment"))
load("00_ora_input_tables.RData")

## ================================
## Final annotated gene universe
## for temperature-response tests
## ================================

## Inputs:
## - dds
## - loc_to_go_MF
##
## Columns in loc_to_go_MF:
## - gene_id   : original gene IDs (LOC_* / transcript-level gene IDs used in DE)
## - GO        : GO term IDs mapped via Xenopus tropicalis
## - ONTOLOGY  : MF

# Use only MF annotations
loc_to_go <- loc_to_go_MF

# -----------------------------
# 1. Define unified universe
# -----------------------------
universe_all_raw <- rownames(dds)

length(universe_all_raw)

## 2. Clean mapping to unique gene-GO pairs only
loc_to_go_clean <- unique(loc_to_go[, c("gene_id", "GO")])

## 3. Restrict to tested genes
loc_to_go_all <- loc_to_go_clean[
    loc_to_go_clean$gene_id %in% universe_all_raw,
]

# TERM2GENE
TERM2GENE_all <- unique(loc_to_go_all[, c("GO", "gene_id")])
colnames(TERM2GENE_all) <- c("term", "gene")

# -----------------------------
# 2. Term-size filtering
# -----------------------------
term_sizes <- table(TERM2GENE_all$term)

min_size <- 10
valid_terms <- names(term_sizes)[term_sizes >= min_size]

TERM2GENE_filtered <- TERM2GENE_all[
    TERM2GENE_all$term %in% valid_terms,
]

# Final annotated universe
universe_final <- sort(unique(TERM2GENE_filtered$gene))

# -----------------------------
# 3. Build ORA-ready lists WITH size filter
# -----------------------------

min_genes <- 20

# Step 1: intersect with universe (your original object)
ora_lists_filtered <- lapply(ora_input_lists, function(gene_vec) {
    intersect(gene_vec, universe_final)
})

# Step 2: diagnostics
ora_list_sizes <- data.frame(
    n_raw = sapply(ora_input_lists, length),
    n_annotated = sapply(ora_lists_filtered, length)
)

ora_list_sizes$prop_retained <- (
    ora_list_sizes$n_annotated / ora_list_sizes$n_raw
)

# View(ora_list_sizes)

# Step 3: apply min_genes filter
ora_lists_final <- ora_lists_filtered[
    sapply(ora_lists_filtered, length) >= min_genes
]

# Step 4: update summary to match final lists
ora_list_sizes <- ora_list_sizes[
    rownames(ora_list_sizes) %in% names(ora_lists_final),
]

View(ora_list_sizes)

save.image("01_prepare_MF_universe.RData")
