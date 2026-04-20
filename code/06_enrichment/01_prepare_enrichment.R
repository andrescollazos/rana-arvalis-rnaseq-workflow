setwd(file.path(Sys.getenv("THESIS_DIR"), "code/04_de"))
load("resultsInteraction.RData")
load("resultsTemperatureEffects.RData")
setwd(file.path(Sys.getenv("THESIS_DIR"), "code/05_GO_annotation"))
load("04_GO_terms_mapping.RData")
setwd(file.path(Sys.getenv("THESIS_DIR"), "code/06_enrichment"))

## ================================
## Final annotated gene universe
## for temperature-response tests
## ================================

## Inputs:
## - dds_group
## - loc_to_go
##
## Columns in loc_to_go:
## - gene_id   : original gene IDs (LOC_* / transcript-level gene IDs used in DE)
## - GO        : GO term IDs mapped via Xenopus tropicalis

## 1. Raw temperature-response universe:
## all genes that were actually tested in the temperature model
universe_temp_raw <- rownames(dds_group)

length(universe_temp_raw)

## 2. Clean mapping to unique gene-GO pairs only
loc_to_go_clean <- unique(loc_to_go[, c("gene_id", "GO")])

## 3. Keep only annotations for genes in the temperature universe
loc_to_go_temp <- loc_to_go_clean[
    loc_to_go_clean$gene_id %in% universe_temp_raw,
]

## 4. Final annotated temperature universe:
## genes tested in dds_group AND represented in GO mapping
universe_temp_annotated <- sort(unique(loc_to_go_temp$gene_id))

length(universe_temp_annotated)

## 5. Diagnostics
n_raw_tested <- length(universe_temp_raw)
n_annotated <- length(universe_temp_annotated)
n_not_annotated <- sum(!universe_temp_raw %in% universe_temp_annotated)

annotation_summary_temp <- data.frame(
    n_tested_genes = n_raw_tested,
    n_annotated_genes = n_annotated,
    n_unannotated_genes = n_not_annotated,
    prop_annotated = n_annotated / n_raw_tested
)

annotation_summary_temp

## 6. TERM2GENE restricted to the temperature universe
## This is the object to be used in enrichment for population up/down lists
TERM2GENE_temp <- unique(loc_to_go_temp[, c("GO", "gene_id")])
colnames(TERM2GENE_temp) <- c("term", "gene")

## 7. Term sizes within this universe
term_sizes_temp <- sort(table(TERM2GENE_temp$term), decreasing = TRUE)

term_sizes_temp_df <- data.frame(
    term = names(term_sizes_temp),
    n_genes = as.integer(term_sizes_temp),
    row.names = NULL
)

head(term_sizes_temp_df)

## 8. Filter very small terms before enrichment
min_size <- 10

valid_terms_temp <- term_sizes_temp_df$term[
    term_sizes_temp_df$n_genes >= min_size
]

TERM2GENE_temp_filtered <- TERM2GENE_temp[
    TERM2GENE_temp$term %in% valid_terms_temp,
]

## 9. Recompute final universe after term-size filtering
## Usually this becomes the practical enrichment universe
universe_temp_annotated_filtered <- sort(unique(TERM2GENE_temp_filtered$gene))

length(universe_temp_annotated_filtered)
