setwd(file.path(Sys.getenv("THESIS_DIR"), "code/04_de"))
load("02_interpretation.RData")
setwd(file.path(Sys.getenv("THESIS_DIR"), "code/05_GO_annotation"))
load("04_GO_terms_mapping.RData")
setwd(file.path(Sys.getenv("THESIS_DIR"), "code/06_enrichment"))

## ================================
## Final annotated gene universe
## for temperature-response tests
## ================================

## Inputs:
## - dds_group
## - loc_to_go_BP
##
## Columns in loc_to_go_BP:
## - gene_id   : original gene IDs (LOC_* / transcript-level gene IDs used in DE)
## - GO        : GO term IDs mapped via Xenopus tropicalis
## - ONTOLOGY  : BP

# Use only BP annotations
loc_to_go <- loc_to_go_BP

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

## ============================================
## Final annotated gene universe
## for concordant / discordant tests
## ============================================

## Inputs:
## - interaction_genes (from LRT_interaction)
## - loc_to_go
##
## Columns in loc_to_go:
## - gene_id : original gene IDs (LOC_*)
## - GO      : GO term IDs

## 1. Raw interaction universe:
## all genes significant in the LRT (G×E signal)
universe_interaction_raw <- interaction_genes

length(universe_interaction_raw)

## 2. Clean mapping to unique gene-GO pairs
loc_to_go_clean <- unique(loc_to_go[, c("gene_id", "GO")])

## 3. Keep only annotations for interaction genes
loc_to_go_interaction <- loc_to_go_clean[
    loc_to_go_clean$gene_id %in% universe_interaction_raw,
]

## 4. Final annotated interaction universe:
## interaction genes AND represented in GO mapping
universe_interaction_annotated <- sort(unique(loc_to_go_interaction$gene_id))

length(universe_interaction_annotated)

## 5. Diagnostics
n_raw_interaction <- length(universe_interaction_raw)
n_annotated_interaction <- length(universe_interaction_annotated)
n_not_annotated_interaction <- sum(
    !universe_interaction_raw %in% universe_interaction_annotated
)

annotation_summary_interaction <- data.frame(
    n_interaction_genes = n_raw_interaction,
    n_annotated_genes = n_annotated_interaction,
    n_unannotated_genes = n_not_annotated_interaction,
    prop_annotated = n_annotated_interaction / n_raw_interaction
)

annotation_summary_interaction

## 6. TERM2GENE restricted to interaction universe
## This is used for concordant / discordant enrichment
TERM2GENE_interaction <- unique(
    loc_to_go_interaction[, c("GO", "gene_id")]
)
colnames(TERM2GENE_interaction) <- c("term", "gene")

## 7. Term sizes within interaction universe
term_sizes_interaction <- sort(
    table(TERM2GENE_interaction$term),
    decreasing = TRUE
)

term_sizes_interaction_df <- data.frame(
    term = names(term_sizes_interaction),
    n_genes = as.integer(term_sizes_interaction),
    row.names = NULL
)

head(term_sizes_interaction_df)

## 8. Filter very small terms (stability)
min_size <- 10

valid_terms_interaction <- term_sizes_interaction_df$term[
    term_sizes_interaction_df$n_genes >= min_size
]

TERM2GENE_interaction_filtered <- TERM2GENE_interaction[
    TERM2GENE_interaction$term %in% valid_terms_interaction,
]

## 9. Final filtered interaction universe
universe_interaction_annotated_filtered <- sort(
    unique(TERM2GENE_interaction_filtered$gene)
)

length(universe_interaction_annotated_filtered)

## ============================================
## Build the 19 enrichment gene lists
## in the same ID space as loc_to_go
## ============================================

## Inputs assumed available:
## - populations
## - lfc_mat
## - padj_mat
## - interaction_genes
## - class_vec
## - universe_temp_annotated_filtered
## - universe_interaction_annotated_filtered

alpha <- 0.05

## --------------------------------------------
## 1. Per-population up/down gene lists
##    Universe: temperature universe
## --------------------------------------------

population_gene_lists <- list()

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

    ## restrict to annotated filtered temperature universe
    sig_up <- sort(intersect(sig_up, universe_temp_annotated_filtered))
    sig_down <- sort(intersect(sig_down, universe_temp_annotated_filtered))

    population_gene_lists[[paste0(pop, "_up")]] <- sig_up
    population_gene_lists[[paste0(pop, "_down")]] <- sig_down
}

## --------------------------------------------
## 2. Concordant up/down gene lists
##    Universe: interaction universe
## --------------------------------------------

concordant_ids <- names(na.omit(class_vec[class_vec == "concordant"]))
length(concordant_ids)

concordant_up <- concordant_ids[
    apply(dir_mat[concordant_ids, , drop = FALSE], 1, function(x) {
        vals <- x[!is.na(x)]
        all(vals == 1)
    })
]
length(concordant_up)

concordant_up <- sort(intersect(
    concordant_up,
    universe_interaction_annotated_filtered
))

length(concordant_up)

concordant_down <- concordant_ids[
    apply(dir_mat[concordant_ids, , drop = FALSE], 1, function(x) {
        vals <- x[!is.na(x)]
        all(vals == -1)
    })
]
length(concordant_down)

concordant_down <- sort(intersect(
    concordant_down,
    universe_interaction_annotated_filtered
))

length(concordant_down)

## --------------------------------------------
## 3. Discordant gene list
##    Universe: interaction universe
## --------------------------------------------

discordant <- names(class_vec)[
    !is.na(class_vec) & class_vec == "discordant"
]
length(discordant)
discordant <- sort(intersect(
    discordant,
    universe_interaction_annotated_filtered
))
length(discordant)

## --------------------------------------------
## 4. Collect all 19 lists
## --------------------------------------------

gene_lists_19 <- c(
    population_gene_lists,
    list(
        concordant_up = concordant_up,
        concordant_down = concordant_down,
        discordant = discordant
    )
)

length(gene_lists_19) # should be 19
names(gene_lists_19)

## --------------------------------------------
## 5. Optional summary table
## --------------------------------------------

gene_list_sizes <- data.frame(
    list_name = names(gene_lists_19),
    n_genes = vapply(gene_lists_19, length, integer(1)),
    universe = c(
        rep("temp", length(population_gene_lists)),
        "interaction", "interaction", "interaction"
    ),
    row.names = NULL
)

View(gene_list_sizes)


save.image("01_prepare_enrichment_BP.RData")
