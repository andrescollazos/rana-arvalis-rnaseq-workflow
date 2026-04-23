setwd(file.path(Sys.getenv("THESIS_DIR"), "code/06_enrichment"))
load("01_prepare_enrichment_MF.RData")

## ============================================
## ORA enrichment analysis
## - run one ORA per gene list
## - BH correction within each run
## - skip undersized lists
## ============================================

## Required packages
library(clusterProfiler)

## --------------------------------------------
## Inputs available from 01_prepare_enrichment.RData:
## - gene_lists_19
## - population_gene_lists
## - universe_temp_annotated_filtered
## - universe_interaction_annotated_filtered
## - TERM2GENE_temp_filtered          # columns: term, gene
## - TERM2GENE_interaction_filtered   # columns: term, gene
## - TERM2NAME                        # optional, if available
## --------------------------------------------

## --------------------------------------------
## Parameters
## --------------------------------------------

p_cutoff <- 0.05
q_cutoff <- 0.05
min_genes <- 10

## --------------------------------------------
## Helper function for ORA
## --------------------------------------------

run_ora <- function(genes, universe, TERM2GENE, TERM2NAME) {
    if (length(genes) < min_genes) {
        return(NULL)
    }

    enricher(
        gene = genes,
        universe = universe,
        TERM2GENE = TERM2GENE,
        TERM2NAME = TERM2NAME,
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1
    )
}

## --------------------------------------------
## Run ORA for all 19 lists
## --------------------------------------------

ora_results_MF <- list()

for (name in names(gene_lists_19)) {
    genes <- gene_lists_19[[name]]

    ## determine correct universe
    if (grepl("^concordant|^discordant", name)) {
        universe <- universe_interaction_annotated_filtered
        TERM2GENE_use <- TERM2GENE_interaction_filtered
    } else {
        universe <- universe_temp_annotated_filtered
        TERM2GENE_use <- TERM2GENE_temp_filtered
    }

    res <- run_ora(
        genes = genes,
        universe = universe,
        TERM2GENE = TERM2GENE_use,
        TERM2NAME = TERM2NAME
    )

    ora_results_MF[[name]] <- res
}

## --------------------------------------------
## Convert results to data frames
## --------------------------------------------

ora_tables_MF <- lapply(ora_results_MF, function(x) {
    if (is.null(x)) {
        return(NULL)
    }
    as.data.frame(x)
})

ora_tables_MF_sig <- lapply(ora_tables_MF, function(df) {
    if (is.null(df)) {
        return(NULL)
    }
    df[df$p.adjust < q_cutoff, ]
})

## --------------------------------------------
## Summary of significant terms
## --------------------------------------------

ora_summary_MF_all <- data.frame(
    list_name = names(ora_tables_MF),
    n_terms = vapply(ora_tables_MF, function(tab) {
        if (is.null(tab)) {
            return(NA_integer_)
        }
        nrow(tab)
    }, integer(1)),
    n_input_genes = vapply(gene_lists_19, length, integer(1)),
    row.names = NULL
)

ora_summary_MF_sig <- data.frame(
    list_name = names(ora_tables_MF_sig),
    n_terms = vapply(ora_tables_MF_sig, function(tab) {
        if (is.null(tab)) {
            return(NA_integer_)
        }
        nrow(tab)
    }, integer(1)),
    n_input_genes = vapply(gene_lists_19, length, integer(1)),
    row.names = NULL
)

ora_summary_MF_sig
ora_summary_MF_all


View(ora_tables_MF_sig[["concordant_down"]])
View(ora_tables_MF_sig[["concordant_down"]][, c("ID", "Description", "p.adjust")])


save.image("02_ORA_enrichment_MF.RData")
