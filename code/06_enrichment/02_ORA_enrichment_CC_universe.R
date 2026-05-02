rm(list = ls())
setwd(file.path(Sys.getenv("THESIS_DIR"), "code/06_enrichment"))
load("01_prepare_CC_universe.RData")

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
## - ora_input_lists
## - universe_final
## - TERM2GENE_all
## --------------------------------------------

## --------------------------------------------
## Parameters
## --------------------------------------------

p_cutoff <- 0.05
q_cutoff <- 0.05

## --------------------------------------------
## Helper function for ORA
## --------------------------------------------

ora_results_CC <- list()

for (name in names(ora_lists_final)) {
    genes <- ora_lists_final[[name]]

    res <- enricher(
        gene = genes,
        universe = universe_final,
        TERM2GENE = TERM2GENE_filtered,
        TERM2NAME = TERM2NAME,
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1
    )

    ora_results_CC[[name]] <- res
}
## --------------------------------------------
## Convert results to data frames
## --------------------------------------------

ora_tables_CC <- lapply(ora_results_CC, function(x) {
    if (is.null(x)) {
        return(NULL)
    }
    as.data.frame(x)
})

## --------------------------------------------
## Significant terms
## --------------------------------------------

ora_tables_CC_sig <- lapply(ora_tables_CC, function(df) {
    if (is.null(df)) {
        return(NULL)
    }
    df[df$p.adjust < q_cutoff, ]
})

## --------------------------------------------
## Summary tables
## --------------------------------------------

ora_summary_CC_all <- data.frame(
    list_name = names(ora_tables_CC),
    n_terms = vapply(ora_tables_CC, function(tab) {
        if (is.null(tab)) {
            return(NA_integer_)
        }
        nrow(tab)
    }, integer(1)),
    n_input_genes = vapply(ora_lists_final, length, integer(1)),
    row.names = NULL
)

ora_summary_CC_sig <- data.frame(
    list_name = names(ora_tables_CC_sig),
    n_terms = vapply(ora_tables_CC_sig, function(tab) {
        if (is.null(tab)) {
            return(NA_integer_)
        }
        nrow(tab)
    }, integer(1)),
    n_input_genes = vapply(ora_lists_final, length, integer(1)),
    row.names = NULL
)

ora_summary_CC_sig
ora_summary_CC_all

names(ora_tables_CC_sig)
View(ora_tables_CC_sig$Ka_full_up)

View(ora_tables_CC_sig[["VF_down"]][, c("ID", "Description", "p.adjust")])


save.image("02_ORA_enrichment_CC_universe.RData")
