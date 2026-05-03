rm(list = ls())
setwd(file.path(Sys.getenv("THESIS_DIR"), "code/06_enrichment"))
load("02_ORA_enrichment_BP_universe.RData")

# 1. Extract ORA tables
df_up <- ora_tables_BP_sig[["concordant_up"]]
df_down <- ora_tables_BP_sig[["concordant_down"]]

# 2. Load libraries
library(rrvgo)
library(GO.db) # required for semantic similarity
library(ggplot2)


# 3. Add score column
df_up$score <- -log10(df_up$p.adjust)
df_down$score <- -log10(df_down$p.adjust)

# 4. Compute semantic similarity (NO orgdb)
simMatrix_down <- calculateSimMatrix(
    df_down$ID,
    orgdb = "org.Hs.eg.db",
    ont = "BP",
    method = "Rel"
)

simMatrix_up <- calculateSimMatrix(
    df_up$ID,
    orgdb = "org.Hs.eg.db",
    ont = "BP",
    method = "Rel"
)

# 5. Reduce to non-redundant set (threshold = 0.7)
reduced_down <- reduceSimMatrix(
    simMatrix_down,
    scores = setNames(df_down$score, df_down$ID),
    threshold = 0.7,
    orgdb = "org.Hs.eg.db"
)

reduced_up <- reduceSimMatrix(
    simMatrix_up,
    scores = setNames(df_up$score, df_up$ID),
    threshold = 0.7,
    orgdb = "org.Hs.eg.db"
)

# 6. Extract representative terms
rep_down <- reduced_down[
    reduced_down$representative == reduced_down$term,
]

rep_up <- reduced_up[
    reduced_up$representative == reduced_up$term,
]

# 7. Merge with ORA stats
rep_down <- merge(rep_down, df_down, by.x = "term", by.y = "ID")
rep_up <- merge(rep_up, df_up, by.x = "term", by.y = "ID")

# 8. Visualization

# 8.1 Treemap — concordant_down (main plot)
pdf(file = "BP/treemap_concordant_down_BP_sig.pdf")
treemapPlot(reduced_down)
dev.off()

# 8.2 Treemap — concordant_up (check redundancy)
pdf(file = "BP/treemap_concordant_up_BP_sig.pdf")
treemapPlot(reduced_up)
dev.off()

# 8.3 Scatterplot — concordant_down
pdf(file = "BP/scatter_concordant_down_BP_sig.pdf")
scatterPlot(
    simMatrix_down,
    reduced_down
)
dev.off()

# 8.4 Scatterplot — concordant_up
pdf(file = "BP/scatter_concordant_up_BP_sig.pdf")
scatterPlot(
    simMatrix_up,
    reduced_up
)
dev.off()

# 8.5 Dotplots
library(enrichplot)

pdf(file = "BP/dotplots_concordant_BP_sig.pdf")
## concordant_down
res_cd <- ora_results_BP[["concordant_down"]]
res_cd@result <- res_cd@result[res_cd@result$p.adjust < 0.05, ]

res_cd_sim <- pairwise_termsim(res_cd)

dotplot(res_cd_sim, showCategory = 15) +
    ggtitle("Concordant genes (downregulated) — GO Biological Process enrichment") +
    theme(plot.title = element_text(hjust = 0.5))

## concordant_up
res_cu <- ora_results_BP[["concordant_up"]]
res_cu@result <- res_cu@result[res_cu@result$p.adjust < 0.05, ]

res_cu_sim <- pairwise_termsim(res_cu)

dotplot(res_cu_sim, showCategory = 15) +
    ggtitle("Concordant genes (upregulated) — GO Biological Process enrichment") +
    theme(plot.title = element_text(hjust = 0.5))

dev.off()


save.image("03_interpret_concordant.RData")
