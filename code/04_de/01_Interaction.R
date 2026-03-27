setwd(file.path(Sys.getenv("THESIS_DIR"), "code/04_de"))

# -----------------------------------------------------------------------------
# Load libraries
library(DESeq2)
library(pheatmap)

# Load workspace
load("dds.RData")
# -----------------------------------------------------------------------------

# 1. LRT for interaction effects
dds_LRT_interaction <- DESeq(
    dds,
    test = "LRT",
    reduced = ~ family + stage + population + temperature
)
res_LRT_interaction <- results(dds_LRT_interaction)

# Convert to data frame
res_LRT_interaction_df <- as.data.frame(res_LRT_interaction)
res_LRT_interaction_df$gene_id <- rownames(res_LRT_interaction_df)

# Define significance filter
sig_idx <- !is.na(res_LRT_interaction_df$padj) &
    res_LRT_interaction_df$padj < 0.05

# Vector of gene IDs
interaction_genes <- res_LRT_interaction_df$gene_id[sig_idx]

# Full table of significant genes
res_LRT_interaction_sig <- res_LRT_interaction_df[sig_idx, ]

# -----------------------------------------------------------------------------

# Basic summary
cat("Tested genes:", nrow(res_LRT_interaction_df), "\n")
cat("Significant interaction genes (padj < 0.05):", length(interaction_genes), "\n")

# -----------------------------------------------------------------------------

save.image("resultsInteraction.RData")

# -----------------------------------------------------------------------------

write.table(
    res_LRT_interaction_df,
    file = file.path("LRT_interaction_all_genes.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
    res_LRT_interaction_sig,
    file = file.path("LRT_interaction_significant_genes.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
    data.frame(gene_id = interaction_genes),
    file = file.path("LRT_interaction_gene_ids.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)
