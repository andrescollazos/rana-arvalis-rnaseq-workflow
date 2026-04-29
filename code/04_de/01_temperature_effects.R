setwd(file.path(Sys.getenv("THESIS_DIR"), "code/04_de"))

# Load libraries
library(DESeq2)

# Load workspace
load("dds.RData")

meta$group <- factor(
    paste(meta$population, meta$temperature, sep = "_")
)

# Set reference level
meta$group <- relevel(meta$group, ref = "VF_15")

# Rebuild DESeqDataSet with ~ group
dds_group <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~ family + stage + group
)

# Apply same low-count filter
keep <- rowSums(counts(dds_group) >= 10) >= (0.5 * min_group_size)
dds_group <- dds_group[keep, ]

# Run Wald model
dds_group <- DESeq(dds_group, test = "Wald")

# Inspect available coefficients if needed
resultsNames(dds_group)

# Initialize matrices
lfc_mat <- matrix(NA, nrow = nrow(dds_group), ncol = length(populations))
padj_mat <- matrix(NA, nrow = nrow(dds_group), ncol = length(populations))

rownames(lfc_mat) <- rownames(dds_group)
rownames(padj_mat) <- rownames(dds_group)
colnames(lfc_mat) <- populations
colnames(padj_mat) <- populations

# Fill matrices
for (pop in populations) {
    res <- results(
        dds_group,
        contrast = c("group", paste0(pop, "_20"), paste0(pop, "_15"))
    )

    lfc_mat[, pop] <- res$log2FoldChange
    padj_mat[, pop] <- res$padj
}

# Per-population significant gene lists
temp_genes_by_pop <- lapply(populations, function(pop) {
    rownames(padj_mat)[!is.na(padj_mat[, pop]) & padj_mat[, pop] < 0.05]
})
names(temp_genes_by_pop) <- populations

save.image("resultsTemperatureEffects.RData")

# Save one file per population
for (pop in populations) {
    write.table(
        data.frame(gene_id = temp_genes_by_pop[[pop]]),
        file = file.path(paste0("results/temperature_genes_", pop, ".tsv")),
        sep = "\t", quote = FALSE, row.names = FALSE
    )
}
