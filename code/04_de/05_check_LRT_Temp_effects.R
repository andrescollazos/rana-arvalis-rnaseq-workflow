rm(list = ls())
setwd(file.path(Sys.getenv("THESIS_DIR"), "code/04_de"))
load("resultsInteraction.RData")

# -----------------------------
# 1. Extract VF response from interaction model
# -----------------------------

resultsNames(dds_LRT_interaction)

res_vf_interaction <- results(
    dds_LRT_interaction,
    name = "temperature_20_vs_15"
)

vf_interaction <- data.frame(
    gene_id = rownames(res_vf_interaction),
    lfc_interaction = res_vf_interaction$log2FoldChange,
    se_interaction = res_vf_interaction$lfcSE,
    stat_interaction = res_vf_interaction$stat,
    pvalue_interaction = res_vf_interaction$pvalue,
    padj_interaction = res_vf_interaction$padj
)

View(vf_interaction)

# -----------------------------
# 2. Extract VF response from group model
# -----------------------------
load("resultsTemperatureEffects.RData")

resultsNames(dds_group)

res_vf_group <- results(
    dds_group,
    contrast = c("group", "VF_20", "VF_15")
)

vf_group <- data.frame(
    gene_id = rownames(res_vf_group),
    lfc_group = res_vf_group$log2FoldChange,
    se_group = res_vf_group$lfcSE,
    stat_group = res_vf_group$stat,
    pvalue_group = res_vf_group$pvalue,
    padj_group = res_vf_group$padj
)

# -----------------------------
# 3. Compare estimates
# -----------------------------

vf_compare <- merge(vf_interaction, vf_group, by = "gene_id")

vf_compare$lfc_diff <- vf_compare$lfc_interaction - vf_compare$lfc_group

summary(vf_compare$lfc_diff)

cor(
    vf_compare$lfc_interaction,
    vf_compare$lfc_group,
    use = "complete.obs"
)

library(ggplot2)

png("vf_temp_lfc_comparison.png", width = 800, height = 800, res = 100)
plot(
    vf_compare$lfc_group,
    vf_compare$lfc_interaction,
    xlab = "VF log2FC from group model",
    ylab = "VF log2FC from interaction model",
    main = "VF temperature response: group contrast vs interaction coefficient"
)
abline(0, 1, col = "red")
dev.off()


vf_compare_abs <- vf_compare[order(abs(vf_compare$lfc_diff), decreasing = TRUE), ]

head(vf_compare_abs[, c(
    "gene_id",
    "lfc_interaction",
    "lfc_group",
    "lfc_diff",
    "se_interaction",
    "se_group",
    "baseMean"
)], 20)


sum(abs(vf_compare$lfc_diff) > 1, na.rm = TRUE)
sum(abs(vf_compare$lfc_diff) > 5, na.rm = TRUE)
sum(abs(vf_compare$lfc_diff) > 10, na.rm = TRUE)

idx <- complete.cases(vf_compare[, c("lfc_interaction", "lfc_group")]) &
    abs(vf_compare$lfc_diff) < 1

cor(vf_compare$lfc_interaction[idx], vf_compare$lfc_group[idx])
