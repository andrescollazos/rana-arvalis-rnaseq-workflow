setwd(file.path(Sys.getenv("THESIS_DIR"), "code/04_de"))

# Load libraries
library(DESeq2)
library(pheatmap)

# Load workspace
load("dds.RData")

# -----------------------------------------------------------------------------
# 1. LRT for temperature-responsive genes
dds_LRT_temp <- DESeq(
    dds,
    test = "LRT",
    reduced = ~ family + stage + population
)
res_LRT_temp <- results(dds_LRT_temp)

# -----------------------------------------------------------------------------
# 2. LRT for interaction effects
dds_LRT_interaction <- DESeq(
    dds,
    test = "LRT",
    reduced = ~ family + stage + population + temperature
)
res_LRT_interaction <- results(dds_LRT_interaction)

# -----------------------------------------------------------------------------
# 3. Gene classification

padj_temp <- res_LRT_temp$padj
padj_interaction <- res_LRT_interaction$padj

temperature_responsive_no_interaction <- rownames(dds_LRT_temp)[
    !is.na(padj_temp) & padj_temp < 0.05 &
        (is.na(padj_interaction) | padj_interaction >= 0.05)
]

differential_plasticity <- rownames(dds_LRT_interaction)[
    !is.na(padj_interaction) & padj_interaction < 0.05
]

# -----------------------------------------------------------------------------
# 4. Post hoc interpretation: population-specific temperature effects

dds_Wald <- DESeq(dds, test = "Wald")
resultsNames(dds_Wald)

# Reference population = VF
# Reference temperature = 15
# Temperature effect being estimated = 20 vs 15

res_VF <- results(dds_Wald, name = "temperature_20_vs_15")

res_Ka <- results(
    dds_Wald,
    contrast = list(c("temperature_20_vs_15", "populationKa.temperature20"))
)

res_L <- results(
    dds_Wald,
    contrast = list(c("temperature_20_vs_15", "populationL.temperature20"))
)

res_E <- results(
    dds_Wald,
    contrast = list(c("temperature_20_vs_15", "populationE.temperature20"))
)

res_Upp <- results(
    dds_Wald,
    contrast = list(c("temperature_20_vs_15", "populationUpp.temperature20"))
)

res_CFin <- results(
    dds_Wald,
    contrast = list(c("temperature_20_vs_15", "populationC.Fin.temperature20"))
)

res_NA <- results(
    dds_Wald,
    contrast = list(c("temperature_20_vs_15", "populationNA.temperature20"))
)

res_NL <- results(
    dds_Wald,
    contrast = list(c("temperature_20_vs_15", "populationNL.temperature20"))
)

lfc_by_pop <- data.frame(
    "VF" = res_VF$log2FoldChange,
    "Ka" = res_Ka$log2FoldChange,
    "L" = res_L$log2FoldChange,
    "E" = res_E$log2FoldChange,
    "Upp" = res_Upp$log2FoldChange,
    "C.Fin" = res_CFin$log2FoldChange,
    "NA" = res_NA$log2FoldChange,
    "NL" = res_NL$log2FoldChange,
    row.names = rownames(dds_Wald)
)

padj_by_pop <- data.frame(
    "VF" = res_VF$padj,
    "Ka" = res_Ka$padj,
    "L" = res_L$padj,
    "E" = res_E$padj,
    "Upp" = res_Upp$padj,
    "C.Fin" = res_CFin$padj,
    "NA" = res_NA$padj,
    "NL" = res_NL$padj,
    row.names = rownames(dds_Wald)
)

# Subset by class for interpretation

shared_concordant <- apply(
    lfc_by_pop[temperature_responsive_no_interaction, , drop = FALSE],
    1,
    function(x) length(unique(sign(na.omit(x)))) == 1
)

# Subset to concordant genes
lfc_shared <- lfc_by_pop[
    temperature_responsive_no_interaction[shared_concordant], ,
    drop = FALSE
]
lfc_diff <- lfc_by_pop[differential_plasticity, , drop = FALSE]

padj_shared <- padj_by_pop[temperature_responsive_no_interaction[shared_concordant], , drop = FALSE]
padj_diff <- padj_by_pop[differential_plasticity, , drop = FALSE]

save.image("resultsConcordant.RData")
