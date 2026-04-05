setwd(file.path(Sys.getenv("THESIS_DIR"), "code/04_de"))

# -----------------------------------------------------------------------------
# Load libraries
library(DESeq2)
library(pheatmap)

# Load workspace
load("dds.RData")


# -----------------------------------------------------------------------------
# Population groupings
# -----------------------------------------------------------------------------

# Latitude grouping
# North: NA, NL, VF, C.Fin
# South: Upp, Ka, E, L
meta$lat_group <- NA_character_
meta$lat_group[meta$population %in% c("NA", "NL", "VF", "C.Fin")] <- "North"
meta$lat_group[meta$population %in% c("Upp", "Ka", "E", "L")] <- "South"
meta$lat_group <- factor(meta$lat_group, levels = c("North", "South"))

# Lineage grouping
# North: NA, NL, VF
# East: C.Fin, E, L
# South: Upp, Ka
meta$lineage <- NA_character_
meta$lineage[meta$population %in% c("NA", "NL", "VF")] <- "North"
meta$lineage[meta$population %in% c("C.Fin", "E", "L")] <- "East"
meta$lineage[meta$population %in% c("Upp", "Ka")] <- "South"
meta$lineage <- factor(meta$lineage, levels = c("North", "East", "South"))

# -----------------------------------------------------------------------------
# Subset: 15C
# -----------------------------------------------------------------------------

meta_15 <- meta[meta$temperature == "15", , drop = FALSE]
counts_15 <- counts[, meta_15$sample, drop = FALSE]

# -----------------------------------------------------------------------------
# Latitude model at 15C
# -----------------------------------------------------------------------------

dds_lat_15 <- DESeqDataSetFromMatrix(
    countData = counts_15,
    colData = meta_15,
    design = ~lat_group
)

dds_lat_15 <- DESeq(dds_lat_15)

# North vs South at 15C
res_lat_15_North_vs_South <- results(
    dds_lat_15,
    contrast = c("lat_group", "North", "South")
)

# -----------------------------------------------------------------------------
# Lineage model at 15C
# -----------------------------------------------------------------------------

dds_lin_15 <- DESeqDataSetFromMatrix(
    countData = counts_15,
    colData = meta_15,
    design = ~lineage
)

dds_lin_15 <- DESeq(dds_lin_15)

# North vs East at 15C
res_lin_15_North_vs_East <- results(
    dds_lin_15,
    contrast = c("lineage", "North", "East")
)

# North vs South at 15C
res_lin_15_North_vs_South <- results(
    dds_lin_15,
    contrast = c("lineage", "North", "South")
)

# East vs South at 15C
res_lin_15_East_vs_South <- results(
    dds_lin_15,
    contrast = c("lineage", "East", "South")
)

# -----------------------------------------------------------------------------
# Subset: 20C
# -----------------------------------------------------------------------------

meta_20 <- meta[meta$temperature == "20", , drop = FALSE]
counts_20 <- counts[, meta_20$sample, drop = FALSE]

# -----------------------------------------------------------------------------
# Latitude model at 20C
# -----------------------------------------------------------------------------

dds_lat_20 <- DESeqDataSetFromMatrix(
    countData = counts_20,
    colData = meta_20,
    design = ~lat_group
)

dds_lat_20 <- DESeq(dds_lat_20)

# North vs South at 20C
res_lat_20_North_vs_South <- results(
    dds_lat_20,
    contrast = c("lat_group", "North", "South")
)

# -----------------------------------------------------------------------------
# Lineage model at 20C
# -----------------------------------------------------------------------------

dds_lin_20 <- DESeqDataSetFromMatrix(
    countData = counts_20,
    colData = meta_20,
    design = ~lineage
)

dds_lin_20 <- DESeq(dds_lin_20)

# North vs East at 20C
res_lin_20_North_vs_East <- results(
    dds_lin_20,
    contrast = c("lineage", "North", "East")
)

# North vs South at 20C
res_lin_20_North_vs_South <- results(
    dds_lin_20,
    contrast = c("lineage", "North", "South")
)

# East vs South at 20C
res_lin_20_East_vs_South <- results(
    dds_lin_20,
    contrast = c("lineage", "East", "South")
)


save.image("resultsPopulationEffects.RData")
