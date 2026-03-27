if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("ExploreModelMatrix")

library(ExploreModelMatrix)
library(ggplot2)

pdf("design_matrix.pdf", width = 10, height = 8)

VisualizeDesign(
    sampleData = as.data.frame(colData(dds)),
    designFormula = ~ population * temperature
)

dev.off()
