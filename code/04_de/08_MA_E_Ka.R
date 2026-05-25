rm(list = ls())
setwd(file.path(Sys.getenv("THESIS_DIR"), "code/04_de"))
load("resultsTemperatureEffects.RData")

library(DESeq2)
library(ashr)
library(ggplot2)

res_list <- list()
res_shrunk_list <- list()

for (pop in populations) {
    group_20 <- paste0(pop, "_20")
    group_15 <- paste0(pop, "_15")

    res <- results(
        dds_group,
        contrast = c("group", group_20, group_15),
        alpha = 0.05
    )

    res_shrunk <- lfcShrink(
        dds_group,
        contrast = c("group", group_20, group_15),
        res = res,
        type = "ashr"
    )

    res_list[[pop]] <- res
    res_shrunk_list[[pop]] <- res_shrunk
}


pdf("Volcano_all_populations_shrunk.pdf", width = 7, height = 6)

for (pop in populations) {
    res_df <- as.data.frame(res_shrunk_list[[pop]])
    res_df$gene <- rownames(res_df)

    res_df <- res_df[
        !is.na(res_df$padj) &
            !is.na(res_df$log2FoldChange),
    ]

    res_df$significance <- "Not significant"

    res_df$significance[
        res_df$padj < 0.05 &
            res_df$log2FoldChange > 0
    ] <- "Up"

    res_df$significance[
        res_df$padj < 0.05 &
            res_df$log2FoldChange < 0
    ] <- "Down"

    p <- ggplot(
        res_df,
        aes(x = log2FoldChange, y = -log10(padj))
    ) +
        geom_point(
            aes(color = significance),
            alpha = 0.6,
            size = 1
        ) +
        geom_hline(
            yintercept = -log10(0.05),
            linetype = "dashed"
        ) +
        geom_vline(
            xintercept = 0,
            linetype = "dashed"
        ) +
        labs(
            title = paste0(pop, ": 20 vs 15"),
            x = "Shrunken log2 fold change",
            y = "-log10 adjusted p-value",
            color = "Direction"
        ) +
        theme_bw()

    print(p)
}

dev.off()
