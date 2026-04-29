rm(list = ls())
setwd(file.path(Sys.getenv("THESIS_DIR"), "code/06_enrichment"))
load("02_ORA_enrichment_BP.RData")

library(enrichplot)
library(ggplot2)

pdf(file = "BP/dotplots_Ka_BP_sig.pdf")
## Ka_down
res_cd <- ora_results_BP[["Ka_down"]]
res_cd@result <- res_cd@result[res_cd@result$p.adjust < 0.05, ]

res_cd_sim <- pairwise_termsim(res_cd)

dotplot(res_cd_sim, showCategory = 20) +
    ggtitle("Ka genes (downregulated) — GO Biological Process enrichment") +
    theme(plot.title = element_text(hjust = 0.5))

## Ka_up
res_cu <- ora_results_BP[["Ka_up"]]
res_cu@result <- res_cu@result[res_cu@result$p.adjust < 0.05, ]

res_cu_sim <- pairwise_termsim(res_cu)

dotplot(res_cu_sim, showCategory = 20) +
    ggtitle("Ka genes (upregulated) — GO Biological Process enrichment") +
    theme(plot.title = element_text(hjust = 0.5))

dev.off()

rm(list = ls())

load("02_ORA_enrichment_MF.RData")

pdf(file = "MF/dotplots_Ka_MF_sig.pdf")
## Ka_down
res_cd <- ora_results_MF[["Ka_down"]]
res_cd@result <- res_cd@result[res_cd@result$p.adjust < 0.05, ]

dotplot(res_cd, showCategory = 20) +
    ggtitle("Ka genes (downregulated) — GO Molecular Function enrichment") +
    theme(plot.title = element_text(hjust = 0.5))

## Ka_up
res_cu <- ora_results_MF[["Ka_up"]]
res_cu@result <- res_cu@result[res_cu@result$p.adjust < 0.05, ]

dotplot(res_cu, showCategory = 20) +
    ggtitle("Ka genes (upregulated) — GO Molecular Function enrichment") +
    theme(plot.title = element_text(hjust = 0.5))
dev.off()

rm(list = ls())

load("02_ORA_enrichment_CC.RData")

pdf(file = "CC/dotplots_Ka_CC_sig.pdf")
## Ka_down
res_cd <- ora_results_CC[["Ka_down"]]
res_cd@result <- res_cd@result[res_cd@result$p.adjust < 0.05, ]

dotplot(res_cd, showCategory = 20) +
    ggtitle("Ka genes (downregulated) — GO Cellular Component enrichment") +
    theme(plot.title = element_text(hjust = 0.5))

## Ka_up
res_cu <- ora_results_CC[["Ka_up"]]
res_cu@result <- res_cu@result[res_cu@result$p.adjust < 0.05, ]

dotplot(res_cu, showCategory = 20) +
    ggtitle("Ka genes (upregulated) — GO Cellular Component enrichment") +
    theme(plot.title = element_text(hjust = 0.5))
dev.off()
