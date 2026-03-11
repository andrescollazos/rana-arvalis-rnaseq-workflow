setwd(file.path(Sys.getenv("THESIS_DIR"), "code/04_de"))

# --------------------------------------------------------------
# Read data
## Inputs
counts_file <- file.path(
    Sys.getenv("THESIS_DIR"),
    "analyses/03_exploratory_qc/gene_counts.tsv"
)
meta_file <- file.path(
    Sys.getenv("THESIS_DIR"),
    "analyses/03_exploratory_qc/metadata.tsv"
)

## Read data
counts <- read.table(counts_file,
    header = TRUE, sep = "\t", check.names = FALSE,
    row.names = 1, quote = "", comment.char = ""
)
counts <- as.matrix(counts)

meta <- read.table(meta_file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    quote = "",
    comment.char = "",
    na.strings = c("", "NaN", "NULL"),
    stringsAsFactors = FALSE
)
stopifnot("sample" %in% colnames(meta))

## Keep only samples present in both
common_samples <- intersect(colnames(counts), as.character(meta$sample))
if (length(common_samples) < 2) stop("Too few overlapping samples between counts and metadata.")
counts <- counts[, common_samples, drop = FALSE]
meta <- meta[match(common_samples, as.character(meta$sample)), , drop = FALSE]
rownames(meta) <- as.character(meta$sample)
meta$family <- factor(meta$family)
meta$stage <- factor(meta$stage)
meta$population <- factor(meta$population)
meta$temperature <- factor(meta$temperature)

# -----------------------------------------------------------------------------
# 2. Dataset and Low-count filtering
library(DESeq2)

dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~ family + stage + population * temperature
)

min_group_size <- min(table(meta$population))
keep <- rowSums(counts(dds) >= 10) >= (0.5 * min_group_size)
dds <- dds[keep, ]

# -----------------------------------------------------------------------------
# 3. Define reference levels
# dds$population <- relevel(dds$population, ref = "VF")
# dds$temperature <- relevel(dds$temperature, ref = "15")

# -----------------------------------------------------------------------------
# 4. Run DESeq2
dds <- DESeq(dds)
resultsNames(dds)

# -----------------------------------------------------------------------------
# 5. Save data
save.image("dds2.RData")
