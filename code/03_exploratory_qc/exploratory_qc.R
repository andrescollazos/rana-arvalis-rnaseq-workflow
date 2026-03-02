setwd('~/Documents/Studies/Thesis/code/03_exploratory_qc')
# QC script

# ---------- Inputs ----------
counts_file <- "../../analyses/03_exploratory_qc/gene_counts.tsv"
meta_file   <- "../../analyses/03_exploratory_qc/metadata.tsv"

# ---------- Read data ----------
counts <- read.table(counts_file, header = TRUE, sep = "\t", check.names = FALSE,
                     row.names = 1, quote = "", comment.char = "")
counts <- as.matrix(counts)

meta <- read.table(meta_file,
                   header = TRUE,
                   sep = "\t",
                   check.names = FALSE,
                   quote = "",
                   comment.char = "",
                   na.strings = c("", "NaN", "NULL"),
                   stringsAsFactors = FALSE)
stopifnot("sample" %in% colnames(meta))

# Keep only samples present in both
common_samples <- intersect(colnames(counts), as.character(meta$sample))
if (length(common_samples) < 2) stop("Too few overlapping samples between counts and metadata.")
counts <- counts[, common_samples, drop = FALSE]
meta <- meta[match(common_samples, as.character(meta$sample)), , drop = FALSE]
rownames(meta) <- as.character(meta$sample)

# ---------- Library Size Plot (Before normalization) ----------

# Compute library sizes (total counts per sample)
library_sizes <- colSums(counts)

# Add library sizes to metadata
meta$library_size <- library_sizes[rownames(meta)]

# Summary statistics
summary(library_sizes)

# Extract vectors from existing meta
lib <- meta$library_size
pop <- meta$population
names(lib) <- meta$sample

# Order by population, then by library size
ord <- order(pop, lib)
lib <- lib[ord]
pop <- pop[ord]

# Convert to millions
lib_m <- lib / 1e6

# Median
median_reads <- median(lib_m)

# Colors by population
pop_levels <- unique(pop)
pop_colors <- setNames(rainbow(length(pop_levels)), pop_levels)
bar_cols <- pop_colors[pop]

pdf("../../analyses/03_exploratory_qc/total_read_counts_before_normalization_barplot.pdf",
    width = 12, height = 6)

par(mar = c(8, 5, 4, 2))

bp <- barplot(
  lib_m,
  col = bar_cols,
  border = NA,
  xaxt = "n",
  main = "Library Size per Sample (Before Normalization)",
  ylab = "Total Reads (Millions)"
)

axis(1, at = bp, labels = names(lib_m), las = 2, cex.axis = 0.6)

abline(h = median_reads, col = "black", lwd = 2, lty = 2)

legend(
  "topleft",
  legend = pop_levels,
  fill = pop_colors,
  border = NA,
  bty = "n",
  cex = 0.8
)

dev.off()


