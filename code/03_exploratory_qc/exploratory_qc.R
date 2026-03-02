setwd("~/Documents/Studies/Thesis/code/03_exploratory_qc")

# --------------------------------------------------------------
# Read data
## Inputs
counts_file <- "../../analyses/03_exploratory_qc/gene_counts.tsv"
meta_file <- "../../analyses/03_exploratory_qc/metadata.tsv"

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
meta$extraction_date <- factor(meta$extraction_date)
meta$family <- factor(meta$family)
meta$stage <- factor(meta$stage)
meta$population <- factor(meta$population)
meta$temperature <- factor(meta$temperature)

# --------------------------------------------------------------
# Library Size Plot (Before normalization)

## Compute library sizes (total counts per sample)
library_sizes <- colSums(counts)

## Add library sizes to metadata
meta$library_size <- library_sizes[rownames(meta)]

## Summary statistics
summary(library_sizes) # min, 1Q, median, mean, 3Q, max
sd(library_sizes) # standard deviation
mean(library_sizes) # mean
median(library_sizes) # median
var(library_sizes) # variance

## Extract vectors from existing meta
lib <- meta$library_size
pop <- meta$population
names(lib) <- meta$sample

## Order by population, then by library size
ord <- order(pop, lib)
lib <- lib[ord]
pop <- pop[ord]

## Convert to millions
lib_m <- lib / 1e6

## Median
median_reads <- median(lib_m)

## Colors by population
pop_levels <- unique(pop)
pop_colors <- setNames(rainbow(length(pop_levels)), pop_levels)
bar_cols <- pop_colors[pop]

pdf("../../analyses/03_exploratory_qc/total_read_counts_before_normalization_barplot.pdf",
  width = 12, height = 6
)

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

## 50% of median threshold (dotted)
half_median_reads <- 0.5 * median_reads
abline(h = half_median_reads, col = "black", lwd = 2, lty = 3)

legend(
  "topleft",
  legend = pop_levels,
  fill = pop_colors,
  border = NA,
  bty = "n",
  cex = 0.8
)

## Labels on the right side
text(
  x = par("usr")[2], y = median_reads,
  labels = paste0("Median (", round(median_reads, 2), "M)"),
  pos = 2, cex = 0.8
)

text(
  x = par("usr")[2], y = half_median_reads,
  labels = paste0("50% Median (", round(half_median_reads, 2), "M)"),
  pos = 2, cex = 0.8
)

dev.off()

# -----------------------------------------------------------------------------
# Dataset and Low-count filtering
library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta,
  design = ~ extraction_date + family + stage + population * temperature
)

min_group_size <- min(table(meta$population))
keep <- rowSums(counts(dds) >= 10) >= (0.5 * min_group_size)
dds <- dds[keep, ]

# Log-count distribution before vs after filtering

## Log10 transform
log_counts_before <- log10(counts + 1)
log_counts_after <- log10(counts(dds) + 1)

## Remove infinite / NA values
vals_before <- as.vector(log_counts_before)
vals_after <- as.vector(log_counts_after)

vals_before <- vals_before[is.finite(vals_before)]
vals_after <- vals_after[is.finite(vals_after)]

## Common x-limits
xlim_range <- range(c(vals_before, vals_after))

## Fixed breaks
breaks_seq <- seq(xlim_range[1], xlim_range[2], length.out = 101)

pdf("../../analyses/03_exploratory_qc/log_count_distribution_before_after_filtering.pdf",
  width = 12, height = 6
)

par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

hist(
  vals_before,
  breaks = breaks_seq,
  freq = FALSE,
  xlim = xlim_range,
  main = "Log10 Gene Counts Before Filtering",
  xlab = "log10(count + 1)",
  ylab = "Density",
  col = "steelblue",
  border = NA
)

hist(
  vals_after,
  breaks = breaks_seq,
  freq = FALSE,
  xlim = xlim_range,
  main = "Log10 Gene Counts After Filtering",
  xlab = "log10(count + 1)",
  ylab = "Density",
  col = "steelblue",
  border = NA
)

dev.off()

# -----------------------------------------------------------------------------
# DENSITY ANALYSIS

# Variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# Extract VST matrix
vst_mat <- assay(vsd)

# Convert to long format
library(reshape2)
library(ggplot2)
library(plotly)
library(htmlwidgets)

vst_df <- melt(vst_mat)
colnames(vst_df) <- c("gene", "sample", "expression")

p <- ggplot(vst_df, aes(x = expression, color = sample, fill = sample)) +
  geom_density(alpha = 0.3) +
  labs(x = "logcount", y = "Density") +
  theme_minimal()

# Convert to interactive plot
p_interactive <- ggplotly(p)

# Save as HTML
saveWidget(p_interactive, "vst_density_plot.html", selfcontained = TRUE)

# Log2 transform raw counts
log_raw_mat <- log2(counts + 1)

# Convert to long format
library(reshape2)
raw_df <- melt(log_raw_mat)
colnames(raw_df) <- c("gene", "sample", "expression")

# Density plot
p_raw <- ggplot(raw_df, aes(x = expression, color = sample, fill = sample)) +
  geom_density(alpha = 0.3) +
  labs(x = "log2(count + 1)", y = "Density") +
  theme_minimal()

# Interactive version
p_raw_interactive <- ggplotly(p_raw)

# Save as HTML
saveWidget(p_raw_interactive, "raw_density_plot.html", selfcontained = TRUE)

# -----------------------------------------------------------------------------
