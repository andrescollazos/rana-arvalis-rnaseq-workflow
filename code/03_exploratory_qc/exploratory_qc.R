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
meta$region <- ifelse(meta$population %in% c("C.Fin", "E", "L"), "East",
  ifelse(meta$population %in% c("Ka", "Upp"), "South",
    ifelse(meta$population %in% c("NA", "NL", "VF"), "North", NA)
  )
)
meta$region <- factor(meta$region)

# --------------------------------------------------------------
# 1. Library Size Plot (Before normalization)

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
# 2. Dataset and Low-count filtering
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
# 3. Gene detection rate per sample (post-filtering)

## 1) Detection counts and percentages (on filtered dds)
detected_genes <- colSums(counts(dds) > 0)
detected_pct <- 100 * detected_genes / nrow(dds)

## Store in metadata
meta$detected_genes <- detected_genes[rownames(meta)]
meta$detected_pct <- detected_pct[rownames(meta)]

## 2) Barplot: detected genes, ordered by population then within-population
pop <- meta$population
ord <- order(pop, meta$detected_genes, decreasing = TRUE)

## Colors by population (reuse your approach)
pop_levels <- unique(pop)
pop_colors <- setNames(rainbow(length(pop_levels)), pop_levels)
bar_cols <- pop_colors[as.character(pop[ord])]

pdf("../../analyses/03_exploratory_qc/gene_detection_rate.pdf", width = 18, height = 10)

par(mar = c(12, 5, 4, 10), xpd = NA)

bp_det_pct <- barplot(
  meta$detected_pct[ord],
  col = bar_cols,
  border = NA,
  xaxt = "n",
  ylab = "Detected genes (%)",
  main = "Percent gene detection per sample"
)

axis(
  side = 1,
  at = bp_det_pct,
  labels = rownames(meta)[ord],
  las = 2,
  cex.axis = 0.7,
  tick = FALSE
)

legend(
  "topleft",
  legend = names(pop_colors),
  fill = pop_colors,
  bty = "n",
  cex = 0.9
)

dev.off()

# -----------------------------------------------------------------------------
# NORMALIZATION


# -----------------------------------------------------------------------------
# 6. Density analysis

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

# Make x-order stable (use the current column order from the matrices)
vst_df$sample <- factor(vst_df$sample, levels = colnames(vst_mat))
raw_df$sample <- factor(raw_df$sample, levels = colnames(log_raw_mat))

p_vst_iqr <- ggplot(vst_df, aes(x = sample, y = expression)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.size = 1.5,
    fill = "#4F9BD4",
    color = "#2C6FA3",
    alpha = 0.6,
    width = 0.6
  ) +
  stat_boxplot(geom = "errorbar", width = 0.25, color = "#2C6FA3") +
  labs(x = "samples", y = "logcount") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.minor = element_blank()
  )

p_raw_iqr <- ggplot(raw_df, aes(x = sample, y = expression)) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.size = 1.5,
    fill = "#4F9BD4",
    color = "#2C6FA3",
    alpha = 0.6,
    width = 0.6
  ) +
  stat_boxplot(geom = "errorbar", width = 0.25, color = "#2C6FA3") +
  labs(x = "samples", y = "log2(count + 1)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.minor = element_blank()
  )

# One PDF with two pages (landscape helps with ~70 samples)
pdf("../../analyses/03_exploratory_qc/IQR_plots.pdf", width = 14, height = 8.5, onefile = TRUE)
print(p_vst_iqr)
print(p_raw_iqr)
dev.off()

# -----------------------------------------------------------------------------
# 7. Correlation analysis

source("correlation_analysis.R")

run_and_save <- function(vst_sub, meta_sub, prefix, plot_title,
                         outdir = "../../analyses/03_exploratory_qc/") {
  # BASELINE
  res_baseline <- run_correlation_analysis(
    vsd_mat = vst_sub,
    meta = meta_sub,
    gene_mode = "baseline",
    palette_mode = "red_white",
    plot_title = paste0(plot_title, " | Pearson correlation heatmap (All genes)")
  )

  res_baseline$plot_heatmap(
    file.path(outdir, paste0(prefix, "_baseline_heatmap.jpg"))
  )

  # TOP VARIABLE
  res_topvar <- run_correlation_analysis(
    vsd_mat = vst_sub,
    meta = meta_sub,
    gene_mode = "top_variable",
    scale_min = 0,
    scale_max = 1,
    palette_mode = "purple_yellow",
    plot_title = paste0(plot_title, " | Pearson correlation heatmap (Top variable genes)")
  )

  res_topvar$plot_heatmap(
    file.path(outdir, paste0(prefix, "_topvar_heatmap.jpg"))
  )

  invisible(NULL)
}
## All samples
run_and_save(vst_mat, meta, "all_samples", plot_title = "All samples")

## South and North (Swedish-only samples)
meta_SN <- meta[meta$region %in% c("South", "North"), , drop = FALSE]
vst_SN <- vst_mat[, rownames(meta_SN), drop = FALSE]
run_and_save(
  vst_SN,
  meta_SN,
  "swedish_samples",
  plot_title = "Swedish samples"
)

## Per population
populations <- c("NA", "NL", "VF", "Ka", "Upp", "C.Fin", "E", "L")

for (pop in populations) {
  meta_sub <- meta[meta$population == pop, , drop = FALSE]
  vst_sub <- vst_mat[, rownames(meta_sub), drop = FALSE]

  run_and_save(
    vst_sub,
    meta_sub,
    paste0("population_", pop),
    plot_title = paste0("Population ", pop)
  )
}

# -----------------------------------------------------------------------------
# 8. PCA

library(ggrepel)
source("PCA.R")

## 8.1 All samples

pca_baseline <- make_pca_plots(
  vsd_mat = vst_mat,
  meta = meta,
  gene_mode = "baseline",
  plot_title = "PCA (All filtered genes)"
)

ggsave(
  filename = "../../analyses/03_exploratory_qc/08_PCA_baseline.pdf",
  plot = gridExtra::marrangeGrob(
    grobs = pca_baseline$plots,
    nrow = 2,
    ncol = 1
  ),
  width = 10,
  height = 12
)

write.table(
  pca_baseline$pca_variance,
  "../../analyses/03_exploratory_qc/08_PCA_baseline_variance.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

## 8.2 Top variable genes

pca_topvar <- make_pca_plots(
  vsd_mat = vst_mat,
  meta = meta,
  gene_mode = "top_variable_100",
  plot_title = "PCA (Top 100 variable genes)"
)

ggsave(
  filename = "../../analyses/03_exploratory_qc/08_PCA_topvar_100.pdf",
  plot = gridExtra::marrangeGrob(
    grobs = pca_topvar$plots,
    nrow = 2,
    ncol = 1
  ),
  width = 10,
  height = 12
)

write.table(
  pca_topvar$pca_variance,
  "../../analyses/03_exploratory_qc/08_PCA_topvar_100_variance.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

## 8.3 Top variable genes

pca_top500 <- make_pca_plots(
  vsd_mat = vst_mat,
  meta = meta,
  gene_mode = "top_variable_500",
  plot_title = "PCA (Top 500 variable genes)"
)

ggsave(
  filename = "../../analyses/03_exploratory_qc/08_PCA_topvar_500.pdf",
  plot = gridExtra::marrangeGrob(
    grobs = pca_top500$plots,
    nrow = 2,
    ncol = 1
  ),
  width = 10,
  height = 12
)

write.table(
  pca_top500$pca_variance,
  "../../analyses/03_exploratory_qc/08_PCA_topvar_500_variance.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

# 8.4 Extraction date
pca_exdate <- make_pca_plots(
  vsd_mat = vst_mat,
  meta = meta,
  gene_mode = "top_variable_100",
  plot_title = "PCA (Top 100 variable genes)",
  mode = "extraction_date"
)

ggsave(
  filename = "../../analyses/03_exploratory_qc/08_PCA_exdate_100.pdf",
  plot = gridExtra::marrangeGrob(
    grobs = pca_exdate$plots,
    nrow = 2,
    ncol = 1
  ),
  width = 10,
  height = 12
)

# -----------------------------------------------------------------------------
# 9. Regression of Top PCs on biological and technical variables

source("PCA_regression.R")

reg_baseline <- run_pca_regression(pca_baseline, n_pcs = 5)
reg_top500 <- run_pca_regression(pca_top500, n_pcs = 5)
reg_top100 <- run_pca_regression(pca_topvar, n_pcs = 5)
