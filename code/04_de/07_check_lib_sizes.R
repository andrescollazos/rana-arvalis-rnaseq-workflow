setwd(file.path(Sys.getenv("THESIS_DIR"), "code/04_de"))
load("resultsTemperatureEffects.RData")

## Filtered counts
filtered_counts <- counts(dds_group)

## Compute library sizes after filtering
library_sizes_filtered <- colSums(filtered_counts)

## Add filtered library sizes to metadata
meta$library_size_filtered <- library_sizes_filtered[rownames(meta)]

## Summary statistics
summary(library_sizes_filtered)
sd(library_sizes_filtered)
mean(library_sizes_filtered)
median(library_sizes_filtered)
var(library_sizes_filtered)

## Extract vectors from existing meta
lib <- meta$library_size_filtered
pop <- meta$population
names(lib) <- meta$sample

## Order by population, then by filtered library size
ord <- order(pop, lib)
lib <- lib[ord]
pop <- pop[ord]

## Convert to millions
lib_m <- lib / 1e6

## Median
median_reads <- median(lib_m)

## Colors by population
population <- c(
    "C.Fin" = "#FF33E4",
    "E" = "#B7FF5E",
    "Ka" = "#FFDE52",
    "L" = "#00CDFF",
    "NA" = "#6f9fe2ff",
    "NL" = "#59FFC0",
    "Upp" = "#FF7070",
    "VF" = "#6a3d9a"
)

bar_cols <- population[pop]

pdf(
    "results/total_read_counts_after_filtering_barplot.pdf",
    width = 12,
    height = 6
)

par(mar = c(8, 5, 4, 2))

bp <- barplot(
    lib_m,
    col = bar_cols,
    border = NA,
    xaxt = "n",
    main = "Library Size per Sample (After Low-count Filtering)",
    ylab = "Total Counts in Retained Genes (Millions)"
)

axis(1, at = bp, labels = names(lib_m), las = 2, cex.axis = 0.6)

abline(h = median_reads, col = "black", lwd = 2, lty = 2)

## 50% of median threshold
half_median_reads <- 0.5 * median_reads
abline(h = half_median_reads, col = "black", lwd = 2, lty = 3)

legend(
    "topleft",
    legend = names(population),
    fill = population,
    border = NA,
    bty = "n",
    cex = 0.8
)

text(
    x = par("usr")[2],
    y = median_reads,
    labels = paste0("Median (", round(median_reads, 2), "M)"),
    pos = 2,
    cex = 0.8
)

text(
    x = par("usr")[2],
    y = half_median_reads,
    labels = paste0("50% Median (", round(half_median_reads, 2), "M)"),
    pos = 2,
    cex = 0.8
)

dev.off()
