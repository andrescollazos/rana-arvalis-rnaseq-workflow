setwd(file.path(Sys.getenv("THESIS_DIR"), "code/04_de"))

# Load libraries
library(DESeq2)
library(dplyr)
library(tibble)
library(tidyr)
library(UpSetR)

# Load workspace
load("dds.RData")

# -----------------------------------------------------------------------------
# 1. Temperature effects: 20 vs 15 within each population
# -----------------------------------------------------------------------------

alpha <- 0.05

temp_contrasts <- list(
    "VF"    = c("temperature", "20", "15"),
    "C.Fin" = list(c("temperature_20_vs_15", "populationC.Fin.temperature20")),
    "E"     = list(c("temperature_20_vs_15", "populationE.temperature20")),
    "L"     = list(c("temperature_20_vs_15", "populationL.temperature20")),
    "Ka"    = list(c("temperature_20_vs_15", "populationKa.temperature20")),
    "Upp"   = list(c("temperature_20_vs_15", "populationUpp.temperature20")),
    "NA"    = list(c("temperature_20_vs_15", "populationNA.temperature20")),
    "NL"    = list(c("temperature_20_vs_15", "populationNL.temperature20"))
)

res_temp <- list()

for (pop in names(temp_contrasts)) {
    contrast_spec <- temp_contrasts[[pop]]

    res_i <- results(dds, contrast = contrast_spec, alpha = alpha)

    res_temp[[pop]] <- as.data.frame(res_i) %>%
        rownames_to_column("gene") %>%
        transmute(
            gene = gene,
            !!paste0(pop, "_log2FC") := log2FoldChange,
            !!paste0(pop, "_padj") := padj,
            !!paste0(pop, "_sig") := !is.na(padj) & padj < alpha
        )
}

temp_merged <- Reduce(function(x, y) full_join(x, y, by = "gene"), res_temp)

# -----------------------------------------------------------------------------
# 2. UpSet plot input
#    This shows overlap of significant temperature-responsive genes across
#    populations. Each population is one set.
# -----------------------------------------------------------------------------

upset_df <- temp_merged %>%
    select(gene, ends_with("_sig")) %>%
    rename_with(~ sub("_sig$", "", .x), ends_with("_sig")) %>%
    filter(rowSums(across(-gene)) > 0)

upset_mat <- upset_df %>%
    select(-gene) %>%
    mutate(across(everything(), as.integer))

pdf(file.path(Sys.getenv("THESIS_DIR"), "analyses/04_de/01_temperature_effects/01_upset.pdf"),
    width = 10, height = 7
)

upset(
    upset_mat,
    sets = c("Ka", "Upp", "L", "E", "NA", "NL", "VF", "C.Fin"),
    order.by = "freq",
    keep.order = TRUE,
    main.bar.color = "black",
    sets.bar.color = "grey40",
    text.scale = 1.2
)
dev.off()

# -----------------------------------------------------------------------------
# 3. Per-gene direction-consistency classification
#    Rules used here:
#    - no_sig: not significant in any population
#    - population_specific: significant in exactly 1 population
#    - consistent_up: significant in >=2 populations and all significant log2FC > 0
#    - consistent_down: significant in >=2 populations and all significant log2FC < 0
#    - mixed_direction: significant in >=2 populations with both positive and negative log2FC
# -----------------------------------------------------------------------------

pop_names <- names(temp_contrasts)

sig_cols <- paste0(pop_names, "_sig")
lfc_cols <- paste0(pop_names, "_log2FC")

direction_per_gene <- temp_merged %>%
    rowwise() %>%
    mutate(
        # Count how many populations show a significant temperature effect
        n_sig = sum(c_across(all_of(sig_cols)), na.rm = TRUE),
        # Count number of populations with significant positive effect
        n_pos = {
            sig_vals <- c_across(all_of(sig_cols))
            lfc_vals <- c_across(all_of(lfc_cols))

            sum(sig_vals & !is.na(lfc_vals) & lfc_vals > 0, na.rm = TRUE)
        },
        # Count number of populations with significant negative effect
        n_neg = {
            sig_vals <- c_across(all_of(sig_cols))
            lfc_vals <- c_across(all_of(lfc_cols))

            sum(sig_vals & !is.na(lfc_vals) & lfc_vals < 0, na.rm = TRUE)
        },
        # Classify the gene according to the pattern of responses
        direction_class = case_when(
            n_sig == 0 ~ "no_sig",
            n_sig <= 2 ~ "population_specific",
            n_sig >= 5 & n_pos > 0 & n_neg == 0 ~ "shared_up",
            n_sig >= 5 & n_neg > 0 & n_pos == 0 ~ "shared_down",
            n_pos >= 2 & n_neg > 0 ~ "mixed_direction",
            TRUE ~ "unclassified"
        )
    ) %>%
    # remove rowwise grouping
    ungroup()

# -----------------------------------------------------------------------------
# 4. Summary table
#    Counts genes in each direction-consistency category
# -----------------------------------------------------------------------------

direction_summary <- direction_per_gene %>%
    count(direction_class, sort = TRUE)

View(direction_summary)

save.image("temperature_effects.RData")
