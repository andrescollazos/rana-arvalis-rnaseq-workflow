# PCA_regression.R

run_pca_regression <- function(pca_obj, n_pcs = 5) {
    # Extract PCA dataframe (contains PC scores + metadata)
    pca_df <- pca_obj$pca_df

    # Ensure factors
    pca_df$extraction_date <- factor(pca_df$extraction_date)
    pca_df$family <- factor(pca_df$family)
    pca_df$stage <- factor(pca_df$stage)
    pca_df$population <- factor(pca_df$population)
    pca_df$temperature <- factor(pca_df$temperature)

    # Storage
    results_list <- list()
    pval_matrix <- matrix(NA, nrow = n_pcs, ncol = 6)

    colnames(pval_matrix) <- c(
        "extraction_date",
        "family",
        "stage",
        "population",
        "temperature",
        "population:temperature"
    )

    for (i in seq_len(n_pcs)) {
        pc_name <- paste0("PC", i)

        formula_str <- paste0(
            pc_name,
            " ~ extraction_date + family + stage + population * temperature"
        )

        fit <- lm(as.formula(formula_str), data = pca_df)

        anova_tab <- anova(fit)

        # Extract p-values for each term
        pvals <- anova_tab$`Pr(>F)`[1:6]
        pval_matrix[i, ] <- pvals

        results_list[[i]] <- list(
            PC = pc_name,
            R2 = summary(fit)$r.squared,
            pvals = pvals
        )
    }

    # Adjust p-values across all tests (BH)
    pval_adj <- p.adjust(as.vector(pval_matrix), method = "BH")
    pval_adj_matrix <- matrix(pval_adj, nrow = n_pcs, byrow = FALSE)
    colnames(pval_adj_matrix) <- colnames(pval_matrix)

    # Build final table
    final_table <- data.frame(
        PC = paste0("PC", seq_len(n_pcs)),
        R2 = sapply(results_list, function(x) x$R2),
        extraction_date = pval_adj_matrix[, "extraction_date"],
        family = pval_adj_matrix[, "family"],
        stage = pval_adj_matrix[, "stage"],
        population = pval_adj_matrix[, "population"],
        temperature = pval_adj_matrix[, "temperature"],
        interaction = pval_adj_matrix[, "population:temperature"],
        row.names = NULL,
        check.names = FALSE
    )

    return(final_table)
}
