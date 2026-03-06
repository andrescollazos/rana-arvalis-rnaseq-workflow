# PCA_regression_pop.R

run_pca_regression_pop <- function(pca_obj, n_pcs = 5) {
    # Extract PCA dataframe (contains PC scores + metadata)
    pca_df <- pca_obj$pca_df

    # Ensure factors
    pca_df$temperature <- factor(pca_df$temperature)

    # Storage
    results_list <- list()
    pval_matrix <- matrix(NA, nrow = n_pcs, ncol = 1)
    colnames(pval_matrix) <- c("temperature")

    for (i in seq_len(n_pcs)) {
        pc_name <- paste0("PC", i)

        formula_str <- paste0(
            pc_name,
            " ~ temperature"
        )

        fit <- lm(as.formula(formula_str), data = pca_df)

        anova_tab <- anova(fit)

        # Extract p-values for each term
        pvals <- anova_tab$`Pr(>F)`[1]
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
        temperature = pval_adj_matrix[, "temperature"],
        row.names = NULL,
        check.names = FALSE
    )

    return(final_table)
}
