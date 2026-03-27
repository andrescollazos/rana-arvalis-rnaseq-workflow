mean_abs_lfc <- apply(abs(lfc_shared), 1, function(x) mean(x, na.rm = TRUE))

df_shared_test <- data.frame(
    mean_abs_lfc = mean_abs_lfc,
    concordant = shared_concordant
)

boxplot(mean_abs_lfc ~ concordant, data = df_shared_test)

wilcox.test(mean_abs_lfc ~ concordant, data = df_shared_test)

df_shared_test$bin <- cut(
    df_shared_test$mean_abs_lfc,
    breaks = quantile(df_shared_test$mean_abs_lfc, probs = seq(0, 1, 0.2), na.rm = TRUE),
    include.lowest = TRUE
)

prop_true <- tapply(df_shared_test$concordant, df_shared_test$bin, mean)

prop_true

plot(
    df_shared_test$mean_abs_lfc,
    as.numeric(df_shared_test$concordant),
    pch = 16,
    col = rgb(0, 0, 0, 0.2)
)
