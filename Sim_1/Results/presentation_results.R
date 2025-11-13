# Ensure we have a readable label for the MNLFA variant
results$variant_label <- dplyr::coalesce(results$variant_label, results$assumed_mnlfa_model)

# Stable ordering for facets
results$moderator_1_type <- factor(results$moderator_1_type,
                                   levels = c("linear","sigmoid","quadratic","noise"))

# ---------- 1) Raw detection outcomes ----------
det_tbl <- results %>%
  dplyr::transmute(
    model_type, N, moderator_1_type, variant_label,
    true_mod      = model_type != "NULL",
    detect_mnlfa  = dplyr::coalesce(as.logical(mnlfa_detect), FALSE),
    detect_treeM  = dplyr::coalesce(as.logical(tree_split_on_M), FALSE),
    detect_treeM2 = dplyr::coalesce(as.logical(tree_split_on_M2), FALSE)
  )

# ---------- 2) Method-long events ----------
det_events <- det_tbl %>%
  tidyr::pivot_longer(
    cols = c(detect_mnlfa, detect_treeM),
    names_to  = "method_key",
    values_to = "detected"
  ) %>%
  dplyr::mutate(
    method = dplyr::recode(method_key,
                           detect_mnlfa = "MNLFA",
                           detect_treeM = "SEMTREE"),
    analytical_model = ifelse(method == "MNLFA", variant_label, "(SEMTREE)")
  ) %>%
  dplyr::select(-method_key)

# ---------- 3) Compute detection rates per method × model ----------
# MNLFA
mnlfa_det <- det_events %>%
  dplyr::filter(method == "MNLFA") %>%
  dplyr::group_by(model_type, N, moderator_1_type, variant_label) %>%
  dplyr::summarise(
    TPR = mean(detected[true_mod],  na.rm = TRUE),
    FPR = mean(detected[!true_mod], na.rm = TRUE),
    .groups = "drop"
  )

# SEMTREE
tree_det <- det_events %>%
  dplyr::filter(method == "SEMTREE") %>%
  # TPR over non-NULL; FPR from NULL slice
  dplyr::group_by(model_type, N, moderator_1_type) %>%
  dplyr::summarise(
    TPR = mean(detected[true_mod],  na.rm = TRUE),
    FPR = mean(detected[!true_mod], na.rm = TRUE),
    .groups = "drop"
  )

# ---------- 4) Long table for plotting both methods ----------
det_by_model_long <- dplyr::bind_rows(
  mnlfa_det %>%
    dplyr::mutate(method = "MNLFA", analytical_model = variant_label) %>%
    tidyr::pivot_longer(c(TPR, FPR), names_to = "metric", values_to = "rate"),
  tree_det %>%
    dplyr::mutate(method = "SEMTREE", analytical_model = "(SEMTREE)") %>%
    tidyr::pivot_longer(c(TPR, FPR), names_to = "metric", values_to = "rate")
)

# Now that det_by_model_long exists, set a stable order for the legend:
det_by_model_long$analytical_model <- factor(
  det_by_model_long$analytical_model,
  levels = c("MNLFA_linear_full","MNLFA_linear_partial","MNLFA_none","(SEMTREE)")
)

# ---------- 5) Counts per cell (for tables) ----------
n_mnlfa <- det_tbl %>%
  dplyr::group_by(model_type, N, moderator_1_type, variant_label) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop")

n_tree <- det_tbl %>%
  dplyr::group_by(model_type, N, moderator_1_type) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop")

# ---------- 6) Build gt tables ----------
library(gt)

mnlfa_table <- mnlfa_det %>%
  dplyr::left_join(n_mnlfa, by = c("model_type","N","moderator_1_type","variant_label")) %>%
  dplyr::mutate(
    TPR = ifelse(is.finite(TPR), TPR, NA_real_),
    FPR = ifelse(is.finite(FPR), FPR, NA_real_)
  ) %>%
  dplyr::arrange(model_type, N, moderator_1_type, variant_label) %>%
  dplyr::mutate(model_type = factor(model_type, levels = c("NULL","1.1","1.2"))) %>%
  gt(rowname_col = NULL, groupname_col = "model_type") |>
  fmt_number(columns = c(N, TPR, FPR), decimals = 2) |>
  cols_label(
    N = html("<i>N</i>"),
    moderator_1_type = "Moderator link",
    variant_label = "MNLFA variant",
    TPR = "Power (TPR)",
    FPR = "Type I (FPR)",
    n = "Reps"
  ) |>
  tab_options(
    table.font.names = "Times New Roman",
    table.font.size = px(12),
    data_row.padding = px(4)
  ) |>
  tab_header(title = md("**Detection Rates for MNLFA**"))

tree_table <- tree_det %>%
  dplyr::left_join(n_tree, by = c("model_type","N","moderator_1_type")) %>%
  dplyr::mutate(
    TPR = ifelse(is.finite(TPR), TPR, NA_real_),
    FPR = ifelse(is.finite(FPR), FPR, NA_real_)
  ) %>%
  dplyr::arrange(model_type, N, moderator_1_type) %>%
  dplyr::mutate(model_type = factor(model_type, levels = c("NULL","1.1","1.2"))) %>%
  gt(rowname_col = NULL, groupname_col = "model_type") |>
  fmt_number(columns = c(N, TPR, FPR), decimals = 2) |>
  cols_label(
    N = html("<i>N</i>"),
    moderator_1_type = "Moderator link",
    TPR = "Power (TPR)",
    FPR = "Type I (FPR)",
    n = "Reps"
  ) |>
  tab_options(
    table.font.names = "Times New Roman",
    table.font.size = px(12),
    data_row.padding = px(4)
  ) |>
  tab_header(title = md("**Detection Rates for SEMTREE**"))

# ---------- 7) Plots ----------
library(ggplot2)
library(rlang)

.keep01 <- function(df, col = "rate") {
  col_sym <- sym(col)
  df %>%
    dplyr::filter(is.finite(!!col_sym)) %>%
    dplyr::mutate(!!col_sym := pmin(pmax(!!col_sym, 0), 1))
}


plot_tpr <- det_by_model_long %>%
  dplyr::filter(metric == "TPR", model_type != "NULL") %>%
  dplyr::mutate(series = interaction(method, analytical_model, sep = " : ")) %>%
  .keep01("rate") %>%
  ggplot(aes(x = factor(N), y = rate,
             color = series, shape = series, linetype = series, group = series)) +
  geom_point(size = 2) +
  geom_line() +
  facet_grid(moderator_1_type ~ model_type) +
  labs(x = "Sample size N", y = "Detection rate (TPR)",
       color = "Analysis model", shape = "Analysis model", linetype = "Analysis model") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw()

p_mnlfa_tpr <- mnlfa_det %>%
  dplyr::filter(model_type != "NULL") %>%
  ggplot(aes(x = factor(N), y = TPR, color = variant_label, group = variant_label)) +
  geom_point(size = 2) + geom_line() +
  facet_grid(moderator_1_type ~ model_type) +
  labs(x = "Sample size N", y = "Detection rate (TPR)", color = "MNLFA variant") +
  coord_cartesian(ylim = c(0,1)) + theme_bw()

p_mnlfa_fpr <- mnlfa_det %>%
  dplyr::filter(model_type == "NULL") %>%
  .keep01("FPR") %>%
  ggplot(aes(x = factor(N), y = FPR, fill = variant_label)) +
  geom_col(position = position_dodge(0.6), width = 0.6, color = "grey30") +
  facet_wrap(~ moderator_1_type) +
  labs(x = "Sample size N", y = "False-positive rate (FPR)", fill = "MNLFA variant") +
  coord_cartesian(ylim = c(0,1)) + theme_bw()

p_tree_tpr <- tree_det %>%
  dplyr::filter(model_type != "NULL") %>%
  ggplot(aes(x = factor(N), y = TPR, group = 1)) +
  geom_point(size = 2) + geom_line() +
  facet_grid(moderator_1_type ~ model_type) +
  labs(x = "Sample size N", y = "Detection rate (TPR)", title = "SEMTREE") +
  coord_cartesian(ylim = c(0,1)) + theme_bw()

p_tree_fpr <- tree_det %>%
  dplyr::filter(model_type == "NULL") %>%
  .keep01("FPR") %>%
  ggplot(aes(x = factor(N), y = FPR, fill = moderator_1_type)) +
  geom_col(position = position_dodge(0.6), width = 0.6, color = "grey30") +
  labs(x = "Sample size N", y = "False-positive rate (FPR)", fill = "Moderator") +
  coord_cartesian(ylim = c(0,1)) + theme_bw()

plot_fpr <- det_by_model_long %>%
  dplyr::filter(metric == "FPR", model_type == "NULL") %>%
  dplyr::mutate(series = interaction(method, analytical_model, sep = " : ")) %>%
  .keep01("rate") %>%
  ggplot(aes(x = factor(N), y = rate,
             color = series, shape = series, linetype = series, group = series)) +
  geom_point(size = 2) +
  geom_line() +
  facet_wrap(~ moderator_1_type) +
  labs(x = "Sample size N", y = "False-positive rate (FPR)",
       color = "Analysis model", shape = "Analysis model", linetype = "Analysis model") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw()



#-----------------------------



library(dplyr)
library(tidyr)
library(ggplot2)
library(binom)

# det_events must already exist (from your earlier code) and contain:
# model_type, N, moderator_1_type, method, analytical_model, detected (0/1), true_mod (logical)

# ---------- TPR (non-NULL) with Wilson CI ----------
tpr_ci <- det_events %>%
  filter(true_mod, model_type != "NULL") %>%
  group_by(model_type, N, moderator_1_type, method, analytical_model) %>%
  summarise(k = sum(detected, na.rm = TRUE),
            n = n(), .groups = "drop") %>%
  mutate(rate = ifelse(n > 0, k / n, NA_real_)) %>%
  rowwise() %>%
  mutate(
    # keep only CI bounds to avoid name clashes
    ci_df = list(binom::binom.confint(k, n, methods = "wilson")[, c("lower","upper")])
  ) %>%
  tidyr::unnest_wider(ci_df) %>%
  ungroup() %>%
  mutate(series = interaction(method, analytical_model, sep = " : "))

# ---------- FPR (NULL only) with Wilson CI ----------
fpr_ci <- det_events %>%
  filter(!true_mod, model_type == "NULL") %>%
  group_by(model_type, N, moderator_1_type, method, analytical_model) %>%
  summarise(k = sum(detected, na.rm = TRUE),
            n = n(), .groups = "drop") %>%
  mutate(rate = ifelse(n > 0, k / n, NA_real_)) %>%
  rowwise() %>%
  mutate(
    ci_df = list(binom::binom.confint(k, n, methods = "wilson")[, c("lower","upper")])
  ) %>%
  tidyr::unnest_wider(ci_df) %>%
  ungroup() %>%
  mutate(series = interaction(method, analytical_model, sep = " : "))


# ---- TPR plot (non-NULL) ----
plot_tpr <- tpr_ci %>%
  ggplot(aes(x = factor(N), y = rate,
             color = series, shape = series, linetype = series, group = series)) +
  geom_point(size = 2) +
  geom_line() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 0.3) +
  facet_grid(moderator_1_type ~ model_type) +
  labs(x = "Sample size N", y = "Detection rate (TPR)",
       color = "Analysis model", shape = "Analysis model", linetype = "Analysis model") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw()

# ---- FPR plot (NULL) ----
plot_fpr <- fpr_ci %>%
  ggplot(aes(x = factor(N), y = rate,
             color = series, shape = series, linetype = series, group = series)) +
  geom_point(size = 2) +
  geom_line() +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.1, linewidth = 0.3) +
  facet_wrap(~ moderator_1_type) +
  labs(x = "Sample size N", y = "False-positive rate (FPR)",
       color = "Analysis model", shape = "Analysis model", linetype = "Analysis model") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw()
# Print
plot_tpr
plot_fpr


# ---------- 8) Save AFTER objects exist ----------
gt::gtsave(mnlfa_table, "results/mnlfa_detection_table.html")
gt::gtsave(tree_table,   "results/semtree_detection_table.html")

ggplot2::ggsave("results/plot_tpr.png", plot_tpr, width = 9, height = 6, dpi = 300)
ggplot2::ggsave("results/p_mnlfa_fpr.png", p_mnlfa_fpr, width = 8, height = 5, dpi = 300)
ggplot2::ggsave("results/p_tree_tpr.png",  p_tree_tpr,  width = 8, height = 5, dpi = 300)
ggplot2::ggsave("results/p_tree_fpr.png",  p_tree_fpr,  width = 8, height = 5, dpi = 300)
ggplot2::ggsave("results/plot_fpr.png", plot_fpr, width = 8, height = 5, dpi = 300)


# Combined plot of TPR across methods
print(plot_tpr)
print(plot_fpr)

# MNLFA only
print(p_mnlfa_tpr)
print(p_mnlfa_fpr)

# SEMTREE only
print(p_tree_tpr)
print(p_tree_fpr)



