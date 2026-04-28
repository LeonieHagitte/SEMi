# Packages
library(dplyr)
library(tidyr)
library(stringr)
library(gt)
library(scales)  # for percent_format()

# ---------- 1) Long per-method events ----------
# Add a single 'method' + 'analytical_model' channel
det_events <- det_tbl %>%
  pivot_longer(
    c(detect_mnlfa, detect_treeM),
    names_to  = "method_key",
    values_to = "detected"
  ) %>%
  mutate(
    method = recode(method_key,
                    detect_mnlfa = "MNLFA",
                    detect_treeM = "SEMTREE"),
    analytical_model = ifelse(method == "MNLFA", variant_label, "(SEMTREE)")
  ) %>%
  select(-method_key)

# For nice, compact column headers (Series = Method + Model)
det_events <- det_events %>%
  mutate(series = paste(method, analytical_model, sep = " : "))

# ---------- 2) TPR (non-NULL DGP only) ----------
tpr_wide <- det_events %>%
  filter(true_mod, model_type != "NULL") %>%
  group_by(model_type, N, moderator_1_type, series) %>%
  summarise(rate = mean(detected, na.rm = TRUE), .groups = "drop") %>%
  mutate(pct = 100 * rate) %>%
  select(-rate) %>%
  pivot_wider(
    names_from = series,
    values_from = pct
  ) %>%
  arrange(model_type, N, moderator_1_type)

# ---------- 3) FPR (NULL DGP only) ----------
fpr_wide <- det_events %>%
  filter(!true_mod, model_type == "NULL") %>%
  group_by(model_type, N, moderator_1_type, series) %>%
  summarise(rate = mean(detected, na.rm = TRUE), .groups = "drop") %>%
  mutate(pct = 100 * rate) %>%
  select(-rate) %>%
  pivot_wider(
    names_from = series,
    values_from = pct
  ) %>%
  arrange(model_type, N, moderator_1_type)

# ---------- 4) Render APA-ish GT tables (percentages) ----------
fmt_pct0 <- function(x) scales::percent(x/100, accuracy = 1)  # e.g., "67%"

# define nicer names for the analysis-model columns
pretty_names <- c(
  "MNLFA : MNLFA_linear_full"    = "MNLFA (linear full)",
  "MNLFA : MNLFA_linear_partial" = "MNLFA (linear partial)",
  "MNLFA : MNLFA_none"           = "MNLFA (null)",
  "SEMTREE : (SEMTREE)"          = "Sem Tree"
)

# TPR table
tpr_cols <- setdiff(names(tpr_wide), c("model_type","N","moderator_1_type"))
tpr_table <- tpr_wide %>%
  gt(rowname_col = NULL, groupname_col = "Model Type") |>
  fmt(
    columns = all_of(tpr_cols),
    fns = fmt_pct0
  ) |>
  fmt_number(columns = "N", decimals = 0) |>
  cols_label(
    N = html("<i>N</i>"),
    moderator_1_type = "Moderator form",
    !!!pretty_names 
  ) |>
  tab_header(title = md("**Detection Rates (Power; TPR)**")) |>
  tab_options(
    table.font.names = "Times New Roman",
    table.font.size = px(12),
    data_row.padding = px(4)
  )

# FPR table
fpr_cols <- setdiff(names(fpr_wide), c("model_type","N","moderator_1_type"))
fpr_table <- fpr_wide %>%
  gt(rowname_col = NULL, groupname_col = "Model Type") |>
  fmt(
    columns = all_of(fpr_cols),
    fns = fmt_pct0
  ) |>
  fmt_number(columns = "N", decimals = 0) |>
  cols_label(
    N = html("<i>N</i>"),
    moderator_1_type = "Moderator Form",
    !!!pretty_names 
  ) |>
  tab_header(title = md("**Detection Rates (Type I; FPR)**")) |>
  tab_options(
    table.font.names = "Times New Roman",
    table.font.size = px(12),
    data_row.padding = px(4)
  )

# View in RStudio:
tpr_table
fpr_table

# Optional: save to disk (HTML/RTF/PNG)
gt::gtsave(tpr_table, "results/apa_tpr_table.html")
gt::gtsave(fpr_table, "results/apa_fpr_table.html")
# gt::gtsave(tpr_table, "results/apa_tpr_table.rtf")
# gt::gtsave(fpr_table, "results/apa_fpr_table.rtf")

t_tpr <- tpr_wide %>% mutate(metric = "TPR")
t_fpr <- fpr_wide %>% mutate(metric = "FPR")
t_both <- bind_rows(t_tpr, t_fpr) %>%
  arrange(metric, model_type, N, moderator_1_type)

both_cols <- setdiff(names(t_both), c("metric","model_type","N","moderator_1_type"))
both_table <- t_both %>%
  gt(rowname_col = NULL, groupname_col = "metric") |>
  fmt(
    columns = all_of(both_cols),
    fns = fmt_pct0
  ) |>
  fmt_number(columns = "N", decimals = 0) |>
  cols_label(
    N = html("<i>N</i>"),
    moderator_1_type = "Moderator Form"
  ) |>
  tab_header(title = md("**Detection Rates (TPR and FPR)**")) |>
  tab_options(
    table.font.names = "Times New Roman",
    table.font.size = px(12),
    data_row.padding = px(4)
  )

both_table
gt::gtsave(both_table, "results/apa_detection_rates.html")

