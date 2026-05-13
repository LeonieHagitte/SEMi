
# TODO: load all results and merge

results %>%
  dplyr::count(popmodel, moderator, true_any_noninvariance)

results %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_errors = sum(!is.na(error_msg)),
    mnlfa_metric_reject_rate = mean(mnlfa_metric_lrt_reject, na.rm = TRUE),
    mnlfa_scalar_reject_rate = mean(mnlfa_scalar_lrt_reject, na.rm = TRUE),
    tree_metric_reject_rate = mean(tree_metric_reject, na.rm = TRUE),
    tree_scalar_reject_rate = mean(tree_scalar_reject, na.rm = TRUE)
  )
