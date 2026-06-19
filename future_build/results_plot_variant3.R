#'
#' Select a population model, then plot MLNFA(linear/quadratic)+SEM Tree
#' results by sample size and effect size
#'
#'

plotdat <- results %>%
  filter(popmodel == 1.22) %>%
  pivot_longer(
    cols = c(mnlfa_metric_lrt_reject, tree_metric_split),
    names_to = "method",
    values_to = "rejected"
  ) %>%
  mutate(
    method = case_when(
      method == "mnlfa_metric_lrt_reject" & analysis_form == "linear" ~
        "MNLFA metric LRT: linear",
      method == "mnlfa_metric_lrt_reject" & analysis_form == "quadratic" ~
        "MNLFA metric LRT: quadratic",
      method == "tree_metric_split" ~
        "Tree metric split",
      TRUE ~ method
    )
  ) %>%
  group_by(N, moderator, delta_lambda, method) %>%
  summarise(
    rate = mean(rejected, na.rm = TRUE),
    .groups = "drop"
  )


ggplot(plotdat, aes(x = delta_lambda, y = rate, color = method, group = method)) +
  geom_hline(yintercept=0.05,lty=2)+
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  facet_grid(moderator ~ N) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "Effect Size",
    y = "MI Detection Rate",
    color = "Method"
  ) + 
  theme_minimal()+
  theme(
    panel.spacing.y = unit(1.2, "lines"),
    panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.5)
  )
