library(ggplot2)

plotdat <- results %>% filter(popmodel == 1.21) %>%
  pivot_longer(
    cols = c(mnlfa_scalar_lrt_reject, tree_metric_split),
    names_to = "method",
    values_to = "rejected"
  ) %>%
  group_by(N, moderator, delta_lambda, method) %>%
  summarise(
    rate = mean(rejected, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    method = recode(
      method,
      mnlfa_scalar_lrt_reject = "MNLFA scalar LRT",
      tree_metric_split = "Tree metric split"
    )
  )

ggplot(plotdat, aes(x = delta_lambda, y = rate, color = method, group = method)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  facet_grid(moderator ~ N) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = expression(delta_lambda),
    y = "MI Detection Rate",
    color = "Method"
  ) +
  theme_minimal()