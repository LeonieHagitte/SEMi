plotdat <- results %>% filter(popmodel == 1.21) %>%
  pivot_longer(
    cols = c(mnlfa_scalar_lrt_reject, tree_metric_split),
    names_to = "method",
    values_to = "rejected"
  ) %>%
  group_by(N, moderator, method) %>%
  summarise(
    rate = mean(rejected, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    method = recode(
      method,
      mnlfa_metric_lrt_reject = "MNLFA metric LRT",
      tree_metric_split = "Tree metric split"
    )
  )

ggplot(plotdat, aes(x = method, y = rate, fill = method)) +
  geom_hline(yintercept=0.05,lty=2)+
  geom_col(width = 0.7) +
  facet_grid(moderator ~ N) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = NULL,
    y = "MI Detection Rate",
    fill = "Method"
  ) +
  theme_minimal()
