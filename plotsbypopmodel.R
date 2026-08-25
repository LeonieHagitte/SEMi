make_metric_plot <- function(pop) {
  
  metric_perf <- results2 %>%
    filter(
      popmodel %in% pop,
      moderator != "noise"
    ) %>%
    mutate(
      method = case_when(
        analysis_form == "linear" ~ "MNLFA linear",
        analysis_form == "quadratic" ~ "MNLFA quadratic",
        TRUE ~ NA_character_
      ),
      detected = mnlfa_metric_lrt_reject_bonf,
      estimand = case_when(
        true_metric_noninvariance == TRUE  ~ "Metric power",
        true_metric_noninvariance == FALSE ~ "Metric rejection under loading invariance",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(method), !is.na(estimand)) %>%
    group_by(popmodel, estimand, method, moderator, N) %>%
    summarise(
      rate = mean(detected, na.rm = TRUE),
      n = sum(!is.na(detected)),
      .groups = "drop"
    ) %>%
    mutate(
      method = factor(method, levels = c("MNLFA linear", "MNLFA quadratic")),
      moderator = factor(
        moderator,
        levels = c("linear", "quadratic", "sigmoid"),
        labels = c("Linear", "Quadratic", "Sigmoid")
      ),
      N = factor(N)
    )
  
  y_lab <- unique(metric_perf$estimand)
  y_lab <- if (length(y_lab) == 1) y_lab else "Metric rejection rate"
  
  hline_dat <- metric_perf %>%
    filter(estimand == "Metric rejection under loading invariance") %>%
    distinct(popmodel, moderator)
  
  ggplot(
    metric_perf,
    aes(
      x = N,
      y = rate,
      colour = method,
      shape = method,
      group = method
    )
  ) +
    geom_hline(
      data = hline_dat,
      aes(yintercept = .05),
      inherit.aes = FALSE,
      linetype = "dashed",
      colour = "grey60"
    ) +
    geom_line(linewidth = .65, alpha = .75) +
    geom_point(size = 2.2) +
    facet_grid(popmodel ~ moderator, drop = TRUE) +
    scale_colour_manual(values = cols) +
    scale_shape_manual(values = shapes) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1),
      breaks = seq(0, 1, .25)
    ) +
    labs(
      title = paste("Population model(s):", paste(pop, collapse = ", ")),
      subtitle = y_lab,
      x = "Sample size",
      y = y_lab,
      colour = NULL,
      shape = NULL
    ) +
    theme_pub() +
    theme(
      axis.text.x = element_text(angle = 0)
    )
}

models <- c(
  "0",
  "1.1",
  "1.11",
  "1.12",
  "1.2",
  "1.21",
  "1.22",
  "1.3",
  "1.32"
)

null_models <- c("0")
scalar_only_models <- c("1.11", "1.21")
metric_models <- c("1.1", "1.12", "1.2", "1.22", "1.3", "1.32")

make_metric_plot(models[1])  # 0
make_metric_plot(models[2])  # 1.1
make_metric_plot(models[3])  # 1.11
make_metric_plot(models[4])  # 1.12
make_metric_plot(models[5])  # 1.2
make_metric_plot(models[6])  # 1.21
make_metric_plot(models[7])  # 1.22
make_metric_plot(models[8])  # 1.3
make_metric_plot(models[9])  # 1.32

make_metric_plot(c("1.11", "1.21"))

make_metric_plot(c("1.1", "1.2"))

make_metric_plot(c("1.12", "1.22"))

make_metric_plot(c("1.3", "1.32"))
