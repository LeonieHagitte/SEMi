library(dplyr)
library(ggplot2)
library(scales)

# ------------------------------------------------------------

fig3_base <- results2_confirmatory %>%
  filter(
    true_metric_noninvariance == FALSE,
    delta_lambda == 0.2,
    delta_nu == 0.5,
    N == 1000,
    popmodel %in% c("1.11", "1.21")
  ) %>%
  mutate(
    model_group = factor(
      popmodel,
      levels = c("1.11", "1.21"),
      labels = c(
        "L0IA",
        "L0I2"
      )
    )
  )

# ------------------------------------------------------------

fig3_dat <- bind_rows(
  
  fig3_base %>%
    filter(analysis_form == "linear") %>%
    mutate(
      method = "MNLFA linear",
      detected = mnlfa_metric_lrt_reject_bonf
    ),
  
  fig3_base %>%
    filter(analysis_form == "quadratic") %>%
    mutate(
      method = "MNLFA quadratic",
      detected = mnlfa_metric_lrt_reject_bonf
    ),
  
  fig3_base %>%
    filter(analysis_form == "linear") %>%
    mutate(
      method = "SEM Tree",
      detected = tree_metric_reject_bonf
    )
  
) %>%
  filter(moderator != "noise") %>%
  group_by(
    reliability,
    moderator,
    method,
    model_group
  ) %>%
  summarise(
    rate = mean(detected, na.rm = TRUE),
    n = sum(!is.na(detected)),
    mcse = sqrt(rate * (1 - rate) / n),
    lo = pmax(0, rate - mcse),
    hi = pmin(1, rate + mcse),
    .groups = "drop"
  ) %>%
  mutate(
    method_num = case_when(
      method == "MNLFA linear" ~ 1.4,
      method == "MNLFA quadratic" ~ 2.0,
      method == "SEM Tree" ~ 2.6
    ),
    
    # Horizontal separation of the two population models
    x_position = method_num +
      if_else(
        model_group == "L0IA",
        -0.075,
        0.075
      )
  )

# ------------------------------------------------------------

method_shapes <- scale_shape_manual(
  name = "Method",
  values = c(
    "MNLFA linear" = 16,
    "MNLFA quadratic" = 15,
    "SEM Tree" = 17
  ),
  breaks = c("MNLFA linear", "MNLFA quadratic", "SEM Tree")
)

model_colours <- scale_colour_manual(
  name = "Population model",
  values = c(
    "L0IA" = "#343231",
    "L0I2" = "#A5A5A5"
  ),
  breaks = c("L0IA", "L0I2")
)

fig3 <- ggplot(
  fig3_dat,
  aes(
    x = x_position,
    y = rate,
    shape = method,
    colour = model_group
  )
) +
  
  # Nominal alpha
  geom_hline(
    yintercept = .05,
    linetype = "dashed",
    colour = "grey50",
    linewidth = .4
  ) +
  
  # Monte Carlo standard-error intervals
  geom_errorbar(
    aes(
      ymin = lo,
      ymax = hi
    ),
    width = .055,
    linewidth = .35
  ) +
  
  # Point estimates
  geom_point(
    size = 2.5
  ) +
  
  method_shapes +
  model_colours + 
  
  scale_x_continuous(
    limits = c(0.9, 3.1),
    breaks = NULL,
    expand = expansion(mult = 0)
  ) +
  
  facet_wrap(
    ~ moderator,
    nrow = 1,
    labeller = labeller(
      moderator = c(
        linear = "Linear",
        quadratic = "Quadratic",
        sigmoid = "Sigmoid"
      )
    )
  ) +
  
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    breaks = c(0, .20, .40, .60, .80)
  ) +
  
  coord_cartesian(
    ylim = c(0, .80)
  ) +
  
  labs(
    x = NULL,
    y = "Metric Rejection Rate"
  ) +
  
  theme_pub() +
  
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_blank(),
    
    axis.line = element_line(
      linewidth = .3,
      colour = "#343231"
    ),
    
    legend.position = "bottom",
    legend.background = element_blank(),
    
    plot.margin = margin(
      t = 15,
      r = 15,
      b = 15,
      l = 15
    )
  )

fig3

#----------------------------------------
fig0_base <- results2_confirmatory %>%
  filter(
    true_metric_noninvariance == FALSE,
    N == 1000,
    popmodel == "0"
  )

# Metric-test false positives
# no loading moderation, but intercept moderation present

fig0_dat <- bind_rows(
  
  fig0_base %>%
    filter(analysis_form == "linear") %>%
    mutate(
      method = "MNLFA linear",
      detected = mnlfa_metric_lrt_reject_bonf
    ),
  
  fig0_base %>%
    filter(analysis_form == "quadratic") %>%
    mutate(
      method = "MNLFA quadratic",
      detected = mnlfa_metric_lrt_reject_bonf
    ),
  
  fig0_base %>%
    filter(analysis_form == "linear") %>%   # one copy only
    mutate(
      method = "SEM Tree",
      detected = tree_metric_reject_bonf
    )
  
) %>%
  filter(moderator != "noise") %>%
  group_by(reliability, moderator, method) %>%
  summarise(
    rate = mean(detected, na.rm = TRUE),
    n = sum(!is.na(detected)),
    mcse = sqrt(rate * (1 - rate) / n),
    lo = pmax(0, rate - mcse),
    hi = pmin(1, rate + mcse),
    .groups = "drop"
  )

fig0_dat <- fig0_dat %>%
  mutate(
    method_num = case_when(
      method == "MNLFA linear" ~ 1.4,
      method == "MNLFA quadratic" ~ 2,
      method == "SEM Tree" ~ 2.6
    )
  )

fig0a <- ggplot(
  fig0_dat,
  aes(
    x = method_num,
    y = rate,
    shape = method
  )
) +
  geom_hline(
    yintercept = .05,
    linetype = "dashed",
    colour = "grey50",
    linewidth = .4
  ) +
  geom_errorbar(
    aes(
      ymin = lo,
      ymax = hi
    ),
    width = .12,
    linewidth = .35,
    colour = "#343231"
  ) +
  geom_point(
    size = 2.5,
    colour = "#343231"
  ) +
  method_shapes +
  scale_x_continuous(
    limits = c(0.9, 3.1),
    breaks = NULL,
    expand = expansion(mult = 0)
  ) +
  facet_wrap(
    ~ moderator,
    nrow = 1,
    labeller = labeller(
      moderator = c(
        linear = "Linear",
        quadratic = "Quadratic",
        sigmoid = "Sigmoid"
      )
    )
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, .80),
    breaks = c(0, .20, .40, .60, .80)
  ) +
  labs(
    x = NULL,
    y = "Metric Rejection Rate",
    shape = NULL
  ) +
  theme_pub() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_blank(),
    axis.line = element_line(linewidth = .3, colour = "#343231"),
    legend.position = "bottom",
    legend.background = element_blank(),
    plot.margin = margin(
      t = 15,
      r = 15,
      b = 15,
      l = 15
    )
  )

fig0a

#--------------------------------------------------------------------------
fig3b_base <- results2_confirmatory %>%
  filter(
    true_scalar_noninvariance == FALSE,
    N == 1000,
    popmodel == "1.1"
  )

fig3b_dat <- bind_rows(
  
  fig3b_base %>%
    filter(analysis_form == "linear") %>%
    mutate(
      method = "MNLFA linear",
      detected = mnlfa_scalar_lrt_reject_bonf
    ),
  
  fig3b_base %>%
    filter(analysis_form == "quadratic") %>%
    mutate(
      method = "MNLFA quadratic",
      detected = mnlfa_scalar_lrt_reject_bonf
    ),
  
  fig3b_base %>%
    filter(analysis_form == "linear") %>%   
    mutate(
      method = "SEM Tree",
      detected = tree_scalar_reject_bonf
    )
  
) %>%
  filter(moderator != "noise") %>%
  group_by(reliability, moderator, method) %>%
  summarise(
    rate = mean(detected, na.rm = TRUE),
    n = sum(!is.na(detected)),
    mcse = sqrt(rate * (1 - rate) / n),
    lo = pmax(0, rate - mcse),
    hi = pmin(1, rate + mcse),
    .groups = "drop"
  )

fig3b_dat <- fig3b_dat %>%
  mutate(
    method_num = case_when(
      method == "MNLFA linear" ~ 1.4,
      method == "MNLFA quadratic" ~ 2,
      method == "SEM Tree" ~ 2.6
    )
  )

fig3b <- ggplot(
  fig4b_dat,
  aes(
    x = method_num,
    y = rate,
    shape = method
  )
) +
  geom_hline(
    yintercept = .05,
    linetype = "dashed",
    colour = "grey50",
    linewidth = .4
  ) +
  geom_segment(
    aes(
      x = method_num - .06,
      xend = method_num + .06,
      y = lo,
      yend = lo
    ),
    linewidth = .35,
    colour = "#343231"
  ) +
  geom_segment(
    aes(
      x = method_num - .06,
      xend = method_num + .06,
      y = hi,
      yend = hi
    ),
    linewidth = .35,
    colour = "#343231"
  ) +
  geom_point(
    size = 2.5,
    colour = "#343231"
  ) +
  scale_shape_manual(values = shapes) +
  scale_x_continuous(
    limits = c(0.9, 3.1),
    breaks = NULL,
    expand = expansion(mult = 0)
  ) +
  facet_wrap(
    ~ moderator,
    nrow = 1,
    labeller = labeller(
      moderator = c(
        linear = "Linear",
        quadratic = "Quadratic",
        sigmoid = "Sigmoid"
      )
    )
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, .80),
    breaks = c(0, .20, .40, .60, .80)
  ) +
  labs(
    x = NULL,
    y = "Scalar Rejection Rate",
    shape = NULL
  ) +
  theme_pub() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_blank(),
    axis.line = element_line(linewidth = .3, colour = "#343231"),
    legend.position = "inside",
    legend.position.inside = c(.98, .98),
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    plot.margin = margin(
      t = 15,
      r = 15,
      b = 15,
      l = 15
    )
  )

fig3b

saveRDS(
  fig3,
  "results/fig3.rds"
)

saveRDS(
  fig0a,
  "results/fig0a.rds"
)

saveRDS(
  fig3b,
  "results/fig3b.rds"
)
