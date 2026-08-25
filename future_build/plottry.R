library(tidyverse)
library(scales)
library(viridisLite)

# -------------------

to_logical <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(as.logical(x))
  x <- as.character(x)
  case_when(
    x %in% c("TRUE", "True", "true", "1") ~ TRUE,
    x %in% c("FALSE", "False", "false", "0") ~ FALSE,
    TRUE ~ NA
  )
}

alpha <- 0.05

results2 <- results %>%
  mutate(
    popmodel = as.character(popmodel),
    N = as.integer(N),
    reliability = as.numeric(reliability),
    delta_lambda = as.numeric(delta_lambda),
    delta_nu = as.numeric(delta_nu),
    moderator = as.character(moderator),
    analysis_form = as.character(analysis_form),
    
    across(
      c(
        true_any_noninvariance,
        true_metric_noninvariance,
        true_scalar_noninvariance,
        true_structured_moderator,
        mnlfa_metric_lrt_reject,
        mnlfa_scalar_lrt_reject,
        mnlfa_omnibus_lrt_reject,
        tree_metric_split,
        tree_scalar_split,
        tree_metric_reject,
        tree_scalar_reject,
        tree_metric_split_on_am1,
        tree_metric_split_on_am2,
        tree_metric_split_on_m0,
        tree_scalar_split_on_am1,
        tree_scalar_split_on_am2,
        tree_scalar_split_on_m0
      ),
      to_logical
    )
  ) %>%
  mutate(
    mnlfa_metric_lrt_reject_bonf = case_when(
      is.na(mnlfa_metric_lrt_p) ~ NA,
      mnlfa_metric_lrt_p <= alpha / 2 ~ TRUE,
      TRUE ~ FALSE
    ),
    mnlfa_scalar_lrt_reject_bonf = case_when(
      is.na(mnlfa_scalar_lrt_p) ~ NA,
      mnlfa_scalar_lrt_p <= alpha / 2 ~ TRUE,
      TRUE ~ FALSE
    ),
    tree_metric_reject_bonf = case_when(
      is.na(tree_metric_p) ~ NA,
      tree_metric_p <= alpha ~ TRUE,
      TRUE ~ FALSE
    ),
    tree_scalar_reject_bonf = case_when(
      is.na(tree_scalar_p) ~ NA,
      tree_scalar_p <= alpha ~ TRUE,
      TRUE ~ FALSE
    ),
    
    tree_metric_correct_split = case_when(
      popmodel %in% c("1.1", "1.12", "1.2", "1.22", "1.32") ~ tree_metric_split_on_am1,
      popmodel == "1.3" ~ tree_metric_split_on_am1 | tree_metric_split_on_am2,
      TRUE ~ FALSE
    ),
    
    tree_scalar_correct_split = case_when(
      popmodel %in% c("1.11", "1.12", "1.21", "1.22") ~ tree_scalar_split_on_am1,
      popmodel == "1.32" ~ tree_scalar_split_on_am2,
      TRUE ~ FALSE
    ),
    
    semtree_metric_power_decision_bonf =
      tree_metric_reject_bonf == TRUE & tree_metric_correct_split == TRUE,
    
    semtree_scalar_power_decision_bonf =
      tree_scalar_reject_bonf == TRUE & tree_scalar_correct_split == TRUE
  )
# --------------------
theme_pub <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", linewidth = .3),
      panel.spacing = unit(.9, "lines"),
      panel.border = element_rect(color = "grey80", fill = NA, linewidth = .35),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 10),
      axis.text = element_text(color = "grey25"),
      axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
      axis.title = element_text(face = "bold", color = "grey15")
    )
}
#   VanGogh1 = list(c("#2c2d54", "#434475", "#6b6ca3", "#969bc7", "#87bcbd", "#89ab7c", "#6f9954"), c(3, 5, 7, 4, 6, 2, 1), colorblind=FALSE),
cols <- c(
  "MNLFA linear" = "#2c2d54",
  "MNLFA quadratic" = "#6b6ca3",
  "SEM Tree" = "#6f9954"
)

shapes <- c(
  "MNLFA linear" = 16,
  "MNLFA quadratic" = 15,
  "SEM Tree" = 17
)
# ----------------------------------------------------------------------

# -------------------------
fig1_dat <- bind_rows(
  results2 %>%
    filter(true_metric_noninvariance == TRUE) %>%
    mutate(
      method = case_when(
        analysis_form == "linear" ~ "MNLFA linear",
        analysis_form == "quadratic" ~ "MNLFA quadratic"
      ),
      detected = mnlfa_metric_lrt_reject_bonf
    ) %>%
    filter(!is.na(method)),
  
  results2 %>%
    filter(true_metric_noninvariance == TRUE) %>%
    mutate(
      method = "SEM Tree",
      detected = semtree_metric_power_decision_bonf
    )
) %>%
  filter(popmodel == "1.22") %>%
  group_by(N, moderator, delta_lambda, method) %>%
  summarise(
    rate = mean(detected, na.rm = TRUE),
    n = sum(!is.na(detected)),
    se = sqrt(rate * (1 - rate) / n),
    lo = pmax(0, rate - 1.96 * se),
    hi = pmin(1, rate + 1.96 * se),
    low_n = n < 30,
    .groups = "drop"
  )

fig1 <- ggplot(
  fig1_dat,
  aes(
    delta_lambda,
    rate,
    colour = method,
    fill = method,
    shape = method,
    group = method
  )
) +
  geom_hline(yintercept = .05, linetype = "dashed", color = "grey70") +
  geom_ribbon(
    data = fig1_dat %>% filter(!low_n),
    aes(ymin = lo, ymax = hi),
    alpha = .10,
    colour = NA
  ) +
  geom_line(linewidth = .70, alpha = .85) +
  geom_point(size = 2) +
  geom_point(
    data = fig1_dat %>% filter(low_n),
    aes(x = delta_lambda, y = rate),
    inherit.aes = FALSE,
    shape = 4,
    size = 3,
    stroke = .8,
    colour = "grey35"
  ) +
  facet_grid(moderator ~ N) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_shape_manual(values = shapes) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(
    x = "Loading effect size",
    y = "Statistical power",
    colour = NULL,
    fill = NULL,
    shape = NULL
  ) +
  theme_pub()

fig1
# ----------------------------------------------
fig2_dat <- bind_rows(
  results2 %>%
    filter(true_scalar_noninvariance == TRUE) %>%
    mutate(
      method = case_when(
        analysis_form == "linear" ~ "MNLFA linear",
        analysis_form == "quadratic" ~ "MNLFA quadratic"
      ),
      detected = mnlfa_scalar_lrt_reject_bonf
    ) %>%
    filter(!is.na(method)),
  
  results2 %>%
    filter(true_scalar_noninvariance == TRUE) %>%
    mutate(
      method = "SEM Tree",
      detected = semtree_scalar_power_decision_bonf
    )
) %>%
  filter(popmodel == "1.22") %>%
  group_by(N, moderator, delta_nu, method) %>%
  summarise(
    rate = mean(detected, na.rm = TRUE),
    n = sum(!is.na(detected)),
    se = sqrt(rate * (1 - rate) / n),
    lo = pmax(0, rate - 1.96 * se),
    hi = pmin(1, rate + 1.96 * se),
    low_n = n < 30,
    .groups = "drop"
  )


fig2 <- ggplot(
  fig2_dat,
  aes(
    delta_nu,
    rate,
    colour = method,
    fill = method,
    shape = method,
    group = method
  )
) +
  geom_hline(
    yintercept = .05,
    linetype = "dashed",
    color = "grey70"
  ) +
  geom_ribbon(
    data = fig2_dat %>% filter(!low_n),
    aes(ymin = lo, ymax = hi),
    alpha = .10,
    colour = NA
  ) +
  geom_line(linewidth = .70, alpha = .85) +
  geom_point(size = 2) +
  geom_point(
    data = fig2_dat %>% filter(low_n),
    aes(x = delta_nu, y = rate),
    inherit.aes = FALSE,
    shape = 4,
    size = 3,
    stroke = .8,
    colour = "grey35"
  ) +
  facet_grid(moderator ~ N) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_shape_manual(values = shapes) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(
    x = "Intercept effect size",
    y = "Statistical power",
    colour = NULL,
    fill = NULL,
    shape = NULL
  ) +
  theme_pub()

fig2

# -------------------------------------------------------------------
# -------------------------------------------------------------------
metric_perf <- bind_rows(
  
  results2 %>%
    filter(true_metric_noninvariance == TRUE) %>%
    mutate(
      estimand = "Power",
      method = case_when(
        analysis_form == "linear" ~ "MNLFA linear",
        analysis_form == "quadratic" ~ "MNLFA quadratic"
      ),
      detected = mnlfa_metric_lrt_reject_bonf
    ),
  
  results2 %>%
    filter(true_metric_noninvariance == FALSE) %>%
    mutate(
      estimand = "Type I error",
      method = case_when(
        analysis_form == "linear" ~ "MNLFA linear",
        analysis_form == "quadratic" ~ "MNLFA quadratic"
      ),
      detected = mnlfa_metric_lrt_reject_bonf
    )
  
) %>%
  filter(
    !is.na(method),
    moderator != "noise"
  ) %>%
  group_by(estimand, method, moderator, N) %>%
  summarise(
    rate = mean(detected, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    estimand = factor(
      estimand,
      levels = c("Power", "Type I error")
    ),
    method = factor(
      method,
      levels = c("MNLFA linear", "MNLFA quadratic")
    ),
    moderator = factor(
      moderator,
      levels = c("linear", "quadratic", "sigmoid")
    ),
    N = factor(N)
  )

fig3 <- ggplot(
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
    data = data.frame(estimand = "Type I error"),
    aes(yintercept = .05),
    inherit.aes = FALSE,
    linetype = "dashed",
    colour = "grey60"
  ) +
  geom_line(linewidth = .65, alpha = .75) +
  geom_point(size = 2.2) +
  facet_grid(estimand ~ moderator) +
  scale_colour_manual(values = cols) +
  scale_shape_manual(values = shapes) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    breaks = seq(0, 1, .25)
  ) +
  labs(
    x = "Sample size",
    y = "",
    colour = NULL,
    shape = NULL
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 0)
  )

fig3
# -------------------------------------------------
scalar_perf <- bind_rows(
  
  results2 %>% filter(popmodel=="1.11") %>%
    filter(true_scalar_noninvariance == TRUE) %>%
    mutate(
      estimand = "Power",
      method = case_when(
        analysis_form == "linear" ~ "MNLFA linear",
        analysis_form == "quadratic" ~ "MNLFA quadratic"
      ),
      detected = mnlfa_scalar_lrt_reject_bonf
    ),
  
  results2 %>% filter(popmodel=="1.11") %>%
    filter(true_scalar_noninvariance == FALSE) %>%
    mutate(
      estimand = "Type I error",
      method = case_when(
        analysis_form == "linear" ~ "MNLFA linear",
        analysis_form == "quadratic" ~ "MNLFA quadratic"
      ),
      detected = mnlfa_scalar_lrt_reject_bonf
    )
  
) %>%
  filter(
    !is.na(method),
    moderator != "noise"
  ) %>%
  group_by(estimand, method, moderator, N) %>%
  summarise(
    rate = mean(detected, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    estimand = factor(
      estimand,
      levels = c("Power", "Type I error")
    ),
    method = factor(
      method,
      levels = c("MNLFA linear", "MNLFA quadratic")
    ),
    moderator = factor(
      moderator,
      levels = c("linear", "quadratic", "sigmoid")
    ),
    N = factor(N)
  )


fig3b <- ggplot(
  scalar_perf,
  aes(
    x = N,
    y = rate,
    colour = method,
    shape = method,
    group = method
  )
) +
  geom_hline(
    data = data.frame(estimand = "Type I error"),
    aes(yintercept = .05),
    inherit.aes = FALSE,
    linetype = "dashed",
    colour = "grey60"
  ) +
  geom_line(linewidth = .65, alpha = .75) +
  geom_point(size = 2.2) +
  facet_grid(estimand ~ moderator) +
  scale_colour_manual(values = cols) +
  scale_shape_manual(values = shapes) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    breaks = seq(0, 1, .25)
  ) +
  labs(
    x = "Sample size",
    y = "",
    colour = NULL,
    shape = NULL
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 0)
  )

fig3b
# --------------------
fig4_dat <- bind_rows(
  results2 %>%
    filter(true_metric_noninvariance == FALSE) %>%
    mutate(method = "MNLFA linear", detected = mnlfa_metric_lrt_reject_bonf) %>%
    filter(analysis_form == "linear"),
  
  results2 %>%
    filter(true_metric_noninvariance == FALSE) %>%
    mutate(method = "MNLFA quadratic", detected = mnlfa_metric_lrt_reject_bonf) %>%
    filter(analysis_form == "quadratic"),
  
  results2 %>%
    filter(true_metric_noninvariance == FALSE) %>%
    mutate(method = "SEM Tree", detected = tree_metric_reject_bonf)
) %>%
  filter(moderator != "noise") %>%
  group_by(N, reliability, moderator, method) %>%
  summarise(
    rate = mean(detected, na.rm = TRUE),
    n = sum(!is.na(detected)),
    se = sqrt(rate * (1 - rate) / n),
    lo = pmax(0, rate - 1.96 * se),
    hi = pmin(1, rate + 1.96 * se),
    .groups = "drop"
  )

fig4 <- ggplot(
  fig4_dat,
  aes(
    factor(N),
    rate,
    colour = method,
    fill = method,
    shape = method,
    group = method
  )
) +
  geom_hline(yintercept = .05, linetype = "dashed", color = "grey70") +
  geom_ribbon(
    aes(ymin = lo, ymax = hi),
    alpha = .10,
    colour = NA,
    show.legend = TRUE
  ) +
  geom_line(linewidth = .65, alpha = .75, position = position_dodge(width = .15)) +
  geom_point(size = 2.1, position = position_dodge(width = .15)) +
  facet_grid(
    moderator ~ reliability,
    labeller = labeller(
      moderator = c(
        linear = "Linear",
        quadratic = "Quadratic",
        sigmoid = "Sigmoid"
      ),
      reliability = label_value
    )
  ) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_shape_manual(values = shapes) +
  guides(fill = "none") +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, .80),
    breaks = c(0, .20, .40, .60, .80)
  ) +
  labs(
    x = "Sample size",
    y = "Metric false-positive rate",
    colour = NULL,
    fill = NULL,
    shape = NULL
  ) +
  theme_pub()

fig4
#----------------
fig4b_dat <- bind_rows(
  results2 %>%
    filter(true_scalar_noninvariance == FALSE) %>%
    mutate(method = "MNLFA linear", detected = mnlfa_scalar_lrt_reject_bonf) %>%
    filter(analysis_form == "linear"),
  
  results2 %>%
    filter(true_scalar_noninvariance == FALSE) %>%
    mutate(method = "MNLFA quadratic", detected = mnlfa_scalar_lrt_reject_bonf) %>%
    filter(analysis_form == "quadratic"),
  
  results2 %>%
    filter(true_scalar_noninvariance == FALSE) %>%
    mutate(method = "SEM Tree", detected = tree_scalar_reject_bonf)
) %>%
  filter(moderator != "noise") %>%
  group_by(N, reliability, moderator, method) %>%
  summarise(
    rate = mean(detected, na.rm = TRUE),
    n = sum(!is.na(detected)),
    se = sqrt(rate * (1 - rate) / n),
    lo = pmax(0, rate - 1.96 * se),
    hi = pmin(1, rate + 1.96 * se),
    .groups = "drop"
  )

fig4b <- ggplot(
  fig4b_dat,
  aes(
    factor(N),
    rate,
    colour = method,
    fill = method,
    shape = method,
    group = method
  )
) +
  geom_hline(yintercept = .05, linetype = "dashed", color = "grey70") +
  geom_ribbon(
    aes(ymin = lo, ymax = hi),
    alpha = .10,
    colour = NA,
    show.legend = TRUE
  ) +
  geom_line(linewidth = .65, alpha = .75, position = position_dodge(width = .15)) +
  geom_point(size = 2.1, position = position_dodge(width = .15)) +
  facet_grid(
    moderator ~ reliability,
    labeller = labeller(
      moderator = c(
        linear = "Linear",
        quadratic = "Quadratic",
        sigmoid = "Sigmoid"
      ),
      reliability = label_value
    )
  ) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_shape_manual(values = shapes) +
  guides(fill = "none") +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, .80),
    breaks = c(0, .20, .40, .60, .80)
  ) +
  labs(
    x = "Sample size",
    y = "Scalar false-positive rate",
    colour = NULL,
    fill = NULL,
    shape = NULL
  ) +
  theme_pub()

fig4b
# -------------------------------

metric_summary_dat <- bind_rows(
  results2 %>%
    filter(true_metric_noninvariance == TRUE) %>%
    mutate(
      method = case_when(
        analysis_form == "linear" ~ "MNLFA linear",
        analysis_form == "quadratic" ~ "MNLFA quadratic"
      ),
      detected = mnlfa_metric_lrt_reject_bonf
    ) %>%
    filter(!is.na(method)),
  
  results2 %>%
    filter(true_metric_noninvariance == TRUE) %>%
    mutate(
      method = "SEM Tree",
      detected = semtree_metric_power_decision_bonf
    )
) %>%
  group_by(popmodel, moderator, method) %>%
  summarise(
    power = mean(detected, na.rm = TRUE),
    n = sum(!is.na(detected)),
    .groups = "drop"
  ) %>%
  mutate(
    popmodel = factor(
      popmodel,
      levels = c("1.1", "1.12", "1.2", "1.22", "1.3", "1.32")
    ),
    method = factor(
      method,
      levels = c("MNLFA linear", "MNLFA quadratic", "SEM Tree")
    )
  )

scalar_summary_dat <- bind_rows(
  results2 %>%
    filter(true_scalar_noninvariance == TRUE) %>%
    mutate(
      method = case_when(
        analysis_form == "linear" ~ "MNLFA linear",
        analysis_form == "quadratic" ~ "MNLFA quadratic"
      ),
      detected = mnlfa_scalar_lrt_reject_bonf
    ) %>%
    filter(!is.na(method)),
  
  results2 %>%
    filter(true_scalar_noninvariance == TRUE) %>%
    mutate(
      method = "SEM Tree",
      detected = semtree_scalar_power_decision_bonf
    )
) %>%
  group_by(popmodel, moderator, method) %>%
  summarise(
    power = mean(detected, na.rm = TRUE),
    n = sum(!is.na(detected)),
    .groups = "drop"
  ) %>%
  mutate(
    popmodel = factor(
      popmodel,
      levels = c("1.11", "1.12", "1.21", "1.22", "1.32")
    ),
    method = factor(
      method,
      levels = c("MNLFA linear", "MNLFA quadratic", "SEM Tree")
    )
  )

# -----------------------------------------------------------------------
fig_metric_popmodel <- ggplot(
  metric_summary_dat,
  aes(
    x = method,
    y = power,
    colour = method,
    shape = method
  )
) +
  geom_hline(
    yintercept = .05,
    linetype = "dashed",
    colour = "grey70"
  ) +
  geom_point(
    size = 2.4,
    alpha = .85
  ) +
  facet_grid(
    moderator ~ popmodel,
    labeller = labeller(
      moderator = c(
        linear = "Linear",
        quadratic = "Quadratic",
        sigmoid = "Sigmoid",
        noise = "Noise"
      )
    )
  ) +
  scale_colour_manual(values = cols) +
  scale_shape_manual(values = shapes) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    breaks = seq(0, 1, .25)
  ) +
  labs(
    x = NULL,
    y = "Average metric power",
    colour = NULL,
    shape = NULL
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank()
  )

fig_metric_popmodel
######

fig_metric_heat <- ggplot(
  metric_summary_dat,
  aes(
    x = popmodel,
    y = method,
    fill = power
  )
) +
  geom_tile(
    color = NA,
    linewidth = 0 #no white borders between the cells
  ) +
  geom_text(
    aes(label = scales::percent(power, accuracy = 1)),
    size = 3,
    colour = "black"
  ) +
  facet_wrap(~ moderator, nrow = 2) +
  scale_fill_viridis_c(
    option = "G",
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1),
    name = "Power"
  ) +
  labs(
    x = "Population model",
    y = NULL,
    fill = "Power"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_blank(),
    legend.position = c(.85, .02),
    legend.justification = c(1, 0),
    panel.grid = element_blank(),
    panel.border = element_blank()   # removes frame
  )


fig_metric_heat
##############

all_popmodels <- c("1.1", "1.11", "1.12", "1.3", "1.2", "1.21", "1.22", "1.32")

heat_dat <- bind_rows(
  metric_summary_dat %>%
    transmute(
      popmodel,
      moderator,
      method,
      power,
      invariance = "Metric"
    ),
  
  scalar_summary_dat %>%
    transmute(
      popmodel,
      moderator,
      method,
      power,
      invariance = "Scalar"
    )
) %>%
  filter(moderator != "noise") %>%
  mutate(
    popmodel = factor(popmodel, levels = all_popmodels),
    moderator = factor(
      moderator,
      levels = c("linear", "quadratic", "sigmoid"),
      labels = c("Linear", "Quadratic", "Sigmoid")
    ),
    method = factor(
      method,
      levels = c("SEM Tree", "MNLFA quadratic", "MNLFA linear")
    ),
    invariance = factor(
      invariance,
      levels = c("Metric", "Scalar")
    )
  )
heat_dat_plot <- heat_dat %>%
  group_by(popmodel, moderator, method) %>%
  summarise(
    power = mean(power, na.rm = TRUE),
    .groups = "drop"
  )

fig_power_heat_popmodel <- ggplot(
  heat_dat_plot,
  aes(
    x = moderator,
    y = method,
    fill = power
  )
) +
  geom_tile(
    color = NA,
    linewidth = 0
  ) +
  geom_text(
    aes(label = scales::percent(power, accuracy = 1)),
    size = 3,
    colour = "black"
  ) +
  facet_wrap(
    ~ popmodel,
    nrow = 2,
    drop = FALSE
  ) + # Hokusai2 = list(c("#abc9c8", "#72aeb6", "#4692b0", "#2f70a1", "#134b73", "#0a3351"), c(5, 2, 4, 1, 6, 3), colorblind=TRUE),
  scale_fill_gradientn(
    colours = c(
      "#aadce0",
      "#abc9c8",
      "#72aeb6",
      "#4692b0",
      "#2f70a1",
      "#134b73",
      "#0a3351"
    ),
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1),
    name = "Statistical Power",
    na.value = "white"
  ) +
  labs(
    x = "Moderator form",
    y = NULL
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = .5),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "bottom"
  )

fig_power_heat_popmodel
# -------------------------------

results %>%
  group_by(N, reliability, popmodel, moderator, analysis_form) %>%
  summarise(
    n = n(),
    mnlfa_metric_power = mean(mnlfa_metric_lrt_reject, na.rm = TRUE),
    mnlfa_scalar_power = mean(mnlfa_scalar_lrt_reject, na.rm = TRUE),
    tree_metric_power  = mean(tree_metric_reject, na.rm = TRUE),
    tree_scalar_power  = mean(tree_scalar_reject, na.rm = TRUE),
    .groups = "drop"
  )

results %>%
  group_by(N, true_any_noninvariance) %>%
  summarise(
    mnlfa_power = mean(mnlfa_det, na.rm = TRUE),
    tree_metric_power = mean(tree_metric_reject, na.rm = TRUE),
    tree_scalar_power = mean(tree_scalar_reject, na.rm = TRUE),
    .groups = "drop"
  )



metric_power_N <- results2 %>%
  filter(true_metric_noninvariance == TRUE) %>%
  mutate(
    method = case_when(
      analysis_form == "linear" ~ "MNLFA linear",
      analysis_form == "quadratic" ~ "MNLFA quadratic"
    ),
    detected = mnlfa_metric_lrt_reject_bonf
  ) %>%
  filter(!is.na(method), !is.na(detected)) %>%
  group_by(
    N, reliability, moderator,
    popmodel, delta_lambda, delta_nu,
    method
  ) %>%
  summarise(
    power = mean(detected),
    n = n(),
    .groups = "drop"
  )

metric_power_N_avg <- metric_power_N %>%
  group_by(N, reliability, moderator, method) %>%
  summarise(
    power = mean(power, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(metric_power_N_avg,
       aes(N, power, colour = method, shape = method, group = method)) +
  geom_line() +
  geom_point(size = 2) +
  facet_grid(moderator ~ reliability) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  labs(x = "Sample size", y = "Average metric power") +
  theme_pub()
# -------------------------------

metric_power_form <- bind_rows(
  results2 %>%
    filter(true_metric_noninvariance == TRUE) %>%
    mutate(
      method = case_when(
        analysis_form == "linear" ~ "MNLFA linear",
        analysis_form == "quadratic" ~ "MNLFA quadratic"
      ),
      detected = mnlfa_metric_lrt_reject_bonf
    ) %>%
    filter(!is.na(method)),
  
  results2 %>%
    filter(true_metric_noninvariance == TRUE) %>%
    mutate(
      method = "SEM Tree",
      detected = semtree_metric_power_decision_bonf
    )
) %>%
  filter(moderator != "noise") %>%
  group_by(moderator, method) %>%
  summarise(
    power = mean(detected, na.rm = TRUE),
    n = sum(!is.na(detected)),
    se = sqrt(power * (1 - power) / n),
    lo = pmax(0, power - 1.96 * se),
    hi = pmin(1, power + 1.96 * se),
    .groups = "drop"
  ) %>%
  mutate(
    moderator = factor(
      moderator,
      levels = c("linear", "quadratic", "sigmoid"),
      labels = c("Linear DGM", "Quadratic DGM", "Sigmoid DGM")
    ),
    method = factor(
      method,
      levels = c("MNLFA linear", "MNLFA quadratic", "SEM Tree")
    )
  )

ggplot(
  metric_power_form,
  aes(x = moderator, y = power, colour = method, shape = method, group = method)
) +
  geom_point(size = 3, position = position_dodge(width = .35)) +
  geom_errorbar(
    aes(ymin = lo, ymax = hi),
    width = .10,
    linewidth = .6,
    position = position_dodge(width = .35)
  ) +
  scale_colour_manual(values = cols) +
  scale_shape_manual(values = shapes) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    breaks = seq(0, 1, .25)
  ) +
  labs(
    x = NULL,
    y = "Metric-level power",
    colour = NULL,
    shape = NULL
  ) +
  theme_pub(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 0),
    legend.position = "bottom"
  )

# -------------------


metric_power_mis <- bind_rows(
  results2 %>%
    filter(true_metric_noninvariance == TRUE) %>%
    mutate(
      method = case_when(
        analysis_form == "linear" ~ "MNLFA linear",
        analysis_form == "quadratic" ~ "MNLFA quadratic"
      ),
      detected = mnlfa_metric_lrt_reject_bonf
    ) %>%
    filter(!is.na(method)),
  
  results2 %>%
    filter(true_metric_noninvariance == TRUE) %>%
    mutate(
      method = "SEM Tree",
      detected = semtree_metric_power_decision_bonf
    )
) %>%
  filter(moderator != "noise") %>%
  group_by(moderator, method) %>%
  summarise(
    power = mean(detected, na.rm = TRUE),
    n = sum(!is.na(detected)),
    se = sqrt(power * (1 - power) / n),
    lo = pmax(0, power - 1.96 * se),
    hi = pmin(1, power + 1.96 * se),
    .groups = "drop"
  ) %>%
  mutate(
    moderator = factor(
      moderator,
      levels = c("linear", "quadratic", "sigmoid"),
      labels = c("Linear DGM", "Quadratic DGM", "Sigmoid DGM")
    ),
    method = factor(
      method,
      levels = c("MNLFA linear", "MNLFA quadratic", "SEM Tree")
    )
  )

ggplot(
  metric_power_mis,
  aes(x = method, y = power, fill = method)
) +
  geom_col(width = .65) +
  geom_errorbar(
    aes(ymin = lo, ymax = hi),
    width = .15,
    linewidth = .55
  ) +
  facet_wrap(~ moderator, nrow = 1) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    breaks = seq(0, 1, .25)
  ) +
  labs(
    x = NULL,
    y = "Metric-level power",
    fill = NULL
  ) +
  theme_pub(base_size = 13) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
