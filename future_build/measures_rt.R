# measures_rt.R

library(tidyverse)
library(dplyr)
library(tidyr)

results_path <- "results.csv"

results <- read.csv(results_path, stringsAsFactors = FALSE) %>%
  tibble::as_tibble()

# ---------------------------
# overall metric performance
# ---------------------------
metric_summary <- results %>%
  summarise(
    mnlfa_metric_type1 = mean(mnlfa_metric_lrt_reject[true_metric_noninvariance == FALSE], na.rm = TRUE),
    mnlfa_metric_power = mean(mnlfa_metric_lrt_reject[true_metric_noninvariance == TRUE], na.rm = TRUE),
    
    tree_metric_type1  = mean(tree_metric_reject[true_metric_noninvariance == FALSE], na.rm = TRUE),
    tree_metric_power  = mean(tree_metric_reject[true_metric_noninvariance == TRUE], na.rm = TRUE)
  )

# ---------------------------
# overall scalar performance
# ---------------------------
scalar_summary <- results %>%
  summarise(
    mnlfa_scalar_type1 = mean(mnlfa_scalar_lrt_reject[true_scalar_noninvariance == FALSE], na.rm = TRUE),
    mnlfa_scalar_power = mean(mnlfa_scalar_lrt_reject[true_scalar_noninvariance == TRUE], na.rm = TRUE),
    
    tree_scalar_type1  = mean(tree_scalar_reject[true_scalar_noninvariance == FALSE], na.rm = TRUE),
    tree_scalar_power  = mean(tree_scalar_reject[true_scalar_noninvariance == TRUE], na.rm = TRUE)
  )

# ---------------------------
# stratified performance
# ---------------------------
condition_summary <- results %>%
  group_by(N, reliability, moderator) %>%
  summarise(
    n_replications = n(),
    
    mnlfa_metric_type1 = mean(mnlfa_metric_lrt_reject[true_metric_noninvariance == FALSE], na.rm = TRUE),
    mnlfa_metric_power = mean(mnlfa_metric_lrt_reject[true_metric_noninvariance == TRUE], na.rm = TRUE),
    
    tree_metric_type1  = mean(tree_metric_reject[true_metric_noninvariance == FALSE], na.rm = TRUE),
    tree_metric_power  = mean(tree_metric_reject[true_metric_noninvariance == TRUE], na.rm = TRUE),
    
    mnlfa_scalar_type1 = mean(mnlfa_scalar_lrt_reject[true_scalar_noninvariance == FALSE], na.rm = TRUE),
    mnlfa_scalar_power = mean(mnlfa_scalar_lrt_reject[true_scalar_noninvariance == TRUE], na.rm = TRUE),
    
    tree_scalar_type1  = mean(tree_scalar_reject[true_scalar_noninvariance == FALSE], na.rm = TRUE),
    tree_scalar_power  = mean(tree_scalar_reject[true_scalar_noninvariance == TRUE], na.rm = TRUE),
    
    tree_scalar_available = mean(!is.na(tree_scalar_reject)),
    
    .groups = "drop"
  )

print(metric_summary)
print(scalar_summary)
print(condition_summary)

######
library(dplyr)
library(tidyr)

full_table <- condition_summary %>%
  select(
    N, reliability, moderator, n_replications,
    
    # metric
    mnlfa_metric_type1, mnlfa_metric_power,
    tree_metric_type1, tree_metric_power,
    
    # scalar
    mnlfa_scalar_type1, mnlfa_scalar_power,
    tree_scalar_type1, tree_scalar_power
  )

long_table <- full_table %>%
  pivot_longer(
    cols = -c(N, reliability, moderator, n_replications),
    names_to = "measure",
    values_to = "value"
  ) %>%
  mutate(
    level  = ifelse(grepl("metric", measure), "Metric", "Scalar"),
    method = ifelse(grepl("mnlfa", measure), "MNLFA", "SEMTREE"),
    type   = ifelse(grepl("type1", measure), "Type I", "Power")
  ) %>%
  select(N, reliability, moderator, level, method, type, value)
######


overall_metric <- metric_summary %>%
  transmute(
    level = "metric",
    mnlfa_type1 = mnlfa_metric_type1,
    mnlfa_power = mnlfa_metric_power,
    tree_type1  = tree_metric_type1,
    tree_power  = tree_metric_power
  )

overall_scalar <- scalar_summary %>%
  transmute(
    level = "scalar",
    mnlfa_type1 = mnlfa_scalar_type1,
    mnlfa_power = mnlfa_scalar_power,
    tree_type1  = tree_scalar_type1,
    tree_power  = tree_scalar_power
  )

overall_table <- bind_rows(overall_metric, overall_scalar) %>%
  pivot_longer(
    cols = -level,
    names_to = "measure",
    values_to = "value"
  ) %>%
  mutate(
    method = ifelse(grepl("^mnlfa", measure), "MNLFA", "SEMTREE"),
    type   = ifelse(grepl("type1", measure), "Type I", "Power")
  ) %>%
  select(level, method, type, value) %>%
  arrange(level, method, type)

overall_table

####

condition_table <- condition_summary %>%
  select(N, reliability, moderator,
         starts_with("mnlfa_"),
         starts_with("tree_")) %>%
  arrange(N, reliability, moderator)

condition_table
#####
library(ggplot2)


ggplot(overall_table, aes(x = method, y = value, fill = type)) +
  geom_col(position = "dodge") +
  facet_wrap(~ level) +
  labs(
    x = NULL,
    y = "Rate"
  ) +
  theme_minimal()

####

metric_plot_data <- condition_summary %>%
  select(N, reliability, moderator,
         mnlfa_metric_type1, mnlfa_metric_power,
         tree_metric_type1, tree_metric_power) %>%
  pivot_longer(
    cols = -c(N, reliability, moderator),
    names_to = "measure",
    values_to = "value"
  ) %>%
  mutate(
    method = ifelse(grepl("mnlfa", measure), "MNLFA", "SEMTREE"),
    stat = ifelse(grepl("type1", measure), "Type I", "Power")
  )

ggplot(metric_plot_data,
       aes(x = factor(N), y = value, fill = method)) +
  geom_col(position = "dodge") +
  facet_grid(stat ~ reliability + moderator) +
  labs(x = "Sample size (N)", y = "Rate") +
  theme_minimal()
####

scalar_plot_data <- condition_summary %>%
  select(N, reliability, moderator,
         mnlfa_scalar_type1, mnlfa_scalar_power,
         tree_scalar_type1, tree_scalar_power) %>%
  pivot_longer(
    cols = -c(N, reliability, moderator),
    names_to = "measure",
    values_to = "value"
  ) %>%
  mutate(
    method = ifelse(grepl("mnlfa", measure), "MNLFA", "SEMTREE"),
    stat = ifelse(grepl("type1", measure), "Type I", "Power")
  )

####
ggplot(scalar_plot_data,
       aes(x = factor(N), y = value, fill = method)) +
  geom_col(position = "dodge") +
  facet_grid(stat ~ reliability + moderator) +
  labs(x = "Sample size (N)", y = "Rate") +
  theme_minimal()
