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
#-------------------------------------------------------------------------------
results2_confirmatory <- results2 %>%
  filter(
    num_noisy_predictors == 0,
    delta_lambda %in% c(0.2, 0.3),
    delta_nu %in% c(0.5, 1)
  )

results2_exploratory <- results2 %>%
  filter(
    num_noisy_predictors == 0,
    !(delta_lambda %in% c(0.1, 0.2) &
        delta_nu %in% c(0.25, 0.5))
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

