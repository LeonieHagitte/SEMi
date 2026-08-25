library(dplyr)
library(ggplot2)
library(ggalluvial)
library(stringr)

tree_flow_dat <- results2 %>%
  filter(moderator != "noise") %>%
  mutate(
    truth = case_when(
      true_metric_noninvariance == TRUE & true_scalar_noninvariance == TRUE ~ "Metric & Scalar",
      true_metric_noninvariance == TRUE & true_scalar_noninvariance == FALSE ~ "Metric only",
      true_metric_noninvariance == FALSE & true_scalar_noninvariance == TRUE ~ "Scalar only",
      true_metric_noninvariance == FALSE & true_scalar_noninvariance == FALSE ~ "Invariant",
      TRUE ~ NA_character_
    ),
    
    metric_tree_decision = case_when(
      tree_metric_reject_bonf == TRUE ~ "Metric tree rejected",
      tree_metric_reject_bonf == FALSE ~ "Metric tree not rejected",
      is.na(tree_metric_reject_bonf) ~ "Metric tree missing"
    ),
    
    scalar_tree_decision = case_when(
      tree_metric_reject_bonf == TRUE ~ "Scalar tree not attempted",
      tree_metric_reject_bonf == FALSE & tree_scalar_reject_bonf == TRUE ~ "Scalar tree rejected",
      tree_metric_reject_bonf == FALSE & tree_scalar_reject_bonf == FALSE ~ "Scalar tree not rejected",
      is.na(tree_scalar_reject_bonf) ~ "Scalar tree missing",
      TRUE ~ "Other"
    ),
    
    final_tree_outcome = case_when(
      true_metric_noninvariance == TRUE &
        tree_metric_reject_bonf == TRUE &
        tree_metric_correct_split == TRUE ~ "Correct metric detection",
      
      true_scalar_noninvariance == TRUE &
        tree_metric_reject_bonf == FALSE &
        tree_scalar_reject_bonf == TRUE &
        tree_scalar_correct_split == TRUE ~ "Correct scalar detection",
      
      true_metric_noninvariance == FALSE &
        tree_metric_reject_bonf == TRUE ~ "Metric false positive",
      
      true_scalar_noninvariance == FALSE &
        tree_metric_reject_bonf == FALSE &
        tree_scalar_reject_bonf == TRUE ~ "Scalar false positive",
      
      true_any_noninvariance == TRUE &
        (tree_metric_reject_bonf == FALSE | is.na(tree_metric_reject_bonf)) &
        (tree_scalar_reject_bonf == FALSE | is.na(tree_scalar_reject_bonf)) ~ "Missed noninvariance",
      
      true_any_noninvariance == FALSE &
        tree_metric_reject_bonf == FALSE &
        tree_scalar_reject_bonf == FALSE ~ "Correct invariance decision",
      
      TRUE ~ "Missing"
    )
  ) %>%
  count(truth, metric_tree_decision, scalar_tree_decision, final_tree_outcome, name = "n") %>%
  mutate(prop = n / sum(n))


#------------------------------------------------------------------------------
truth_cols <- c(
  "Metric & Scalar" = "#b38711",
  "Metric only"     = "#af4f2f",
  "Scalar only"     = "#732f30",
  "Invariant"       = "#59385c"
)
axis_levels <- c(
  "truth",
  "metric_tree_decision",
  "scalar_tree_decision",
  "final_tree_outcome"
)

axis_labels <- c(
  truth = "Truth",
  metric_tree_decision = "Metric tree",
  scalar_tree_decision = "Scalar tree",
  final_tree_outcome = "Final outcome"
)

#Redon = list(c("#5b859e", "#1e395f", "#75884b", "#1e5a46", "#df8d71", "#af4f2f", "#d48f90", "#732f30", "#ab84a5", "#59385c", "#d8b847", "#b38711"), c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), colorblind=FALSE),

tree_flow_plot <- tree_flow_dat %>%
  mutate(
    alluvium_id = row_number(),
    flow_truth = truth
  )

tree_lodes <- to_lodes_form(
  tree_flow_plot,
  axes = all_of(axis_levels),
  id = "alluvium_id",
  key = "axis",
  value = "stratum"
) %>%
  mutate(
    axis = factor(axis, levels = axis_levels, labels = axis_labels)
  )

fig_tree_flow <- ggplot(
  tree_lodes,
  aes(
    x = axis,
    stratum = stratum,
    alluvium = alluvium_id,
    y = n
  )
) +
  
  ## ribbons first, so nodes are drawn on top and cover ribbon ends
  geom_alluvium(
    aes(fill = flow_truth),
    alpha = 0.60,
    width = 0.24,
    knot.pos = 0.45,
    curve_type = "sigmoid"
  ) +
  
  ## grey decision / outcome nodes
  geom_stratum(
    data = tree_lodes %>% filter(axis != "Truth"),
    fill = "#dcdcdc",
    colour = "white",
    linewidth = 1,
    width = 0.10,
    alpha = 1
  ) +
  
  ## colored truth nodes
  geom_stratum(
    data = tree_lodes %>% filter(axis == "Truth"),
    aes(fill = stratum),
    colour = "white",
    linewidth = 1,
    width = 0.24,
    alpha = 1
  ) +
  
  ## labels inside first column
  geom_text(
    data = tree_lodes %>% filter(axis == "Truth"),
    stat = "stratum",
    aes(label = after_stat(str_wrap(stratum, width = 14))),
    hjust = 0.5,
    size = 3.2,
    fontface = "bold",
    lineheight = 0.9,
    colour = "black"
  ) +
  
  ## labels outside middle and final columns
  geom_text(
    data = tree_lodes %>% filter(axis != "Truth"),
    stat = "stratum",
    aes(label = after_stat(str_wrap(stratum, width = 20))),
    hjust = 0,
    nudge_x = 0.09,
    size = 3.0,
    lineheight = 0.9,
    colour = "black"
  ) +
  
  scale_x_discrete(
    expand = expansion(mult = c(0.00, 0.00))
  ) +
  
  scale_fill_manual(
    values = truth_cols,
    breaks = names(truth_cols)
  ) +
  
  labs(
    x = NULL,
    y = NULL,
    fill = NULL
  ) +
  
  coord_cartesian(clip = "off") +
  
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    plot.margin = margin(10, 90, 10, 10)
  )

fig_tree_flow
#------------------------------------------

tree_flow_plot_perc <- tree_flow_dat %>%
  mutate(
    alluvium_id = row_number(),
    flow_truth = truth
  )

tree_lodes <- to_lodes_form(
  tree_flow_plot_perc,
  axes = all_of(axis_levels),
  id = "alluvium_id",
  key = "axis",
  value = "stratum"
) %>%
  mutate(
    axis = factor(axis, levels = axis_levels, labels = axis_labels)
  )

stratum_props <- tree_lodes %>%
  group_by(axis, stratum) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  group_by(axis) %>%
  mutate(
    pct = n / sum(n),
    pct_label = scales::percent(pct, accuracy = 0.1)
  ) %>%
  ungroup()

fig_tree_flow_perc <- ggplot(
  tree_lodes,
  aes(
    x = axis,
    stratum = stratum,
    alluvium = alluvium_id,
    y = n
  )
) +
  
  geom_alluvium(
    aes(fill = flow_truth),
    alpha = 0.60,
    width = 0.24,
    knot.pos = 0.45,
    curve_type = "sigmoid",
    colour = NA
  ) +
  
  geom_stratum(
    data = tree_lodes %>% filter(axis != "Truth"),
    fill = "#dcdcdc",
    colour = "white",
    linewidth = 1,
    width = 0.10,
    alpha = 1
  ) +
  
  geom_stratum(
    data = tree_lodes %>% filter(axis == "Truth"),
    aes(fill = stratum),
    colour = "white",
    linewidth = 1,
    width = 0.24,
    alpha = 1
  ) +
  
  ## Truth labels: category + percentage inside first column
  geom_text(
    data = stratum_props %>% filter(axis == "Truth"),
    inherit.aes = FALSE,
    stat = "stratum",
    aes(
      x = axis,
      stratum = stratum,
      y = n,
      label = paste0(
        str_wrap(stratum, width = 14),
        "\n",
        pct_label
      )
    ),
    hjust = 0.5,
    size = 3.0,
    fontface = "bold",
    lineheight = 0.9,
    colour = "black"
  ) +
  
  ## Percentages only inside decision / outcome panels
  geom_text(
    data = stratum_props %>% filter(axis != "Truth"),
    inherit.aes = FALSE,
    stat = "stratum",
    aes(
      x = axis,
      stratum = stratum,
      y = n,
      label = pct_label
    ),
    size = 2.7,
    fontface = "bold",
    colour = "grey25"
  ) +
  ## Category labels outside decision / outcome panels
  geom_text(
    data = stratum_props %>% filter(axis != "Truth"),
    inherit.aes = FALSE,
    stat = "stratum",
    aes(
      x = axis,
      stratum = stratum,
      y = n,
      label = str_wrap(stratum, width = 20)
    ),
    hjust = 0,
    nudge_x = 0.08,
    size = 2.9,
    lineheight = 0.9,
    colour = "black"
  ) +
  
  scale_x_discrete(
    expand = expansion(mult = c(0.00, 0.00))
  ) +
  
  scale_fill_manual(
    values = truth_cols,
    breaks = names(truth_cols)
  ) +
  
  labs(
    x = NULL,
    y = NULL,
    fill = NULL
  ) +
  
  coord_cartesian(clip = "off") +
  
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    plot.margin = margin(10, 90, 10, 10)
  )

fig_tree_flow_perc
