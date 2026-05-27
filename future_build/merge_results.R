library(tidyverse)

files <- list.files("rds", pattern = "\\.rds$", full.names = TRUE)

stopifnot(length(files) > 0)

results <- data.frame(do.call(rbind, lapply(files, function(f) {
  x <- readRDS(f)
  stopifnot(is.matrix(x))
  x <- as.data.frame(lapply(as.data.frame(x), unlist))
  x
})))

############
# display results
#

# Sanity check on the number of simulated conditions:
# Count number of cases based on some simulation parameters
results %>%
  dplyr::count(popmodel, moderator, true_any_noninvariance)


results %>% ggplot(aes(x=popmodel,fill=moderator,y=tree_metric_p))+  stat_summary(
  fun = mean,
  geom = "bar",
  position = position_dodge(width = 0.9)
)

results %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_errors = sum(!is.na(error_msg)),
    mnlfa_metric_reject_rate = mean(mnlfa_metric_lrt_reject, na.rm = TRUE),
    mnlfa_scalar_reject_rate = mean(mnlfa_scalar_lrt_reject, na.rm = TRUE),
    tree_metric_reject_rate = mean(tree_metric_reject, na.rm = TRUE),
    tree_scalar_reject_rate = mean(tree_scalar_reject, na.rm = TRUE)
  )


