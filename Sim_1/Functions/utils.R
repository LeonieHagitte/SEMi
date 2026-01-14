# ---------- DGP truth from simulator output -----------------------------------
truth_from_sim <- function(sim, link_label) {
  # sim$true_deltas is a named vector like delta_lambda_y1, ..., delta_theta_y4
  td <- sim$true_deltas
  pick_idx <- function(prefix) {
    nm <- names(td)[startsWith(names(td), prefix)]
    idx <- as.integer(gsub("^.*_y", "", nm)[td[nm] != 0])
    sort(idx)
  }
  list(
    true_link          = link_label,                    # "linear"/"sigmoid"/"quadratic"/"noise"
    true_load_items    = pick_idx("delta_lambda"),
    true_theta_items   = pick_idx("delta_theta"),
    # optional: whether reference y4 was moderated in the DGP simulation
    true_ref_moderated = any(pick_idx("delta_lambda") == 4)
  )
}

# ---------- Parse analysis model syntax to infer assumptions ------------------
# Works for mxsem builders and CFA baseline.
assumptions_from_syntax <- function(model_syntax, method = c("MNLFA","SEMTREE")) {
  method <- match.arg(method)
  
  # helper: return sorted unique item indices from regex matches like "*y1" etc.
  grab_items <- function(re_vec) {
    if (length(re_vec) == 0L) return(integer(0))
    m <- regmatches(model_syntax, gregexpr(re_vec, model_syntax, perl = TRUE))[[1]]
    if (length(m) == 0L) return(integer(0))
    as.integer(sort(unique(gsub("^.*y([1-9][0-9]*)$", "\\1", m))))
  }
  
  if (method == "SEMTREE") {
    # CFA baseline: no parametric moderation; splits are nonparametric by tree.
    return(list(
      model_label   = "SEMTREE_CFA_baseline",
      parametric    = FALSE,
      link          = NA_character_,
      load_items    = integer(0),
      theta_items   = integer(0),
      ref_moderated = FALSE
    ))
  }
  
  # MNLFA (mxsem): detect moderation via presence of "*data.M"
  # Loadings: pieces like "{... delta_lambda_y1*data.M}*y1"
  # Loadings moderation (λ): detect delta_lambda_y<k> * data.<var>
  load_moderated <- {
    m <- regmatches(
      model_syntax,
      gregexpr("delta_lambda_y([0-9]+)\\s*\\*\\s*data\\.[A-Za-z_][A-Za-z0-9_]*", model_syntax, perl = TRUE)
    )[[1]]
    if (length(m) == 0L) integer(0) else as.integer(sort(unique(sub(".*delta_lambda_y([0-9]+).*", "\\1", m))))
  }
  
  # Residual variance moderation (Θ): detect delta_theta_y<k> * data.<var>
  theta_moderated <- {
    m <- regmatches(
      model_syntax,
      gregexpr("delta_theta_y([0-9]+)\\s*\\*\\s*data\\.[A-Za-z_][A-Za-z0-9_]*", model_syntax, perl = TRUE)
    )[[1]]
    if (length(m) == 0L) integer(0) else as.integer(sort(unique(sub(".*delta_theta_y([0-9]+).*", "\\1", m))))
  }
  
  
  # Reference y4 fixed as marker if "1*y4" appears and no "*data.M" with y4 in loadings
  ref_fixed_1   <- grepl("\\b1\\*y4\\b", model_syntax, perl = TRUE)
  ref_moderated <- 4L %in% load_moderated
  
  # Link: builders are linear in M for loadings, log link for Θ when present.
  link_label <- if (length(theta_moderated)) "linear (λ), log-linear (Θ)" else "linear (λ)"
  
  list(
    model_label   = if (length(theta_moderated)) "MNLFA_linear_full" else "MNLFA_linear_partial",
    parametric    = TRUE,
    link          = link_label,
    load_items    = sort(load_moderated),
    theta_items   = sort(theta_moderated),
    ref_moderated = ref_moderated || (!ref_fixed_1 && (4L %in% load_moderated))
  )
}

# ---------- small pretty-printer used in the row tibble -----------------------
.items_str <- function(idx) if (length(idx)) paste0("y", idx, collapse = ",") else "none"
