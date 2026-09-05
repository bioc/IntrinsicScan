# Assisted-by: OpenAI ChatGPT (GPT-5.6 Sol).
# AI assistance included method/code review, package engineering, testing
# strategy, documentation, and scientific drafting. The human authors
# reviewed and validated the implementation and retain responsibility.

# ============================================================================
# LABEL-BLIND SINGLE-FEATURE STATE CONCORDANCE v11
# ============================================================================
# Purpose
# -------
# For each continuous feature in a binary-phenotype dataset:
#   1. Infer a two-state threshold from the feature values alone, without using
#      phenotype labels. The threshold is the exact global minimum of the
#      two-cluster 1D L1 medoid objective over all valid contiguous splits.
#   2. After the expression/state partition is fixed, quantify how strongly it
#      concords with the binary phenotype using orientation-free Youden J.
#   3. Optionally assess INTERNAL STABILITY by repeatedly refitting the
#      label-blind threshold on stratified training resamples, choosing the
#      phenotype orientation from training labels only, and applying the locked
#      threshold/orientation to held-out or out-of-bag observations. Each
#      feature receives a deterministic RNG stream derived from the master seed
#      and feature ID, so its resampling history does not depend on rank/order.
#   4. Provide an exact feature-wise conditional null test for J. Conditional
#      on the observed high-state size and phenotype counts, the number of
#      phenotype-positive observations in the high state is hypergeometric.
#   5. Optionally use label permutations for genome-wide maxT/global-null
#      calibration of the FIXED full-sample label-blind partitions.
#
# What this code does NOT claim
# -----------------------------
#   - It does not claim that exact 1D medoid clustering is itself novel.
#   - It does not optimize thresholds using phenotype labels.
#   - Internal resampling is a stability assessment, not independent external
#     validation and not an unbiased estimate of a full feature-selection
#     pipeline's predictive performance.
#   - Permuting labels after fixing the full-sample partitions calibrates
#     phenotype concordance of those fixed partitions. It does NOT rerun the
#     resampling/feature-selection procedure and therefore is not a permutation
#     test of the "whole pipeline".
#   - No single categorical label is used to merge effect size, stability and
#     inferential evidence. These are returned as separate axes.
#   - Internal stability RNG is feature-specific: reranking/reordering other
#     features cannot silently change a selected feature's resampling draws.
#
# Data model
# ----------
# Rows are generic FEATURES (probes, genes, transcripts, proteins, metabolites,
# CpGs, etc.) and columns are observations/samples. Feature row names and sample
# column names must be unique.
#
# Dependencies: base R / stats only.
# ============================================================================


# ----------------------------------------------------------------------------
# Small utilities
# ----------------------------------------------------------------------------
.validate_seed_v11 <- function(seed, name = "seed") {
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != floor(seed) || seed < 0 || seed > 2147483646) {
    stop(name, " must be one integer between 0 and 2147483646")
  }
  as.integer(seed)
}


# Deterministically derive a feature-specific RNG seed from a master seed and
# the unique feature ID. This deliberately does NOT depend on rank, row order,
# or the number/order of other features that happen to be resampled. Therefore
# the same feature ID + master seed receives the same resampling stream whenever
# its values/labels are unchanged. The hash is deterministic, not cryptographic.
.feature_stability_seed_v11 <- function(master_seed, feature_id) {
  master_seed <- .validate_seed_v11(master_seed, "master seed")
  if (length(feature_id) != 1L || is.na(feature_id) || !nzchar(feature_id)) {
    stop("feature_id must be one non-missing, non-empty string")
  }

  modulus <- 2147483646
  h <- (as.double(master_seed) + 1) %% modulus
  bytes <- utf8ToInt(enc2utf8(as.character(feature_id)))

  if (length(bytes)) {
    for (b in bytes) {
      # All intermediate integers remain far below 2^53, so double-precision
      # arithmetic represents them exactly before the modulo operation.
      h <- (h * 131 + as.double(b) + 17) %% modulus
    }
  }

  as.integer((h %% (modulus - 1)) + 1)
}


.validate_binary_inputs_v11 <- function(feature_matrix, labels, positive_label) {
  if (!is.matrix(feature_matrix) && !is.data.frame(feature_matrix)) {
    stop("feature_matrix must be a matrix or data.frame")
  }

  feature_matrix <- as.matrix(feature_matrix)
  suppressWarnings(storage.mode(feature_matrix) <- "double")

  if (!is.numeric(feature_matrix)) {
    stop("feature_matrix must be coercible to a numeric matrix")
  }
  if (nrow(feature_matrix) < 1L || ncol(feature_matrix) < 2L) {
    stop("feature_matrix must contain at least 1 feature and 2 samples")
  }
  if (is.null(rownames(feature_matrix))) {
    stop("feature_matrix must have feature row names")
  }
  if (is.null(colnames(feature_matrix))) {
    stop("feature_matrix must have sample column names")
  }
  if (anyNA(rownames(feature_matrix)) || any(!nzchar(rownames(feature_matrix)))) {
    stop("feature row names must be non-missing and non-empty")
  }
  if (anyNA(colnames(feature_matrix)) || any(!nzchar(colnames(feature_matrix)))) {
    stop("sample column names must be non-missing and non-empty")
  }
  if (anyDuplicated(rownames(feature_matrix))) {
    stop("feature row names must be unique; duplicated feature IDs are not allowed")
  }
  if (anyDuplicated(colnames(feature_matrix))) {
    stop("sample column names must be unique")
  }
  if (anyNA(feature_matrix) || any(!is.finite(feature_matrix))) {
    stop("v11 currently requires a complete finite feature matrix (no NA/Inf)")
  }

  if (length(labels) != ncol(feature_matrix)) {
    stop("labels length must equal ncol(feature_matrix)")
  }
  if (anyNA(labels)) {
    stop("labels must not contain NA")
  }

  labels_chr <- as.character(labels)
  classes <- unique(labels_chr)
  if (length(classes) != 2L) {
    stop("v11 supports exactly two phenotype classes")
  }

  if (length(positive_label) != 1L || is.na(positive_label)) {
    stop("positive_label must be explicitly supplied as one non-missing value")
  }
  positive_label <- as.character(positive_label)
  if (!positive_label %in% classes) {
    stop("positive_label is not present in labels")
  }

  negative_label <- setdiff(classes, positive_label)
  # Force a deterministic two-level factor while preserving the user's stated
  # positive class rather than relying on alphabetical factor-level order.
  labels_factor <- factor(labels_chr, levels = c(negative_label, positive_label))

  class_counts <- table(labels_factor)
  if (any(class_counts < 2L)) {
    stop("each phenotype class must contain at least 2 samples")
  }
  if (min(class_counts) < 5L) {
    warning(
      "At least one phenotype class has fewer than 5 samples. ",
      "Effect sizes can still be calculated, but held-out stability estimates ",
      "will be highly discrete and inferential power may be limited.",
      call. = FALSE
    )
  }

  list(
    feature_matrix = feature_matrix,
    labels = labels_factor,
    positive_label = positive_label,
    negative_label = negative_label,
    n_positive = unname(class_counts[positive_label]),
    n_negative = unname(class_counts[negative_label])
  )
}

.quantile_finite_v11 <- function(z, prob) {
  z <- z[is.finite(z)]
  if (length(z) == 0L) return(NA_real_)
  as.numeric(stats::quantile(z, probs = prob, names = FALSE, type = 7))
}


# ----------------------------------------------------------------------------
# Exact globally optimal 1D two-medoids partition under L1 loss
# ----------------------------------------------------------------------------
# In one dimension, an optimum of the two-cluster L1 medoid objective can be
# represented by two contiguous groups after sorting. Therefore every valid
# split between distinct adjacent values can be enumerated exactly.
#
# Interval L1 costs are evaluated in O(1) after sorting using prefix sums, so
# the per-feature complexity is O(N log N) for sorting plus O(N) split scans.
# Labels are not an argument and cannot influence the fitted threshold.
# ----------------------------------------------------------------------------
exact_two_medoids_1d_v11 <- function(x) {
  x <- as.numeric(x)

  if (length(x) < 2L) {
    return(list(valid = FALSE, reason = "fewer than 2 observations"))
  }
  if (anyNA(x) || any(!is.finite(x))) {
    return(list(valid = FALSE, reason = "non-finite values"))
  }
  if (length(unique(x)) < 2L) {
    return(list(valid = FALSE, reason = "zero-variance feature"))
  }

  ord <- order(x)
  xs <- x[ord]
  n <- length(xs)

  # Thresholds are considered only between distinct observed values.
  candidate_splits <- which(diff(xs) > 0)
  if (length(candidate_splits) == 0L) {
    return(list(valid = FALSE, reason = "no valid split"))
  }

  cs <- c(0, cumsum(xs))

  segment_l1_cost <- function(l, r) {
    # Any median minimizes L1 loss. For even segments use the lower median
    # deterministically; the minimum L1 cost is unchanged by this tie choice.
    m <- floor((l + r) / 2)
    med <- xs[m]

    left_sum <- cs[m + 1L] - cs[l]
    right_sum <- cs[r + 1L] - cs[m + 1L]

    left_cost <- med * (m - l + 1L) - left_sum
    right_cost <- right_sum - med * (r - m)
    left_cost + right_cost
  }

  costs <- vapply(candidate_splits, function(s) {
    segment_l1_cost(1L, s) + segment_l1_cost(s + 1L, n)
  }, numeric(1))

  min_cost <- min(costs)
  tol <- sqrt(.Machine$double.eps) * max(1, abs(min_cost))
  best <- candidate_splits[costs <= min_cost + tol]

  # Deterministic and label-blind tie breaking:
  #   1) prefer the largest empty interval between the two states;
  #   2) if still tied, prefer the split closest to the middle of the ordering.
  if (length(best) > 1L) {
    gaps <- xs[best + 1L] - xs[best]
    best <- best[gaps == max(gaps)]
  }
  if (length(best) > 1L) {
    best <- best[which.min(abs(best - n / 2))]
  }

  s <- best[1L]
  threshold <- xs[s] + (xs[s + 1L] - xs[s]) / 2
  low_medoid <- xs[floor((1L + s) / 2)]
  high_medoid <- xs[floor((s + 1L + n) / 2)]
  gap <- xs[s + 1L] - xs[s]

  # Because threshold lies strictly between distinct observations, these state
  # assignments reproduce the chosen contiguous split exactly.
  high_state <- x > threshold

  iqr_x <- stats::IQR(x)
  relative_gap <- if (is.finite(iqr_x) && iqr_x > 0) gap / iqr_x else NA_real_

  list(
    valid = TRUE,
    threshold = threshold,
    high_state = high_state,
    low_medoid = low_medoid,
    high_medoid = high_medoid,
    medoid_distance = high_medoid - low_medoid,
    gap = gap,
    relative_gap = relative_gap,
    within_l1_cost = min_cost,
    split_index_sorted = s,
    sorted_order = ord
  )
}


# ----------------------------------------------------------------------------
# Phenotype concordance after the label-blind partition is fixed
# ----------------------------------------------------------------------------
# The two expression states have no phenotype meaning during fitting. Labels
# are revealed only here. Orientation-free Youden J is
#
#   J = |TPR_high - FPR_high| = sensitivity + specificity - 1
#
# after choosing whether phenotype-positive corresponds to the high or low
# state. Thus J lies in [0, 1], with J=1 meaning exact full-sample concordance.
# ----------------------------------------------------------------------------
score_fixed_partition_v11 <- function(high_state, labels, positive_label) {
  labels_chr <- as.character(labels)
  positive_label <- as.character(positive_label)

  if (length(high_state) != length(labels_chr)) {
    stop("high_state and labels lengths differ")
  }
  if (length(unique(labels_chr)) != 2L) {
    stop("exactly two phenotype classes are required")
  }
  if (!positive_label %in% labels_chr) {
    stop("positive_label not found in labels")
  }

  y <- labels_chr == positive_label
  n_pos <- sum(y)
  n_neg <- sum(!y)
  if (n_pos == 0L || n_neg == 0L) {
    stop("both phenotype classes must be present")
  }

  high_state <- as.logical(high_state)
  tpr_high <- sum(y & high_state) / n_pos
  fpr_high <- sum(!y & high_state) / n_neg
  signed_J_high <- tpr_high - fpr_high

  if (signed_J_high >= 0) {
    positive_state <- "high"
    predicted_positive <- high_state
  } else {
    positive_state <- "low"
    predicted_positive <- !high_state
  }

  sensitivity <- sum(predicted_positive & y) / n_pos
  specificity <- sum(!predicted_positive & !y) / n_neg
  J <- sensitivity + specificity - 1

  tp <- sum(predicted_positive & y)
  fn <- sum(!predicted_positive & y)
  fp <- sum(predicted_positive & !y)
  tn <- sum(!predicted_positive & !y)

  tab <- table(
    state = factor(high_state, levels = c(FALSE, TRUE), labels = c("low", "high")),
    phenotype = factor(y, levels = c(FALSE, TRUE), labels = c("negative", "positive"))
  )

  fisher_p <- tryCatch(
    stats::fisher.test(tab, alternative = "two.sided")$p.value,
    error = function(e) NA_real_
  )

  list(
    J = as.numeric(J),
    signed_J_high = as.numeric(signed_J_high),
    positive_state = positive_state,
    sensitivity = as.numeric(sensitivity),
    specificity = as.numeric(specificity),
    balanced_accuracy = as.numeric((sensitivity + specificity) / 2),
    accuracy = as.numeric(mean(predicted_positive == y)),
    tp = tp,
    fn = fn,
    fp = fp,
    tn = tn,
    fisher_p = fisher_p,
    table = tab
  )
}


# ----------------------------------------------------------------------------
# Exact conditional null test aligned directly to orientation-free J
# ----------------------------------------------------------------------------
# Once the label-blind state partition is fixed, condition on:
#   n_positive, n_negative, and high_state_size.
# Under exchangeable phenotype labels, the number K of positive samples in the
# high state follows Hypergeometric(n_positive, n_negative, high_state_size).
#
# We sum hypergeometric probabilities for all feasible K producing J >= J_obs.
# This gives an exact finite-sample p-value for the observed fixed-partition J
# and avoids the p-value resolution problem of using only ~1000 permutations
# for tens of thousands of features.
# ----------------------------------------------------------------------------
exact_J_hypergeom_p_v11 <- function(high_state_size,
                                    n_positive,
                                    n_negative,
                                    observed_J,
                                    tol = 1e-12) {
  h <- as.integer(high_state_size)
  np <- as.integer(n_positive)
  nn <- as.integer(n_negative)

  if (np < 1L || nn < 1L) return(NA_real_)
  if (h < 0L || h > np + nn) return(NA_real_)
  if (!is.finite(observed_J)) return(NA_real_)

  k_min <- max(0L, h - nn)
  k_max <- min(h, np)
  k <- seq.int(k_min, k_max)

  tpr_high <- k / np
  fpr_high <- (h - k) / nn
  null_J <- abs(tpr_high - fpr_high)

  probs <- stats::dhyper(k, m = np, n = nn, k = h)
  p <- sum(probs[null_J >= observed_J - tol])
  min(1, max(0, as.numeric(p)))
}


# ----------------------------------------------------------------------------
# Evaluate a locked threshold and locked phenotype orientation on new samples
# ----------------------------------------------------------------------------
evaluate_locked_rule_v11 <- function(x,
                                     labels,
                                     threshold,
                                     positive_state,
                                     positive_label) {
  labels_chr <- as.character(labels)
  y <- labels_chr == as.character(positive_label)

  if (length(y) == 0L || sum(y) == 0L || sum(!y) == 0L) {
    return(list(
      J = NA_real_, sensitivity = NA_real_, specificity = NA_real_,
      balanced_accuracy = NA_real_, accuracy = NA_real_,
      n_positive = sum(y), n_negative = sum(!y)
    ))
  }

  high_state <- x > threshold
  predicted_positive <- if (positive_state == "high") high_state else !high_state

  sensitivity <- sum(predicted_positive & y) / sum(y)
  specificity <- sum(!predicted_positive & !y) / sum(!y)

  list(
    J = as.numeric(sensitivity + specificity - 1),
    sensitivity = as.numeric(sensitivity),
    specificity = as.numeric(specificity),
    balanced_accuracy = as.numeric((sensitivity + specificity) / 2),
    accuracy = as.numeric(mean(predicted_positive == y)),
    n_positive = sum(y),
    n_negative = sum(!y)
  )
}


# ----------------------------------------------------------------------------
# Orientation-free rank AUC benchmark
# ----------------------------------------------------------------------------
# This is max(AUC, 1-AUC), so it is deliberately orientation-free and is named
# accordingly in outputs to avoid confusion with a conventional directional AUC.
# ----------------------------------------------------------------------------
orientation_free_auc_v11 <- function(x, labels, positive_label) {
  y <- as.character(labels) == as.character(positive_label)
  n_pos <- sum(y)
  n_neg <- sum(!y)
  if (n_pos == 0L || n_neg == 0L) return(NA_real_)

  r <- rank(x, ties.method = "average")
  auc_raw <- (sum(r[y]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
  max(auc_raw, 1 - auc_raw)
}


# ----------------------------------------------------------------------------
# One stratified training/held-out resample
# ----------------------------------------------------------------------------
make_stratified_resample_v11 <- function(labels,
                                         method = c("subsample", "bootstrap"),
                                         fraction = 0.80) {
  method <- match.arg(method)
  labels <- as.factor(labels)
  groups <- split(seq_along(labels), labels)

  if (method == "subsample") {
    train <- unlist(lapply(groups, function(idx) {
      ng <- length(idx)
      if (ng == 1L) return(idx)

      # Keep at least one held-out observation per class whenever possible.
      n_train <- floor(fraction * ng)
      n_train <- max(1L, n_train)
      if (ng >= 3L) n_train <- max(2L, n_train)
      n_train <- min(ng - 1L, n_train)

      sample(idx, size = n_train, replace = FALSE)
    }), use.names = FALSE)

    test <- setdiff(seq_along(labels), train)
    return(list(train = train, heldout = test))
  }

  # Stratified bootstrap: each class is sampled at its original class size.
  # The held-out set is the out-of-bag set; some iterations may lack one class
  # and are then ineligible for held-out J calculation.
  train <- unlist(lapply(groups, function(idx) {
    sample(idx, size = length(idx), replace = TRUE)
  }), use.names = FALSE)

  heldout <- setdiff(seq_along(labels), unique(train))
  list(train = train, heldout = heldout)
}


# ----------------------------------------------------------------------------
# Internal stability assessment for one feature
# ----------------------------------------------------------------------------
# IMPORTANT: this is intentionally called stability assessment, not validation.
# When only top-ranked full-sample features are resampled, the procedure is
# post-selection and is not an unbiased estimate of an end-to-end discovery
# pipeline's predictive performance.
# ----------------------------------------------------------------------------
resample_feature_stability_v11 <- function(x,
                                           labels,
                                           positive_label,
                                           full_fit,
                                           full_positive_state,
                                           n_resamples = 200L,
                                           method = c("subsample", "bootstrap"),
                                           fraction = 0.80,
                                           stability_seed) {
  method <- match.arg(method)
  stability_seed <- .validate_seed_v11(stability_seed, "stability_seed")
  withr::local_seed(stability_seed)
  n <- length(x)

  thresholds <- rep(NA_real_, n_resamples)
  training_J <- rep(NA_real_, n_resamples)
  heldout_J <- rep(NA_real_, n_resamples)
  training_sensitivity <- rep(NA_real_, n_resamples)
  training_specificity <- rep(NA_real_, n_resamples)
  heldout_sensitivity <- rep(NA_real_, n_resamples)
  heldout_specificity <- rep(NA_real_, n_resamples)
  positive_state <- rep(NA_character_, n_resamples)
  state_partition_agreement <- rep(NA_real_, n_resamples)
  heldout_n_positive <- rep(NA_real_, n_resamples)
  heldout_n_negative <- rep(NA_real_, n_resamples)

  high_state_counts <- numeric(n)
  predicted_positive_counts <- numeric(n)
  valid_runs <- 0L

  full_high_state <- full_fit$high_state

  for (b in seq_len(n_resamples)) {
    rs <- make_stratified_resample_v11(
      labels,
      method = method,
      fraction = fraction
    )
    train <- rs$train
    heldout <- rs$heldout

    fit <- exact_two_medoids_1d_v11(x[train])
    if (!isTRUE(fit$valid)) next

    sc_train <- score_fixed_partition_v11(
      fit$high_state,
      labels[train],
      positive_label
    )

    thresholds[b] <- fit$threshold
    training_J[b] <- sc_train$J
    training_sensitivity[b] <- sc_train$sensitivity
    training_specificity[b] <- sc_train$specificity
    positive_state[b] <- sc_train$positive_state

    # Apply the training-derived threshold to all original observations only to
    # measure state-assignment stability relative to the full-sample partition.
    high_all <- x > fit$threshold
    state_partition_agreement[b] <- mean(high_all == full_high_state)

    high_state_counts <- high_state_counts + as.numeric(high_all)
    pred_all <- if (sc_train$positive_state == "high") high_all else !high_all
    predicted_positive_counts <- predicted_positive_counts + as.numeric(pred_all)
    valid_runs <- valid_runs + 1L

    # Held-out/OOB assessment: both threshold and phenotype orientation are
    # locked using training data before held-out samples are evaluated.
    if (length(heldout) > 0L && length(unique(as.character(labels[heldout]))) == 2L) {
      ev <- evaluate_locked_rule_v11(
        x = x[heldout],
        labels = labels[heldout],
        threshold = fit$threshold,
        positive_state = sc_train$positive_state,
        positive_label = positive_label
      )

      heldout_J[b] <- ev$J
      heldout_sensitivity[b] <- ev$sensitivity
      heldout_specificity[b] <- ev$specificity
      heldout_n_positive[b] <- ev$n_positive
      heldout_n_negative[b] <- ev$n_negative
    }
  }

  valid_fit <- is.finite(thresholds)
  valid_heldout <- is.finite(heldout_J)
  iqr_x <- stats::IQR(x)

  threshold_iqr <- if (any(valid_fit)) stats::IQR(thresholds[valid_fit]) else NA_real_
  threshold_iqr_scaled <- if (
    is.finite(threshold_iqr) && is.finite(iqr_x) && iqr_x > 0
  ) threshold_iqr / iqr_x else NA_real_

  list(
    stability_seed = stability_seed,
    n_resamples_requested = n_resamples,
    n_resamples_valid = sum(valid_fit),
    n_heldout_valid = sum(valid_heldout),
    threshold_median = if (any(valid_fit)) stats::median(thresholds[valid_fit]) else NA_real_,
    threshold_q025 = .quantile_finite_v11(thresholds, 0.025),
    threshold_q975 = .quantile_finite_v11(thresholds, 0.975),
    threshold_iqr = threshold_iqr,
    threshold_iqr_scaled = threshold_iqr_scaled,
    median_training_J = if (any(is.finite(training_J))) {
      stats::median(training_J, na.rm = TRUE)
    } else NA_real_,
    median_heldout_J = if (any(valid_heldout)) {
      stats::median(heldout_J[valid_heldout])
    } else NA_real_,
    heldout_J_q25 = .quantile_finite_v11(heldout_J, 0.25),
    heldout_J_q75 = .quantile_finite_v11(heldout_J, 0.75),
    heldout_J_q025 = .quantile_finite_v11(heldout_J, 0.025),
    heldout_J_q975 = .quantile_finite_v11(heldout_J, 0.975),
    positive_state_stability = if (any(valid_fit)) {
      mean(positive_state[valid_fit] == full_positive_state)
    } else NA_real_,
    state_partition_stability = if (any(is.finite(state_partition_agreement))) {
      mean(state_partition_agreement, na.rm = TRUE)
    } else NA_real_,
    exact_training_fraction = if (any(is.finite(training_J))) {
      mean(training_J[is.finite(training_J)] >= 1 - 1e-12)
    } else NA_real_,
    exact_heldout_fraction = if (any(valid_heldout)) {
      mean(heldout_J[valid_heldout] >= 1 - 1e-12)
    } else NA_real_,
    median_heldout_n_positive = .quantile_finite_v11(heldout_n_positive, 0.5),
    median_heldout_n_negative = .quantile_finite_v11(heldout_n_negative, 0.5),
    high_state_frequency = if (valid_runs > 0L) {
      high_state_counts / valid_runs
    } else rep(NA_real_, n),
    predicted_positive_frequency = if (valid_runs > 0L) {
      predicted_positive_counts / valid_runs
    } else rep(NA_real_, n),
    thresholds = thresholds,
    training_J = training_J,
    heldout_J = heldout_J,
    positive_state = positive_state,
    state_partition_agreement = state_partition_agreement,
    heldout_n_positive = heldout_n_positive,
    heldout_n_negative = heldout_n_negative
  )
}


# ----------------------------------------------------------------------------
# Main analysis
# ----------------------------------------------------------------------------
label_blind_state_concordance_v11 <- function(
    feature_matrix,
    labels,
    positive_label,
    n_resamples = 200L,
    resample_method = c("subsample", "bootstrap", "none"),
    resample_fraction = 0.80,
    resample_top_n = 500L,
    seed = 1L,
    keep_high_state_matrix = TRUE,
    verbose = TRUE) {

  validated <- .validate_binary_inputs_v11(
    feature_matrix = feature_matrix,
    labels = labels,
    positive_label = positive_label
  )

  feature_matrix <- validated$feature_matrix
  labels <- validated$labels
  positive_label <- validated$positive_label
  negative_label <- validated$negative_label
  n_positive <- validated$n_positive
  n_negative <- validated$n_negative

  resample_method <- match.arg(resample_method)

  if (!is.numeric(n_resamples) || length(n_resamples) != 1L ||
      !is.finite(n_resamples) || n_resamples < 0 || n_resamples != floor(n_resamples)) {
    stop("n_resamples must be a non-negative integer")
  }
  n_resamples <- as.integer(n_resamples)

  if (!is.numeric(resample_fraction) || length(resample_fraction) != 1L ||
      !is.finite(resample_fraction) || resample_fraction <= 0 || resample_fraction >= 1) {
    stop("resample_fraction must be strictly between 0 and 1")
  }

  if (!is.null(resample_top_n) && !is.infinite(resample_top_n)) {
    if (!is.numeric(resample_top_n) || length(resample_top_n) != 1L ||
        !is.finite(resample_top_n) || resample_top_n < 1) {
      stop("resample_top_n must be NULL, Inf, or a positive integer")
    }
    resample_top_n <- as.integer(resample_top_n)
  }

  seed <- .validate_seed_v11(seed, "seed")

  feature_ids <- rownames(feature_matrix)
  G <- nrow(feature_matrix)
  N <- ncol(feature_matrix)

  if (verbose) {
    message(
      "v11: fitting label-blind 1D two-state partitions for ", G,
      " features across ", N, " samples (", n_positive, " ", positive_label,
      " vs ", n_negative, " ", negative_label, ")."
    )
  }

  out <- vector("list", G)
  full_details <- vector("list", G)

  high_state_matrix <- if (isTRUE(keep_high_state_matrix)) {
    matrix(
      FALSE,
      nrow = G,
      ncol = N,
      dimnames = list(feature_ids, colnames(feature_matrix))
    )
  } else NULL

  for (g in seq_len(G)) {
    feature_id <- feature_ids[g]
    x <- feature_matrix[g, ]

    fit <- exact_two_medoids_1d_v11(x)

    if (!isTRUE(fit$valid)) {
      out[[g]] <- data.frame(
        feature_id = feature_id,
        valid = FALSE,
        reason = fit$reason,
        threshold = NA_real_,
        low_medoid = NA_real_,
        high_medoid = NA_real_,
        medoid_distance = NA_real_,
        gap = NA_real_,
        relative_gap = NA_real_,
        within_l1_cost = NA_real_,
        high_state_size = NA_integer_,
        J = 0,
        balanced_accuracy = 0.5,
        sensitivity = NA_real_,
        specificity = NA_real_,
        accuracy = NA_real_,
        positive_state = NA_character_,
        exact_full_concordance = FALSE,
        exact_J_p = 1,
        fisher_p = 1,
        orientation_free_auc = 0.5,
        wilcox_p = 1,
        median_difference_positive_minus_negative = 0,
        stringsAsFactors = FALSE
      )
      next
    }

    sc <- score_fixed_partition_v11(
      high_state = fit$high_state,
      labels = labels,
      positive_label = positive_label
    )

    y <- as.character(labels) == positive_label
    high_state_size <- sum(fit$high_state)

    exact_J_p <- exact_J_hypergeom_p_v11(
      high_state_size = high_state_size,
      n_positive = n_positive,
      n_negative = n_negative,
      observed_J = sc$J
    )

    auc <- orientation_free_auc_v11(x, labels, positive_label)

    wilcox_p <- tryCatch(
      suppressWarnings(
        stats::wilcox.test(x[y], x[!y], exact = FALSE)$p.value
      ),
      error = function(e) NA_real_
    )

    median_difference <- stats::median(x[y]) - stats::median(x[!y])
    exact_full <- (
      sc$J >= 1 - 1e-12 &&
      sc$sensitivity >= 1 - 1e-12 &&
      sc$specificity >= 1 - 1e-12
    )

    if (isTRUE(keep_high_state_matrix)) {
      high_state_matrix[g, ] <- fit$high_state
    }

    # Store by integer row position. This avoids name-based accidental reuse and
    # remains correct even if future code changes the displayed feature labels.
    full_details[[g]] <- list(
      threshold = fit$threshold,
      high_state = fit$high_state,
      low_medoid = fit$low_medoid,
      high_medoid = fit$high_medoid,
      gap = fit$gap,
      score = sc
    )

    out[[g]] <- data.frame(
      feature_id = feature_id,
      valid = TRUE,
      reason = NA_character_,
      threshold = fit$threshold,
      low_medoid = fit$low_medoid,
      high_medoid = fit$high_medoid,
      medoid_distance = fit$medoid_distance,
      gap = fit$gap,
      relative_gap = fit$relative_gap,
      within_l1_cost = fit$within_l1_cost,
      high_state_size = high_state_size,
      J = sc$J,
      balanced_accuracy = sc$balanced_accuracy,
      sensitivity = sc$sensitivity,
      specificity = sc$specificity,
      accuracy = sc$accuracy,
      positive_state = sc$positive_state,
      exact_full_concordance = exact_full,
      exact_J_p = exact_J_p,
      fisher_p = sc$fisher_p,
      orientation_free_auc = auc,
      wilcox_p = wilcox_p,
      median_difference_positive_minus_negative = median_difference,
      stringsAsFactors = FALSE
    )

    if (verbose && (g %% 2000L == 0L || g == G)) {
      message("  label-blind partitioning: ", g, "/", G, " features")
    }
  }

  summary_table <- do.call(rbind, out)

  # Feature-wise statistical evidence. exact_J_p is the primary inferential
  # measure aligned to J. Fisher is retained only as a familiar 2x2 benchmark
  # (and may be numerically identical to exact_J_p for many observed tables);
  # Wilcoxon is retained as a distributional benchmark.
  summary_table$exact_J_fdr <- stats::p.adjust(summary_table$exact_J_p, method = "BH")
  summary_table$fisher_fdr <- stats::p.adjust(summary_table$fisher_p, method = "BH")
  summary_table$wilcox_fdr <- stats::p.adjust(summary_table$wilcox_p, method = "BH")

  # Primary ranking is the interpretable concordance effect size J. Exact-J FDR
  # and the label-blind relative state gap are tie breakers, not a composite score.
  ord <- order(
    -summary_table$J,
    summary_table$exact_J_fdr,
    -summary_table$relative_gap,
    na.last = TRUE
  )
  summary_table <- summary_table[ord, , drop = FALSE]
  rownames(summary_table) <- NULL
  summary_table$rank_J <- seq_len(nrow(summary_table))

  # Stability columns remain NA for features not selected for resampling.
  stability_cols <- c(
    "n_resamples_valid", "n_heldout_valid",
    "threshold_median", "threshold_q025", "threshold_q975",
    "threshold_iqr", "threshold_iqr_scaled",
    "median_training_J", "median_heldout_J",
    "heldout_J_q25", "heldout_J_q75", "heldout_J_q025", "heldout_J_q975",
    "positive_state_stability", "state_partition_stability",
    "exact_training_fraction", "exact_heldout_fraction",
    "median_heldout_n_positive", "median_heldout_n_negative"
  )
  for (nm in stability_cols) summary_table[[nm]] <- NA_real_
  # The exact feature-specific RNG seed is exposed for audit/reproduction.
  summary_table$stability_seed <- NA_integer_

  resampling_details <- list()

  if (resample_method != "none" && n_resamples > 0L) {
    valid_rows <- which(summary_table$valid)
    if (is.null(resample_top_n) || is.infinite(resample_top_n)) {
      chosen_rows <- valid_rows
    } else {
      chosen_rows <- head(valid_rows, min(resample_top_n, length(valid_rows)))
    }

    if (verbose) {
      message(
        "v11: internal ", resample_method, " stability assessment for ",
        length(chosen_rows), " top-ranked features x ", n_resamples,
        " resamples. This is not external validation."
      )
    }

    for (ii in seq_along(chosen_rows)) {
      sr <- chosen_rows[ii]
      feature_id <- summary_table$feature_id[sr]
      original_idx <- match(feature_id, feature_ids)
      if (is.na(original_idx)) stop("internal feature-index mapping failure")

      x <- feature_matrix[original_idx, ]
      fd <- full_details[[original_idx]]
      stability_seed <- .feature_stability_seed_v11(seed, feature_id)

      rr <- resample_feature_stability_v11(
        x = x,
        labels = labels,
        positive_label = positive_label,
        full_fit = list(
          high_state = fd$high_state,
          threshold = fd$threshold
        ),
        full_positive_state = summary_table$positive_state[sr],
        n_resamples = n_resamples,
        method = resample_method,
        fraction = resample_fraction,
        stability_seed = stability_seed
      )

      summary_table$stability_seed[sr] <- stability_seed
      summary_table$n_resamples_valid[sr] <- rr$n_resamples_valid
      summary_table$n_heldout_valid[sr] <- rr$n_heldout_valid
      summary_table$threshold_median[sr] <- rr$threshold_median
      summary_table$threshold_q025[sr] <- rr$threshold_q025
      summary_table$threshold_q975[sr] <- rr$threshold_q975
      summary_table$threshold_iqr[sr] <- rr$threshold_iqr
      summary_table$threshold_iqr_scaled[sr] <- rr$threshold_iqr_scaled
      summary_table$median_training_J[sr] <- rr$median_training_J
      summary_table$median_heldout_J[sr] <- rr$median_heldout_J
      summary_table$heldout_J_q25[sr] <- rr$heldout_J_q25
      summary_table$heldout_J_q75[sr] <- rr$heldout_J_q75
      summary_table$heldout_J_q025[sr] <- rr$heldout_J_q025
      summary_table$heldout_J_q975[sr] <- rr$heldout_J_q975
      summary_table$positive_state_stability[sr] <- rr$positive_state_stability
      summary_table$state_partition_stability[sr] <- rr$state_partition_stability
      summary_table$exact_training_fraction[sr] <- rr$exact_training_fraction
      summary_table$exact_heldout_fraction[sr] <- rr$exact_heldout_fraction
      summary_table$median_heldout_n_positive[sr] <- rr$median_heldout_n_positive
      summary_table$median_heldout_n_negative[sr] <- rr$median_heldout_n_negative

      rr$sample_history <- data.frame(
        sample_id = colnames(feature_matrix),
        label = as.character(labels),
        full_state = ifelse(fd$high_state, "high", "low"),
        high_state_frequency = rr$high_state_frequency,
        predicted_positive_frequency = rr$predicted_positive_frequency,
        stringsAsFactors = FALSE
      )

      resampling_details[[feature_id]] <- rr

      if (verbose && (ii %% 50L == 0L || ii == length(chosen_rows))) {
        message("  stability assessment: ", ii, "/", length(chosen_rows), " features")
      }
    }
  }

  result <- list(
    method = paste(
      "Label-blind single-feature 1D two-medoids state concordance v11:",
      "fixed-partition J + exact conditional J inference + optional internal stability"
    ),
    positive_label = positive_label,
    negative_label = negative_label,
    labels = labels,
    sample_ids = colnames(feature_matrix),
    n_positive = n_positive,
    n_negative = n_negative,
    parameters = list(
      n_resamples = n_resamples,
      resample_method = resample_method,
      resample_fraction = resample_fraction,
      resample_top_n = resample_top_n,
      seed = seed,
      stability_rng_scheme = "deterministic feature-ID seed derived from master seed"
    ),
    summary = summary_table,
    resampling_details = resampling_details,
    high_state_matrix = high_state_matrix
  )

  class(result) <- c("lbsc_v11", "list")
  result
}


# ----------------------------------------------------------------------------
# Genome-wide permutation calibration for fixed full-sample partitions
# ----------------------------------------------------------------------------
# This procedure preserves each feature's label-blind high/low state assignment
# and permutes phenotype labels. It estimates:
#   - maxT FWER-adjusted p-values for each observed J;
#   - the null distribution of the genome-wide maximum J;
#   - null distributions for the number of features above chosen J cutoffs.
#
# It does NOT rerun feature partitioning (unnecessary because partitioning is
# label-blind) and does NOT rerun resampling or feature selection. Therefore it
# calibrates the fixed-partition full-sample concordance statistics, not the
# complete post-selection/resampling workflow.
# ----------------------------------------------------------------------------
permutation_global_calibration_v11 <- function(v11_result,
                                               n_permutations = 1000L,
                                               cutoffs = c(0.80, 1.00),
                                               seed = 1L,
                                               verbose = TRUE) {
  if (!inherits(v11_result, "lbsc_v11")) {
    stop("v11_result must be an lbsc_v11 object")
  }
  if (is.null(v11_result$high_state_matrix)) {
    stop("v11_result must have keep_high_state_matrix=TRUE")
  }
  if (!is.numeric(n_permutations) || length(n_permutations) != 1L ||
      !is.finite(n_permutations) || n_permutations < 1 ||
      n_permutations != floor(n_permutations)) {
    stop("n_permutations must be a positive integer")
  }
  n_permutations <- as.integer(n_permutations)

  if (!is.numeric(cutoffs) || any(!is.finite(cutoffs)) ||
      any(cutoffs < 0 | cutoffs > 1)) {
    stop("cutoffs must be finite values between 0 and 1")
  }
  cutoffs <- sort(unique(as.numeric(cutoffs)))
  seed <- .validate_seed_v11(seed, "seed")

  summary_table <- v11_result$summary
  matrix_feature_ids <- rownames(v11_result$high_state_matrix)
  valid_feature <- summary_table$valid[match(matrix_feature_ids, summary_table$feature_id)]
  valid_feature[is.na(valid_feature)] <- FALSE

  M <- v11_result$high_state_matrix[valid_feature, , drop = FALSE]
  feature_ids <- rownames(M)
  observed_J <- summary_table$J[match(feature_ids, summary_table$feature_id)]

  labels <- v11_result$labels
  y <- as.character(labels) == v11_result$positive_label
  n_pos <- sum(y)
  n_neg <- sum(!y)

  G <- nrow(M)
  high_size <- rowSums(M)
  max_J <- numeric(n_permutations)
  null_counts <- matrix(0L, nrow = n_permutations, ncol = length(cutoffs))
  colnames(null_counts) <- paste0("J_ge_", format(cutoffs, trim = TRUE))
  observed_counts <- vapply(cutoffs, function(z) sum(observed_J >= z - 1e-12), integer(1))
  names(observed_counts) <- colnames(null_counts)

  withr::local_seed(seed)
  Mnum <- M * 1

  if (verbose) {
    message(
      "v11 fixed-partition permutation calibration: ", G, " valid features x ",
      n_permutations, " label permutations."
    )
  }

  for (b in seq_len(n_permutations)) {
    yp <- sample(y, replace = FALSE)

    positive_high <- as.vector(Mnum %*% as.numeric(yp))
    negative_high <- high_size - positive_high
    Jp <- abs(positive_high / n_pos - negative_high / n_neg)

    max_J[b] <- max(Jp)
    for (j in seq_along(cutoffs)) {
      null_counts[b, j] <- sum(Jp >= cutoffs[j] - 1e-12)
    }

    if (verbose && (b %% 100L == 0L || b == n_permutations)) {
      message("  permutations: ", b, "/", n_permutations)
    }
  }

  maxT_fwer_p <- vapply(observed_J, function(z) {
    (sum(max_J >= z - 1e-12) + 1) / (n_permutations + 1)
  }, numeric(1))

  count_global_p <- vapply(seq_along(cutoffs), function(j) {
    (sum(null_counts[, j] >= observed_counts[j]) + 1) / (n_permutations + 1)
  }, numeric(1))
  names(count_global_p) <- colnames(null_counts)

  per_feature <- data.frame(
    feature_id = feature_ids,
    J = observed_J,
    maxT_fwer_p = maxT_fwer_p,
    stringsAsFactors = FALSE
  )

  list(
    per_feature = per_feature,
    null_max_J = max_J,
    null_counts = null_counts,
    observed_counts = observed_counts,
    count_global_p = count_global_p,
    cutoffs = cutoffs,
    n_permutations = n_permutations,
    monte_carlo_min_p = 1 / (n_permutations + 1),
    seed = seed,
    interpretation = paste(
      "Label permutations calibrate full-sample phenotype concordance of the",
      "fixed label-blind state partitions; they do not rerun internal stability resampling."
    )
  )
}


# ----------------------------------------------------------------------------
# Attach genome-wide permutation results to the main object
# ----------------------------------------------------------------------------
add_permutation_calibration_v11 <- function(v11_result, permutation_result) {
  p <- permutation_result$per_feature
  idx <- match(v11_result$summary$feature_id, p$feature_id)

  v11_result$summary$maxT_fwer_p <- p$maxT_fwer_p[idx]
  v11_result$permutation <- permutation_result
  v11_result
}


# ----------------------------------------------------------------------------
# Publication-oriented filtering without collapsing separate evidence axes
# ----------------------------------------------------------------------------
# This helper only filters/ranks. It deliberately does not create a label such
# as "robust discriminator" that could be mistaken for inferential validation.
# Users can choose thresholds appropriate to their study and report the actual
# numerical effect, stability and evidence metrics.
# ----------------------------------------------------------------------------
top_features_v11 <- function(v11_result,
                             n = 50L,
                             min_J = 0.80,
                             max_exact_J_fdr = 0.05,
                             require_stability = FALSE,
                             min_median_heldout_J = NULL,
                             min_positive_state_stability = NULL,
                             min_state_partition_stability = NULL,
                             min_exact_heldout_fraction = NULL,
                             max_maxT_fwer_p = NULL) {
  if (!inherits(v11_result, "lbsc_v11")) {
    stop("v11_result must be an lbsc_v11 object")
  }

  d <- v11_result$summary
  keep <- d$valid & is.finite(d$J) & d$J >= min_J &
    is.finite(d$exact_J_fdr) & d$exact_J_fdr <= max_exact_J_fdr

  if (require_stability) {
    keep <- keep & is.finite(d$median_heldout_J)
  }
  if (!is.null(min_median_heldout_J)) {
    keep <- keep & is.finite(d$median_heldout_J) &
      d$median_heldout_J >= min_median_heldout_J
  }
  if (!is.null(min_positive_state_stability)) {
    keep <- keep & is.finite(d$positive_state_stability) &
      d$positive_state_stability >= min_positive_state_stability
  }
  if (!is.null(min_state_partition_stability)) {
    keep <- keep & is.finite(d$state_partition_stability) &
      d$state_partition_stability >= min_state_partition_stability
  }
  if (!is.null(min_exact_heldout_fraction)) {
    keep <- keep & is.finite(d$exact_heldout_fraction) &
      d$exact_heldout_fraction >= min_exact_heldout_fraction
  }
  if (!is.null(max_maxT_fwer_p)) {
    if (!"maxT_fwer_p" %in% names(d)) {
      stop("maxT_fwer_p is unavailable; run and attach permutation calibration first")
    }
    keep <- keep & is.finite(d$maxT_fwer_p) & d$maxT_fwer_p <= max_maxT_fwer_p
  }

  ans <- d[keep, , drop = FALSE]
  ans <- ans[order(
    -ans$J,
    ans$exact_J_fdr,
    -ans$median_heldout_J,
    na.last = TRUE
  ), , drop = FALSE]

  head(ans, as.integer(n))
}


# ----------------------------------------------------------------------------
# Align a metadata field to feature-matrix columns
# ----------------------------------------------------------------------------
build_labels_v11 <- function(metadata, group_column, feature_matrix) {
  if (is.null(rownames(metadata))) {
    stop("metadata must have row names corresponding to sample IDs")
  }
  if (anyDuplicated(rownames(metadata))) {
    stop("metadata row names/sample IDs must be unique")
  }
  if (!group_column %in% colnames(metadata)) {
    stop("group_column not found in metadata")
  }
  if (is.null(colnames(feature_matrix))) {
    stop("feature_matrix must have sample column names")
  }

  idx <- match(colnames(feature_matrix), rownames(metadata))
  if (anyNA(idx)) {
    stop("some feature-matrix sample IDs are absent from metadata")
  }

  metadata[[group_column]][idx]
}


# ----------------------------------------------------------------------------
# Compact print method
# ----------------------------------------------------------------------------
print.lbsc_v11 <- function(x, ...) {
  cat("Label-blind single-feature state concordance v11\n")
  cat("  samples: ", length(x$sample_ids), "\n", sep = "")
  cat("  positive: ", x$n_positive, " [", x$positive_label, "]\n", sep = "")
  cat("  negative: ", x$n_negative, " [", x$negative_label, "]\n", sep = "")
  cat("  features: ", nrow(x$summary), "\n", sep = "")
  cat("  exact full-sample concordance: ",
      sum(x$summary$exact_full_concordance, na.rm = TRUE), "\n", sep = "")
  cat("  exact-J FDR <= 0.05: ",
      sum(x$summary$exact_J_fdr <= 0.05, na.rm = TRUE), "\n", sep = "")
  cat("  stability assessed: ",
      sum(is.finite(x$summary$median_heldout_J)), "\n", sep = "")
  if ("maxT_fwer_p" %in% names(x$summary)) {
    cat("  maxT FWER <= 0.05: ",
        sum(x$summary$maxT_fwer_p <= 0.05, na.rm = TRUE), "\n", sep = "")
  }
  invisible(x)
}


# ============================================================================
# Minimal usage example
# ============================================================================
# labels <- build_labels_v11(metadata, "group", feature_matrix)
#
# fit <- label_blind_state_concordance_v11(
#   feature_matrix = feature_matrix,
#   labels = labels,
#   positive_label = "case",            # explicit by design
#   n_resamples = 500,
#   resample_method = "subsample",
#   resample_fraction = 0.80,
#   resample_top_n = 1000,
#   seed = 123,
#   keep_high_state_matrix = TRUE
# )
#
# # Optional genome-wide calibration of the FIXED full-sample partitions.
# # 1000 permutations gives a minimum Monte Carlo p of 1/1001 ~= 0.001;
# # increase the number for finer maxT resolution if computationally feasible.
# perm <- permutation_global_calibration_v11(
#   fit,
#   n_permutations = 1000,
#   cutoffs = c(0.80, 1.00),
#   seed = 456
# )
# fit <- add_permutation_calibration_v11(fit, perm)
#
# head(fit$summary, 30)
# top_features_v11(
#   fit,
#   n = 30,
#   min_J = 0.80,
#   max_exact_J_fdr = 0.05,
#   require_stability = TRUE,
#   min_median_heldout_J = 0.60,
#   min_positive_state_stability = 0.90
# )
# ============================================================================
