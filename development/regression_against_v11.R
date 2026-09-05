# Run from the IntrinsicScan package root after installing/loading the package.
# This checks that the Bioconductor wrapper preserves the frozen v11 numerical engine.

if (!requireNamespace("IntrinsicScan", quietly = TRUE)) {
    stop("Install IntrinsicScan first, e.g. R CMD INSTALL .")
}
source(file.path("development", "label_blind_state_concordance_v11_reference.R"))

set.seed(111)
mat <- matrix(rnorm(30 * 20), nrow = 30,
              dimnames = list(paste0("F", 1:30), paste0("S", 1:20)))
# Inject several clear partitions and a constant feature.
mat[1, 11:20] <- mat[1, 11:20] + 4
mat[2, 11:20] <- mat[2, 11:20] + 2
mat[30, ] <- 7
labs <- rep(c("control", "case"), each = 10)

legacy <- label_blind_state_concordance_v11(
    mat, labs, positive_label = "case",
    n_resamples = 25, resample_method = "subsample",
    resample_fraction = 0.8, resample_top_n = 10,
    seed = 123, keep_high_state_matrix = TRUE, verbose = FALSE
)

modern <- IntrinsicScan::intrinsicScan(
    mat, labs, positive = "case",
    nResamples = 25, resampleMethod = "subsample",
    resampleFraction = 0.8, resampleTop = 10,
    seed = 123, keepStates = TRUE, verbose = FALSE
)
md <- as.data.frame(IntrinsicScan::intrinsicResults(modern))
ld <- legacy$summary

cols <- c("feature_id", "valid", "threshold", "low_medoid", "high_medoid",
          "medoid_distance", "gap", "relative_gap", "within_l1_cost",
          "high_state_size", "J", "balanced_accuracy", "sensitivity",
          "specificity", "accuracy", "positive_state", "exact_full_concordance",
          "exact_J_p", "exact_J_fdr", "rank_J", "stability_seed",
          "median_training_J", "median_heldout_J", "heldout_J_q25",
          "positive_state_stability", "state_partition_stability")
stopifnot(identical(ld$feature_id, md$feature_id))
for (nm in cols) {
    if (is.numeric(ld[[nm]]) || is.integer(ld[[nm]])) {
        stopifnot(isTRUE(all.equal(ld[[nm]], md[[nm]], tolerance = 0, check.attributes = FALSE)))
    } else {
        stopifnot(identical(ld[[nm]], md[[nm]]))
    }
}
state_modern <- SummarizedExperiment::assay(modern, "state")
state_modern <- state_modern[rownames(legacy$high_state_matrix), , drop = FALSE]
stopifnot(identical(unname(legacy$high_state_matrix), unname(state_modern)))

cat("IntrinsicScan wrapper is numerically identical to frozen v11 on the regression dataset.\n")
