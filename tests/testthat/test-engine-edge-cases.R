legacy_toy_fit <- function(n_resamples = 0L,
                           resample_method = "none",
                           keep_states = TRUE,
                           seed = 31L) {
    x <- rbind(
        A = c(
            1, 1.1, 1.2, 1.3, 1.4,
            5, 5.1, 5.2, 5.3, 5.4
        ),
        B = c(
            2, 2.1, 2.2, 2.3, 2.4,
            4, 4.1, 4.2, 4.3, 4.4
        ),
        C = c(
            1, 2, 3, 4, 5,
            6, 7, 8, 9, 10
        ),
        D = rep(3, 10)
    )
    colnames(x) <- paste0("S", seq_len(ncol(x)))

    y <- rep(c("control", "case"), each = 5)

    IntrinsicScan:::label_blind_state_concordance_v11(
        feature_matrix = x,
        labels = y,
        positive_label = "case",
        n_resamples = n_resamples,
        resample_method = resample_method,
        resample_fraction = 0.8,
        resample_top_n = Inf,
        seed = seed,
        keep_high_state_matrix = keep_states,
        verbose = FALSE
    )
}


test_that("legacy filtering helpers preserve separate evidence axes", {
    fit <- legacy_toy_fit(
        n_resamples = 12L,
        resample_method = "subsample"
    )

    expect_s3_class(fit, "lbsc_v11")

    top <- IntrinsicScan:::top_features_v11(
        fit,
        n = 10L,
        min_J = 0,
        max_exact_J_fdr = 1
    )

    expect_true(is.data.frame(top))
    expect_true(all(top$valid))
    expect_true(all(top$J >= 0))

    stable <- IntrinsicScan:::top_features_v11(
        fit,
        n = 10L,
        min_J = 0,
        max_exact_J_fdr = 1,
        require_stability = TRUE,
        min_median_heldout_J = 0,
        min_positive_state_stability = 0,
        min_state_partition_stability = 0,
        min_exact_heldout_fraction = 0
    )

    expect_true(all(is.finite(stable$median_heldout_J)))

    expect_error(
        IntrinsicScan:::top_features_v11(
            fit,
            min_J = 0,
            max_exact_J_fdr = 1,
            max_maxT_fwer_p = 1
        ),
        "permutation calibration"
    )

    perm <- IntrinsicScan:::permutation_global_calibration_v11(
        fit,
        n_permutations = 10L,
        cutoffs = c(0.5, 1),
        seed = 7L,
        verbose = FALSE
    )

    fit_perm <- IntrinsicScan:::add_permutation_calibration_v11(
        fit,
        perm
    )

    expect_true("maxT_fwer_p" %in% names(fit_perm$summary))

    top_perm <- IntrinsicScan:::top_features_v11(
        fit_perm,
        n = 10L,
        min_J = 0,
        max_exact_J_fdr = 1,
        max_maxT_fwer_p = 1
    )

    expect_true(all(top_perm$maxT_fwer_p <= 1))
})


test_that("legacy print method reports analysis summaries", {
    fit <- legacy_toy_fit(
        n_resamples = 5L,
        resample_method = "subsample"
    )

    txt <- capture.output(
        IntrinsicScan:::print.lbsc_v11(fit)
    )

    expect_true(any(grepl(
        "Label-blind single-feature state concordance v11",
        txt,
        fixed = TRUE
    )))
    expect_true(any(grepl("samples:", txt, fixed = TRUE)))
    expect_true(any(grepl("features:", txt, fixed = TRUE)))
    expect_true(any(grepl("stability assessed:", txt, fixed = TRUE)))

    perm <- IntrinsicScan:::permutation_global_calibration_v11(
        fit,
        n_permutations = 5L,
        seed = 9L,
        verbose = FALSE
    )

    fit <- IntrinsicScan:::add_permutation_calibration_v11(
        fit,
        perm
    )

    txt_perm <- capture.output(
        IntrinsicScan:::print.lbsc_v11(fit)
    )

    expect_true(any(grepl("maxT FWER", txt_perm, fixed = TRUE)))
})


test_that("build_labels_v11 aligns metadata by sample ID", {
    x <- matrix(
        1:12,
        nrow = 2,
        dimnames = list(
            c("A", "B"),
            c("S3", "S1", "S2", "S6", "S4", "S5")
        )
    )

    md <- data.frame(
        group = c("case", "control", "control",
                  "case", "case", "control"),
        row.names = c("S3", "S1", "S2", "S6", "S4", "S5")
    )

    labs <- IntrinsicScan:::build_labels_v11(
        md,
        "group",
        x
    )

    expect_identical(
        as.character(labs),
        as.character(md[colnames(x), "group"])
    )

    md_no_names <- S4Vectors::DataFrame(
    group = md$group
)

expect_null(rownames(md_no_names))

expect_error(
    IntrinsicScan:::build_labels_v11(
        md_no_names,
        "group",
        x
    ),
    "row names"
)

md_dup <- S4Vectors::DataFrame(
    group = md$group,
    row.names = c(
        "S1", "S1", "S2", "S3", "S4", "S5"
    )
)

expect_error(
    IntrinsicScan:::build_labels_v11(
        md_dup,
        "group",
        x
    ),
    "unique"
)

    expect_error(
        IntrinsicScan:::build_labels_v11(
            md,
            "missing_group",
            x
        ),
        "group_column"
    )

    x_no_names <- x
    colnames(x_no_names) <- NULL

    expect_error(
        IntrinsicScan:::build_labels_v11(
            md,
            "group",
            x_no_names
        ),
        "sample column names"
    )

    x_missing <- x
    colnames(x_missing)[1] <- "ABSENT"

    expect_error(
        IntrinsicScan:::build_labels_v11(
            md,
            "group",
            x_missing
        ),
        "absent from metadata"
    )
})


test_that("two-medoid helper handles invalid features explicitly", {
    one <- IntrinsicScan:::exact_two_medoids_1d_v11(1)
    expect_false(one$valid)
    expect_identical(
        one$reason,
        "fewer than 2 observations"
    )

    nonfinite <- IntrinsicScan:::exact_two_medoids_1d_v11(
        c(1, Inf)
    )
    expect_false(nonfinite$valid)
    expect_identical(
        nonfinite$reason,
        "non-finite values"
    )

    constant <- IntrinsicScan:::exact_two_medoids_1d_v11(
        rep(2, 5)
    )
    expect_false(constant$valid)
    expect_identical(
        constant$reason,
        "zero-variance feature"
    )

    valid <- IntrinsicScan:::exact_two_medoids_1d_v11(
        c(1, 1.1, 1.2, 5, 5.1, 5.2)
    )
    expect_true(valid$valid)
    expect_equal(valid$threshold, 3.1)
})


test_that("fixed-partition scoring validates its inputs", {
    score <- IntrinsicScan:::score_fixed_partition_v11(
        high_state = c(FALSE, FALSE, TRUE, TRUE),
        labels = c("control", "control", "case", "case"),
        positive_label = "case"
    )

    expect_equal(score$J, 1)

    expect_error(
        IntrinsicScan:::score_fixed_partition_v11(
            c(FALSE, TRUE),
            c("control"),
            "case"
        ),
        "lengths differ"
    )

    expect_error(
        IntrinsicScan:::score_fixed_partition_v11(
            c(FALSE, TRUE, TRUE),
            c("A", "B", "C"),
            "A"
        ),
        "exactly two"
    )

    expect_error(
        IntrinsicScan:::score_fixed_partition_v11(
            c(FALSE, FALSE, TRUE, TRUE),
            c("control", "control", "case", "case"),
            "absent"
        ),
        "not found"
    )
})


test_that("exact-J and locked-rule helpers handle undefined cases", {
    expect_true(is.na(
        IntrinsicScan:::exact_J_hypergeom_p_v11(
            2, 0, 4, 0.5
        )
    ))

    expect_true(is.na(
        IntrinsicScan:::exact_J_hypergeom_p_v11(
            -1, 2, 2, 0.5
        )
    ))

    expect_true(is.na(
        IntrinsicScan:::exact_J_hypergeom_p_v11(
            2, 2, 2, NA_real_
        )
    ))

    locked <- IntrinsicScan:::evaluate_locked_rule_v11(
        x = c(1, 2, 3),
        labels = rep("control", 3),
        threshold = 2,
        positive_state = "high",
        positive_label = "case"
    )

    expect_true(is.na(locked$J))
    expect_identical(locked$n_positive, 0L)

    auc <- IntrinsicScan:::orientation_free_auc_v11(
        x = c(1, 2, 3),
        labels = rep("control", 3),
        positive_label = "case"
    )

    expect_true(is.na(auc))

    expect_true(is.na(
        IntrinsicScan:::.quantile_finite_v11(
            c(NA_real_, Inf, -Inf),
            0.5
        )
    ))
})


test_that("stratified resampling covers subsample and bootstrap paths", {
    set.seed(1)

    sub <- IntrinsicScan:::make_stratified_resample_v11(
        labels = c("A", "B", "B", "B"),
        method = "subsample",
        fraction = 0.8
    )

    expect_true(1L %in% sub$train)
    expect_true(length(sub$heldout) >= 1L)

    labels <- rep(c("control", "case"), each = 5)

    set.seed(2)

    boot <- IntrinsicScan:::make_stratified_resample_v11(
        labels = labels,
        method = "bootstrap",
        fraction = 0.8
    )

    expect_length(boot$train, 10L)
    expect_true(all(boot$train %in% seq_along(labels)))
    expect_true(all(boot$heldout %in% seq_along(labels)))

    fit <- legacy_toy_fit(
        n_resamples = 10L,
        resample_method = "bootstrap",
        seed = 44L
    )

    expect_true(any(
        is.finite(fit$summary$median_training_J)
    ))
})


test_that("engine input validation rejects malformed data", {
    valid_x <- rbind(
        A = seq_len(10),
        B = seq_len(10) + 1
    )
    colnames(valid_x) <- paste0("S", 1:10)

    y <- rep(c("control", "case"), each = 5)

    expect_error(
        IntrinsicScan:::.validate_binary_inputs_v11(
            list(A = 1),
            y,
            "case"
        ),
        "matrix or data.frame"
    )

    too_small <- matrix(
        1,
        nrow = 1,
        dimnames = list("A", "S1")
    )

    expect_error(
        IntrinsicScan:::.validate_binary_inputs_v11(
            too_small,
            "case",
            "case"
        ),
        "at least 1 feature and 2 samples"
    )

    no_rows <- valid_x
    rownames(no_rows) <- NULL

    expect_error(
        IntrinsicScan:::.validate_binary_inputs_v11(
            no_rows,
            y,
            "case"
        ),
        "feature row names"
    )

    no_cols <- valid_x
    colnames(no_cols) <- NULL

    expect_error(
        IntrinsicScan:::.validate_binary_inputs_v11(
            no_cols,
            y,
            "case"
        ),
        "sample column names"
    )

    dup_rows <- valid_x
    rownames(dup_rows) <- c("A", "A")

    expect_error(
        IntrinsicScan:::.validate_binary_inputs_v11(
            dup_rows,
            y,
            "case"
        ),
        "unique"
    )

    dup_cols <- valid_x
    colnames(dup_cols)[2] <- colnames(dup_cols)[1]

    expect_error(
        IntrinsicScan:::.validate_binary_inputs_v11(
            dup_cols,
            y,
            "case"
        ),
        "unique"
    )

    bad_value <- valid_x
    bad_value[1, 1] <- Inf

    expect_error(
        IntrinsicScan:::.validate_binary_inputs_v11(
            bad_value,
            y,
            "case"
        ),
        "complete finite"
    )

    expect_error(
        IntrinsicScan:::.validate_binary_inputs_v11(
            valid_x,
            y[-1],
            "case"
        ),
        "labels length"
    )

    y_na <- y
    y_na[1] <- NA

    expect_error(
        IntrinsicScan:::.validate_binary_inputs_v11(
            valid_x,
            y_na,
            "case"
        ),
        "must not contain NA"
    )

    expect_error(
        IntrinsicScan:::.validate_binary_inputs_v11(
            valid_x,
            rep(c("A", "B", "C", "A", "B"), 2),
            "A"
        ),
        "exactly two"
    )

    expect_error(
        IntrinsicScan:::.validate_binary_inputs_v11(
            valid_x,
            y,
            "absent"
        ),
        "not present"
    )
})


test_that("engine analysis validates resampling parameters", {
    x <- rbind(
        A = c(
            1, 1.1, 1.2, 1.3, 1.4,
            5, 5.1, 5.2, 5.3, 5.4
        )
    )
    colnames(x) <- paste0("S", 1:10)
    y <- rep(c("control", "case"), each = 5)

    expect_error(
        IntrinsicScan:::label_blind_state_concordance_v11(
            x, y, "case",
            n_resamples = -1,
            verbose = FALSE
        ),
        "non-negative integer"
    )

    expect_error(
        IntrinsicScan:::label_blind_state_concordance_v11(
            x, y, "case",
            n_resamples = 1,
            resample_fraction = 1,
            verbose = FALSE
        ),
        "strictly between 0 and 1"
    )

    expect_error(
        IntrinsicScan:::label_blind_state_concordance_v11(
            x, y, "case",
            n_resamples = 1,
            resample_top_n = 0,
            verbose = FALSE
        ),
        "positive integer"
    )

    expect_error(
        IntrinsicScan:::label_blind_state_concordance_v11(
            x, y, "case",
            seed = -1,
            verbose = FALSE
        ),
        "integer between"
    )
})


test_that("missing phenotype metadata column is rejected", {
    x <- matrix(
        seq_len(20),
        nrow = 2,
        dimnames = list(
            c("A", "B"),
            paste0("S", 1:10)
        )
    )

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(expression = x),
        colData = S4Vectors::DataFrame(
            group = rep(c("control", "case"), each = 5),
            row.names = colnames(x)
        )
    )

    expect_error(
        intrinsicScan(
            se,
            phenotype = "missing_column",
            positive = "case",
            verbose = FALSE
        ),
        "phenotype column not found"
    )
})