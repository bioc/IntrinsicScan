toy_data <- function() {
    x <- rbind(
        feature_A = c(
            1, 1.1, 1.2, 1.3, 1.4,
            5, 5.1, 5.2, 5.3, 5.4
        ),
        feature_B = c(
            1, 2, 3, 2, 3,
            2, 3, 4, 3, 4
        ),
        feature_C = rep(7, 10),
        feature_D = c(
            6, 6.1, 6.2, 6.3, 6.4,
            2, 2.1, 2.2, 2.3, 2.4
        )
    )
    colnames(x) <- paste0("S", seq_len(ncol(x)))
    y <- rep(c("control", "case"), each = 5)
    list(x = x, y = y)
}

test_that("topIntrinsic filters and ranks public results", {
    z <- toy_data()

    fit <- intrinsicScan(
        z$x,
        z$y,
        positive = "case",
        nResamples = 12L,
        resampleMethod = "subsample",
        resampleFraction = 0.75,
        resampleTop = Inf,
        seed = 11L,
        verbose = FALSE
    )

    top <- topIntrinsic(
        fit,
        n = 2L,
        minJ = 0,
        maxFDR = 1
    )

    expect_s4_class(top, "DataFrame")
    expect_lte(nrow(top), 2L)
    expect_true(all(top$valid))
    expect_true(all(diff(top$J) <= 0))

    stable <- topIntrinsic(
        fit,
        n = 10L,
        minJ = 0,
        maxFDR = 1,
        requireStability = TRUE,
        minHeldoutJ = 0,
        minOrientationStability = 0,
        minPartitionStability = 0,
        minExactHeldout = 0
    )

    expect_s4_class(stable, "DataFrame")
    expect_true(all(is.finite(stable$median_heldout_J)))
    expect_true(all(stable$median_heldout_J >= 0))
    expect_true(all(stable$positive_state_stability >= 0))
    expect_true(all(stable$state_partition_stability >= 0))
    expect_true(all(stable$exact_heldout_fraction >= 0))

    expect_error(
        topIntrinsic(
            fit,
            minJ = 0,
            maxFDR = 1,
            maxFWER = 1
        ),
        "run intrinsicPermutation"
    )
})


test_that("intrinsicPermutation adds reproducible permutation calibration", {
    z <- toy_data()

    fit <- intrinsicScan(
        z$x,
        z$y,
        positive = "case",
        nResamples = 0L,
        resampleMethod = "none",
        keepStates = TRUE,
        verbose = FALSE
    )

    p1 <- intrinsicPermutation(
        fit,
        nPermutations = 20L,
        cutoffs = c(0.5, 1),
        seed = 23L,
        verbose = FALSE
    )

    p2 <- intrinsicPermutation(
        fit,
        nPermutations = 20L,
        cutoffs = c(0.5, 1),
        seed = 23L,
        verbose = FALSE
    )

    d1 <- as.data.frame(intrinsicResults(p1))
    d2 <- as.data.frame(intrinsicResults(p2))

    expect_true("maxT_fwer_p" %in% names(d1))
    expect_equal(d1$maxT_fwer_p, d2$maxT_fwer_p)

    expect_true(all(
        d1$maxT_fwer_p[is.finite(d1$maxT_fwer_p)] > 0 &
        d1$maxT_fwer_p[is.finite(d1$maxT_fwer_p)] <= 1
    ))

    rd <- SummarizedExperiment::rowData(p1)
    expect_true("scan_maxT_fwer_p" %in% names(rd))

    perm <- S4Vectors::metadata(p1)$IntrinsicScan$permutation
    expect_identical(perm$n_permutations, 20L)
    expect_length(perm$null_max_J, 20L)
    expect_equal(perm$cutoffs, c(0.5, 1))
    expect_equal(perm$seed, 23L)

    top <- topIntrinsic(
        p1,
        n = 10L,
        minJ = 0,
        maxFDR = 1,
        maxFWER = 1
    )

    expect_s4_class(top, "DataFrame")
    expect_true(all(top$maxT_fwer_p <= 1))
})


test_that("intrinsicPermutation validates its inputs", {
    z <- toy_data()

    fit_no_states <- intrinsicScan(
        z$x,
        z$y,
        positive = "case",
        keepStates = FALSE,
        verbose = FALSE
    )

    expect_error(
        intrinsicPermutation(
            fit_no_states,
            nPermutations = 5L,
            verbose = FALSE
        ),
        "keepStates=TRUE"
    )

    fit <- intrinsicScan(
        z$x,
        z$y,
        positive = "case",
        keepStates = TRUE,
        verbose = FALSE
    )

    expect_error(
        intrinsicPermutation(
            fit,
            nPermutations = 0L,
            verbose = FALSE
        ),
        "positive integer"
    )

    expect_error(
        intrinsicPermutation(
            fit,
            nPermutations = 5L,
            cutoffs = c(-0.1, 1),
            verbose = FALSE
        ),
        "between 0 and 1"
    )

    expect_error(
        intrinsicPermutation(
            fit,
            nPermutations = 5L,
            seed = -1,
            verbose = FALSE
        ),
        "integer between"
    )
})


test_that("ExpressionSet input matches matrix input", {
    z <- toy_data()

    pd <- data.frame(
        group = z$y,
        row.names = colnames(z$x)
    )

    eset <- Biobase::ExpressionSet(
        assayData = z$x,
        phenoData = Biobase::AnnotatedDataFrame(pd)
    )

    a <- intrinsicScan(
        z$x,
        z$y,
        positive = "case",
        verbose = FALSE
    )

    b <- intrinsicScan(
        eset,
        phenotype = "group",
        positive = "case",
        verbose = FALSE
    )

    da <- as.data.frame(intrinsicResults(a))
    db <- as.data.frame(intrinsicResults(b))

    cols <- c(
        "feature_id",
        "threshold",
        "J",
        "positive_state",
        "exact_J_p",
        "exact_J_fdr",
        "rank_J"
    )

    expect_equal(da[, cols], db[, cols])

    meta <- S4Vectors::metadata(b)$IntrinsicScan
    expect_identical(meta$source_class, "ExpressionSet")
    expect_identical(meta$analyzed_assay, "exprs")
})


test_that("named phenotype vectors are aligned to sample IDs", {
    z <- toy_data()

    named_y <- stats::setNames(z$y, colnames(z$x))
    named_y <- named_y[rev(names(named_y))]

    a <- intrinsicScan(
        z$x,
        z$y,
        positive = "case",
        verbose = FALSE
    )

    b <- intrinsicScan(
        z$x,
        named_y,
        positive = "case",
        verbose = FALSE
    )

    da <- as.data.frame(intrinsicResults(a))
    db <- as.data.frame(intrinsicResults(b))

    expect_equal(da$threshold, db$threshold)
    expect_equal(da$J, db$J)
    expect_equal(da$exact_J_p, db$exact_J_p)
})


test_that("public input validation gives informative errors", {
    z <- toy_data()

    expect_error(
        intrinsicScan(z$x, positive = "case"),
        "phenotype must be supplied"
    )

    expect_error(
        intrinsicScan(z$x, phenotype = z$y),
        "positive must be supplied"
    )

    expect_error(
        intrinsicScan(
            z$x,
            phenotype = z$y[-1],
            positive = "case",
            verbose = FALSE
        ),
        "one value per sample"
    )

    dup_y <- stats::setNames(z$y, colnames(z$x))
    names(dup_y)[2] <- names(dup_y)[1]

    expect_error(
        intrinsicScan(
            z$x,
            phenotype = dup_y,
            positive = "case",
            verbose = FALSE
        ),
        "unique names"
    )

    incomplete_y <- stats::setNames(
        z$y[-1],
        colnames(z$x)[-1]
    )

    expect_error(
        intrinsicScan(
            z$x,
            phenotype = incomplete_y,
            positive = "case",
            verbose = FALSE
        ),
        "every sample ID"
    )

    expect_error(
        intrinsicScan(
            list(a = 1),
            phenotype = z$y,
            positive = "case",
            verbose = FALSE
        ),
        "matrix/data.frame"
    )
})


test_that("metadata access validates IntrinsicScan result objects", {
    expect_error(
        intrinsicResults(matrix(1:4, nrow = 2)),
        "IntrinsicScan SummarizedExperiment"
    )

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(x = matrix(
            1:4,
            nrow = 2,
            dimnames = list(
                c("A", "B"),
                c("S1", "S2")
            )
        ))
    )

    expect_error(
        intrinsicResults(se),
        "does not contain IntrinsicScan metadata"
    )
})

test_that("small phenotype classes trigger the intended warning", {
    x <- rbind(
        A = c(1, 1.1, 1.2, 1.3, 5, 5.1, 5.2, 5.3)
    )
    colnames(x) <- paste0("S", 1:8)
    y <- rep(c("control", "case"), each = 4)

    expect_warning(
        intrinsicScan(
            x,
            y,
            positive = "case",
            verbose = FALSE
        ),
        "fewer than 5 samples"
    )
})