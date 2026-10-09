test_that("clean two-state feature gives exact concordance", {
    x <- rbind(
        feature_A = c(
            1, 1.1, 1.2, 1.3, 1.4,
            5, 5.1, 5.2, 5.3, 5.4
        ),
        feature_B = c(
            1, 2, 3, 2, 3,
            2, 3, 4, 3, 4
        ),
        feature_C = rep(7, 10)
    )
    colnames(x) <- paste0("S", 1:10)
    y <- rep(c("control", "case"), each = 5)

    z <- intrinsicScan(
        x,
        y,
        positive = "case",
        nResamples = 0,
        resampleMethod = "none",
        verbose = FALSE
    )

    d <- as.data.frame(intrinsicResults(z))
    a <- d[d$feature_id == "feature_A", ]

    expect_equal(a$threshold, 3.2, tolerance = 1e-12)
    expect_equal(a$J, 1, tolerance = 1e-12)
    expect_identical(a$positive_state, "high")
    expect_true(a$exact_full_concordance)
    expect_false(d$valid[d$feature_id == "feature_C"])
})


test_that("phenotype orientation changes without changing partition", {
    x <- rbind(
        feature_A = c(
            1, 1.1, 1.2, 1.3, 1.4,
            5, 5.1, 5.2, 5.3, 5.4
        )
    )
    colnames(x) <- paste0("S", 1:10)
    y <- rep(c("control", "case"), each = 5)

    a <- intrinsicScan(
        x,
        y,
        positive = "case",
        verbose = FALSE
    )

    b <- intrinsicScan(
        x,
        y,
        positive = "control",
        verbose = FALSE
    )

    da <- as.data.frame(intrinsicResults(a))
    db <- as.data.frame(intrinsicResults(b))

    expect_equal(da$threshold, db$threshold)
    expect_equal(da$J, db$J)
    expect_identical(da$positive_state, "high")
    expect_identical(db$positive_state, "low")
})


test_that("SummarizedExperiment input matches matrix input", {
    x <- rbind(
        feature_A = c(
            1, 1.1, 1.2, 1.3, 1.4,
            5, 5.1, 5.2, 5.3, 5.4
        ),
        feature_B = c(
            1, 2, 3, 2, 3,
            2, 3, 4, 3, 4
        )
    )
    colnames(x) <- paste0("S", 1:10)
    y <- rep(c("control", "case"), each = 5)

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(logexpr = x),
        colData = S4Vectors::DataFrame(
            group = y,
            row.names = colnames(x)
        )
    )

    a <- intrinsicScan(
        x,
        y,
        positive = "case",
        verbose = FALSE
    )

    b <- intrinsicScan(
        se,
        "group",
        positive = "case",
        assay = "logexpr",
        verbose = FALSE
    )

    da <- as.data.frame(intrinsicResults(a))
    db <- as.data.frame(intrinsicResults(b))

    cols <- c(
        "feature_id",
        "threshold",
        "J",
        "exact_J_p",
        "exact_J_fdr",
        "rank_J"
    )

    expect_equal(da[, cols], db[, cols])
})


test_that("feature-specific stability survives row reorder", {
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
        )
    )
    colnames(x) <- paste0("S", 1:10)
    y <- rep(c("control", "case"), each = 5)

    a <- intrinsicScan(
        x,
        y,
        positive = "case",
        nResamples = 20,
        resampleMethod = "subsample",
        resampleFraction = 0.8,
        resampleTop = Inf,
        seed = 17,
        verbose = FALSE
    )

    b <- intrinsicScan(
        x[c("C", "A", "B"), ],
        y,
        positive = "case",
        nResamples = 20,
        resampleMethod = "subsample",
        resampleFraction = 0.8,
        resampleTop = Inf,
        seed = 17,
        verbose = FALSE
    )

    da <- as.data.frame(intrinsicResults(a))
    db <- as.data.frame(intrinsicResults(b))

    ia <- match(c("A", "B", "C"), da$feature_id)
    ib <- match(c("A", "B", "C"), db$feature_id)

    expect_equal(
        da$stability_seed[ia],
        db$stability_seed[ib]
    )
    expect_equal(
        da$median_heldout_J[ia],
        db$median_heldout_J[ib]
    )
    expect_equal(
        da$state_partition_stability[ia],
        db$state_partition_stability[ib]
    )
})