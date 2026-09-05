# Assisted-by: OpenAI ChatGPT (GPT-5.6 Sol).
# AI assistance included method/code review, package engineering, testing
# strategy, documentation, and scientific drafting. The human authors
# reviewed and validated the implementation and retain responsibility.

#' Scan features for intrinsic binary partitions
#'
#' `intrinsicScan()` is the Bioconductor-facing interface to the frozen v11
#' label-blind partitioning engine. For each feature, it fits the exact global
#' optimum of a one-dimensional two-medoid L1 objective using feature values
#' only. Phenotype labels are used only after that partition has been fixed.
#'
#' @param x A numeric matrix/data.frame with features in rows and samples in
#'   columns, a `SummarizedExperiment`, or an `ExpressionSet`.
#' @param phenotype Either a vector of binary phenotype labels aligned to the
#'   columns of `x`, or (for Bioconductor containers) the name of a column in
#'   `colData(x)` / `pData(x)`.
#' @param positive Explicit phenotype level to treat as positive when orienting
#'   the already-fixed low/high partition.
#' @param assay Assay name or index when `x` is a `SummarizedExperiment`.
#' @param nResamples Number of internal stability resamples. The default is zero
#'   because stability is optional and can be computationally expensive.
#' @param resampleMethod One of `"none"`, `"subsample"`, or `"bootstrap"`.
#' @param resampleFraction Training fraction for subsampling.
#' @param resampleTop Number of top-ranked features for which stability is
#'   assessed. May be `Inf`.
#' @param seed Master random seed. Feature-specific stability seeds are derived
#'   deterministically from this seed and the feature identifier.
#' @param keepStates Logical; store the full high/low state assignment matrix as
#'   assay `"state"` in the returned object.
#' @param verbose Logical; print progress messages.
#'
#' @return A `SummarizedExperiment`. The analyzed expression matrix is stored in
#'   assay `"expression"`; optional high-state assignments are stored in assay
#'   `"state"`. Features are ordered by the frozen v11 ranking. Scan statistics
#'   are available through `intrinsicResults()` and are also appended to
#'   `rowData()` with prefix `scan_`.
#'
#' @details No normalization, transformation, missing-value imputation, or
#' feature filtering is performed. Inputs must be complete, finite,
#' preprocessed continuous measurements. A high phenotype-concordance score
#' does not by itself establish bimodality or discrete biological states.
#'
#' @examples
#' x <- rbind(
#'   feature_A = c(1, 1.1, 1.2, 1.3, 1.4, 5, 5.1, 5.2, 5.3, 5.4),
#'   feature_B = c(1, 2, 3, 2, 3, 2, 3, 4, 3, 4)
#' )
#' colnames(x) <- paste0("S", seq_len(ncol(x)))
#' y <- rep(c("control", "case"), each = 5)
#' fit <- intrinsicScan(x, y, positive = "case", verbose = FALSE)
#' intrinsicResults(fit)
#' @export
intrinsicScan <- function(x,
                          phenotype,
                          positive,
                          assay = 1L,
                          nResamples = 0L,
                          resampleMethod = c("subsample", "bootstrap", "none"),
                          resampleFraction = 0.80,
                          resampleTop = 500L,
                          seed = 1L,
                          keepStates = TRUE,
                          verbose = TRUE) {
    if (missing(phenotype)) stop("phenotype must be supplied")
    if (missing(positive)) stop("positive must be supplied explicitly")
    resampleMethod <- match.arg(resampleMethod)

    inp <- .intrinsic_input(x, phenotype = phenotype, assay = assay)

    eng <- label_blind_state_concordance_v11(
        feature_matrix = inp$matrix,
        labels = inp$labels,
        positive_label = positive,
        n_resamples = nResamples,
        resample_method = resampleMethod,
        resample_fraction = resampleFraction,
        resample_top_n = resampleTop,
        seed = seed,
        keep_high_state_matrix = keepStates,
        verbose = verbose
    )

    ranked_ids <- eng$summary$feature_id
    idx <- match(ranked_ids, rownames(inp$matrix))
    if (anyNA(idx)) stop("internal feature alignment failure")
    expr_ranked <- inp$matrix[idx, , drop = FALSE]

    rd0 <- inp$rowData[idx, , drop = FALSE]
    scan_df <- S4Vectors::DataFrame(eng$summary)
    rownames(scan_df) <- ranked_ids
    scan_rd <- scan_df
    names(scan_rd) <- paste0("scan_", names(scan_rd))
    rd <- cbind(rd0, scan_rd)
    rownames(rd) <- ranked_ids

    cd <- inp$colData
    cd$IntrinsicScanLabel <- eng$labels

    assays_out <- list(expression = expr_ranked)
    if (isTRUE(keepStates)) {
        st <- eng$high_state_matrix[ranked_ids, , drop = FALSE]
        storage.mode(st) <- "logical"
        assays_out$state <- st
    }

    out <- SummarizedExperiment::SummarizedExperiment(
        assays = assays_out,
        rowData = rd,
        colData = cd
    )

    md <- S4Vectors::metadata(out)
    md$IntrinsicScan <- list(
        engine = "v11",
        method = eng$method,
        positive_label = eng$positive_label,
        negative_label = eng$negative_label,
        parameters = eng$parameters,
        analyzed_assay = inp$assay_name,
        results = scan_df,
        resampling_details = eng$resampling_details,
        source_class = inp$source_class
    )
    metadata(out) <- md
    out
}

#' Extract IntrinsicScan feature results
#' @param x A result returned by `intrinsicScan()`.
#' @return An `S4Vectors::DataFrame`, ordered by the frozen v11 ranking.
#' @examples
#' x <- rbind(
#'   feature_A = c(1, 1.1, 1.2, 1.3, 1.4, 5, 5.1, 5.2, 5.3, 5.4),
#'   feature_B = c(1, 2, 3, 2, 3, 2, 3, 4, 3, 4)
#' )
#' colnames(x) <- paste0("S", seq_len(ncol(x)))
#' y <- rep(c("control", "case"), each = 5)
#' fit <- intrinsicScan(x, y, positive = "case", verbose = FALSE)
#' intrinsicResults(fit)
#' @export
intrinsicResults <- function(x) {
    meta <- .intrinsic_metadata(x)
    meta$results
}

#' Filter and rank IntrinsicScan results
#'
#' This helper does not create a composite classifier label. Effect size,
#' statistical evidence, and stability remain separate result axes.
#'
#' @param x A result returned by `intrinsicScan()`.
#' @param n Maximum number of rows to return.
#' @param minJ Minimum orientation-free Youden J.
#' @param maxFDR Maximum exact-J Benjamini-Hochberg FDR.
#' @param requireStability Require a finite held-out stability estimate.
#' @param minHeldoutJ Optional minimum median held-out J.
#' @param minOrientationStability Optional minimum positive-state stability.
#' @param minPartitionStability Optional minimum state-partition stability.
#' @param minExactHeldout Optional minimum exact-heldout fraction.
#' @param maxFWER Optional maximum maxT FWER P value; requires prior
#'   `intrinsicPermutation()`.
#' @return A filtered `S4Vectors::DataFrame`.
#' @examples
#' x <- rbind(
#'   feature_A = c(1, 1.1, 1.2, 1.3, 1.4, 5, 5.1, 5.2, 5.3, 5.4),
#'   feature_B = c(1, 2, 3, 2, 3, 2, 3, 4, 3, 4)
#' )
#' colnames(x) <- paste0("S", seq_len(ncol(x)))
#' y <- rep(c("control", "case"), each = 5)
#' fit <- intrinsicScan(x, y, positive = "case", verbose = FALSE)
#' topIntrinsic(fit, n = 2, minJ = 0, maxFDR = 1)
#' @export
topIntrinsic <- function(x, n = 50L, minJ = 0.80, maxFDR = 0.05,
                         requireStability = FALSE,
                         minHeldoutJ = NULL,
                         minOrientationStability = NULL,
                         minPartitionStability = NULL,
                         minExactHeldout = NULL,
                         maxFWER = NULL) {
    d <- as.data.frame(intrinsicResults(x))
    keep <- d$valid & is.finite(d$J) & d$J >= minJ &
        is.finite(d$exact_J_fdr) & d$exact_J_fdr <= maxFDR
    if (requireStability) keep <- keep & is.finite(d$median_heldout_J)
    if (!is.null(minHeldoutJ)) keep <- keep & is.finite(d$median_heldout_J) & d$median_heldout_J >= minHeldoutJ
    if (!is.null(minOrientationStability)) keep <- keep & is.finite(d$positive_state_stability) & d$positive_state_stability >= minOrientationStability
    if (!is.null(minPartitionStability)) keep <- keep & is.finite(d$state_partition_stability) & d$state_partition_stability >= minPartitionStability
    if (!is.null(minExactHeldout)) keep <- keep & is.finite(d$exact_heldout_fraction) & d$exact_heldout_fraction >= minExactHeldout
    if (!is.null(maxFWER)) {
        if (!"maxT_fwer_p" %in% names(d)) stop("maxT_fwer_p unavailable; run intrinsicPermutation() first")
        keep <- keep & is.finite(d$maxT_fwer_p) & d$maxT_fwer_p <= maxFWER
    }
    ans <- d[keep, , drop = FALSE]
    ans <- ans[order(-ans$J, ans$exact_J_fdr, -ans$median_heldout_J, na.last = TRUE), , drop = FALSE]
    S4Vectors::DataFrame(utils::head(ans, as.integer(n)))
}

#' Fixed-partition genome-wide permutation calibration
#'
#' Permutes phenotype labels while keeping the full-sample label-blind state
#' partitions fixed. This calibrates fixed-partition phenotype concordance; it
#' is not a permutation test of feature selection or the stability workflow.
#'
#' @param x An `intrinsicScan()` result created with `keepStates = TRUE`.
#' @param nPermutations Number of label permutations.
#' @param cutoffs J cutoffs for global-null feature-count summaries.
#' @param seed Random seed.
#' @param verbose Logical; print progress.
#' @return The input `SummarizedExperiment` with `scan_maxT_fwer_p` appended to
#'   `rowData()` and permutation details stored in `metadata()`.
#' @examples
#' x <- rbind(
#'   feature_A = c(1, 1.1, 1.2, 1.3, 1.4, 5, 5.1, 5.2, 5.3, 5.4),
#'   feature_B = c(1, 2, 3, 2, 3, 2, 3, 4, 3, 4)
#' )
#' colnames(x) <- paste0("S", seq_len(ncol(x)))
#' y <- rep(c("control", "case"), each = 5)
#' fit <- intrinsicScan(
#'   x, y, positive = "case", keepStates = TRUE, verbose = FALSE
#' )
#' fit <- intrinsicPermutation(
#'   fit, nPermutations = 10L, seed = 1L, verbose = FALSE
#' )
#' intrinsicResults(fit)
#' @export
intrinsicPermutation <- function(x, nPermutations = 1000L,
                                 cutoffs = c(0.80, 1.00), seed = 1L,
                                 verbose = TRUE) {
    meta <- .intrinsic_metadata(x)
    if (!"state" %in% SummarizedExperiment::assayNames(x)) {
        stop("x must have been created with keepStates=TRUE")
    }

    legacy <- list(
        summary = as.data.frame(meta$results),
        high_state_matrix = SummarizedExperiment::assay(x, "state"),
        labels = SummarizedExperiment::colData(x)$IntrinsicScanLabel,
        positive_label = meta$positive_label,
        negative_label = meta$negative_label
    )
    class(legacy) <- c("lbsc_v11", "list")

    p <- permutation_global_calibration_v11(
        legacy,
        n_permutations = nPermutations,
        cutoffs = cutoffs,
        seed = seed,
        verbose = verbose
    )

    idx <- match(rownames(x), p$per_feature$feature_id)
    fwer <- p$per_feature$maxT_fwer_p[idx]
    rd <- SummarizedExperiment::rowData(x)
    rd$scan_maxT_fwer_p <- fwer
    rowData(x) <- rd

    rr <- meta$results
    rr$maxT_fwer_p <- fwer[match(rr$feature_id, rownames(x))]
    md <- S4Vectors::metadata(x)
    md$IntrinsicScan$results <- rr
    md$IntrinsicScan$permutation <- p
    metadata(x) <- md
    x
}

.intrinsic_metadata <- function(x) {
    if (!inherits(x, "SummarizedExperiment")) stop("x must be an IntrinsicScan SummarizedExperiment result")
    z <- S4Vectors::metadata(x)$IntrinsicScan
    if (is.null(z)) stop("x does not contain IntrinsicScan metadata")
    z
}

.intrinsic_labels <- function(phenotype, sample_ids, data = NULL) {
    if (length(phenotype) == 1L && is.character(phenotype) && !is.null(data)) {
        if (!phenotype %in% colnames(data)) stop("phenotype column not found in sample metadata")
        return(data[[phenotype]])
    }
    if (!is.null(names(phenotype))) {
        if (anyDuplicated(names(phenotype))) stop("named phenotype vector must have unique names")
        if (!all(sample_ids %in% names(phenotype))) {
            stop("named phenotype vector must contain every sample ID")
        }
        phenotype <- phenotype[sample_ids]
    }
    if (length(phenotype) != length(sample_ids)) {
        stop("phenotype must have one value per sample, or name a sample-metadata column")
    }
    phenotype
}

.intrinsic_input <- function(x, phenotype, assay = 1L) {
    if (inherits(x, "SummarizedExperiment")) {
        mat <- as.matrix(SummarizedExperiment::assay(x, assay))
        cd <- S4Vectors::DataFrame(SummarizedExperiment::colData(x))
        rd <- S4Vectors::DataFrame(SummarizedExperiment::rowData(x))
        labs <- .intrinsic_labels(phenotype, colnames(mat), data = as.data.frame(cd))
        an <- SummarizedExperiment::assayNames(x)
        assay_name <- if (is.character(assay)) assay else if (length(an) >= assay) an[[assay]] else as.character(assay)
        return(list(matrix = mat, labels = labs, colData = cd, rowData = rd,
                    source_class = class(x)[1L], assay_name = assay_name))
    }

    if (inherits(x, "ExpressionSet")) {
        mat <- as.matrix(Biobase::exprs(x))
        pd <- Biobase::pData(x)
        fd <- Biobase::fData(x)
        cd <- S4Vectors::DataFrame(pd, row.names = rownames(pd))
        rd <- S4Vectors::DataFrame(fd, row.names = rownames(fd))
        labs <- .intrinsic_labels(phenotype, colnames(mat), data = pd)
        return(list(matrix = mat, labels = labs, colData = cd, rowData = rd,
                    source_class = class(x)[1L], assay_name = "exprs"))
    }

    if (is.matrix(x) || is.data.frame(x)) {
        mat <- as.matrix(x)
        cd <- S4Vectors::DataFrame(row.names = colnames(mat))
        rd <- S4Vectors::DataFrame(row.names = rownames(mat))
        labs <- .intrinsic_labels(phenotype, colnames(mat), data = NULL)
        return(list(matrix = mat, labels = labs, colData = cd, rowData = rd,
                    source_class = class(x)[1L], assay_name = "matrix"))
    }

    stop("x must be a matrix/data.frame, SummarizedExperiment, or ExpressionSet")
}
