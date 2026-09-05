\name{intrinsicScan}
\alias{intrinsicScan}
\title{Scan features for intrinsic binary partitions}
\usage{
intrinsicScan(x, phenotype, positive, assay = 1L, nResamples = 0L,
  resampleMethod = c("subsample", "bootstrap", "none"),
  resampleFraction = 0.8, resampleTop = 500L, seed = 1L,
  keepStates = TRUE, verbose = TRUE)
}
\arguments{
\item{x}{Numeric feature-by-sample matrix/data.frame, SummarizedExperiment, or ExpressionSet.}
\item{phenotype}{Binary phenotype vector or sample-metadata column name.}
\item{positive}{Explicit positive phenotype level.}
\item{assay}{Assay name or index for SummarizedExperiment input.}
\item{nResamples}{Number of optional internal stability resamples.}
\item{resampleMethod}{Resampling method.}
\item{resampleFraction}{Training fraction for subsampling.}
\item{resampleTop}{Number of top-ranked features to resample.}
\item{seed}{Master random seed.}
\item{keepStates}{Store full high-state assignment matrix.}
\item{verbose}{Print progress messages.}
}
\value{A SummarizedExperiment ordered by the frozen v11 ranking.}
\description{Fits a label-blind exact one-dimensional two-medoid L1 partition to each feature, then evaluates phenotype concordance after the partition is fixed.}
\examples{
x <- rbind(
    feature_A = c(1, 1.1, 1.2, 1.3, 1.4, 5, 5.1, 5.2, 5.3, 5.4),
    feature_B = c(1, 2, 3, 2, 3, 2, 3, 4, 3, 4)
)
colnames(x) <- paste0("S", seq_len(ncol(x)))
y <- rep(c("control", "case"), each = 5)
fit <- intrinsicScan(x, y, positive = "case", verbose = FALSE)
intrinsicResults(fit)
}
