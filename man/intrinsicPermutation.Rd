\name{intrinsicPermutation}
\alias{intrinsicPermutation}
\title{Fixed-partition permutation calibration}
\usage{
intrinsicPermutation(
  x,
  nPermutations = 1000L,
  cutoffs = c(0.8, 1),
  seed = 1L,
  verbose = TRUE
)
}
\arguments{
\item{x}{An IntrinsicScan result created with keepStates=TRUE.}
\item{nPermutations}{Number of label permutations.}
\item{cutoffs}{J cutoffs for global-null counts.}
\item{seed}{Random seed.}
\item{verbose}{Print progress.}
}
\value{The updated SummarizedExperiment with maxT FWER results.}
\description{Permutes phenotype labels while keeping label-blind full-sample partitions fixed.}
\examples{
x <- rbind(
    feature_A = c(1, 1.1, 1.2, 1.3, 1.4, 5, 5.1, 5.2, 5.3, 5.4),
    feature_B = c(1, 2, 3, 2, 3, 2, 3, 4, 3, 4)
)
colnames(x) <- paste0("S", seq_len(ncol(x)))
y <- rep(c("control", "case"), each = 5)
fit <- intrinsicScan(
    x, y, positive = "case", keepStates = TRUE, verbose = FALSE
)
fit <- intrinsicPermutation(
    fit, nPermutations = 10L, seed = 1L, verbose = FALSE
)
intrinsicResults(fit)
}
