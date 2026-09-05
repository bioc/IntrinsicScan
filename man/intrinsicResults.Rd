\name{intrinsicResults}
\alias{intrinsicResults}
\title{Extract IntrinsicScan results}
\usage{intrinsicResults(x)}
\arguments{\item{x}{An IntrinsicScan result.}}
\value{An S4Vectors DataFrame of feature-level results.}
\description{Returns the ranked feature-level statistics from intrinsicScan.}
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
