\name{topIntrinsic}
\alias{topIntrinsic}
\title{Filter top IntrinsicScan features}
\usage{
topIntrinsic(x, n = 50L, minJ = 0.8, maxFDR = 0.05,
  requireStability = FALSE, minHeldoutJ = NULL,
  minOrientationStability = NULL, minPartitionStability = NULL,
  minExactHeldout = NULL, maxFWER = NULL)
}
\arguments{
\item{x}{An IntrinsicScan result.} \item{n}{Maximum rows.}
\item{minJ}{Minimum Youden J.} \item{maxFDR}{Maximum exact-J FDR.}
\item{requireStability}{Require stability estimates.}
\item{minHeldoutJ}{Optional minimum median held-out J.}
\item{minOrientationStability}{Optional minimum orientation stability.}
\item{minPartitionStability}{Optional minimum state-partition stability.}
\item{minExactHeldout}{Optional minimum exact-heldout fraction.}
\item{maxFWER}{Optional maximum maxT FWER P value.}
}
\value{A filtered S4Vectors DataFrame.}
\description{Filters and ranks results without collapsing effect, evidence, and stability into a composite category.}
\examples{
x <- rbind(
    feature_A = c(1, 1.1, 1.2, 1.3, 1.4, 5, 5.1, 5.2, 5.3, 5.4),
    feature_B = c(1, 2, 3, 2, 3, 2, 3, 4, 3, 4)
)
colnames(x) <- paste0("S", seq_len(ncol(x)))
y <- rep(c("control", "case"), each = 5)
fit <- intrinsicScan(x, y, positive = "case", verbose = FALSE)
topIntrinsic(fit, n = 2, minJ = 0, maxFDR = 1)
}
