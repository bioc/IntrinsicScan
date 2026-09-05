# IntrinsicScan prototype design notes

## Frozen method

The numerical v11 engine is copied unchanged into `R/engine-v11.R`. Public
package work should not modify the partition objective, exact-J calculation,
ranking rule, or feature-specific RNG scheme without creating a new method
version and repeating regression validation.

## Public API

`intrinsicScan()` is a thin adapter around v11. It accepts matrix,
SummarizedExperiment, and ExpressionSet inputs and returns a
SummarizedExperiment. Stability is optional and defaults to off for speed.

## Deliberately excluded

StepMiner, GMM, Youden, Messina, limma benchmarking, and external-validation
pipelines are paper-development/reproducibility machinery and are not package
functionality.

## Next technical tasks

1. Run the strict legacy-v11 regression under R.
2. Run `R CMD build`, `R CMD check`, and `BiocCheck::BiocCheck()`.
3. Decide whether copying an assay to an in-memory matrix is acceptable for the
   first release or whether DelayedArray/block processing is required.
4. Add a single-feature diagnostic plot only after the numerical API is frozen.
5. Replace placeholder maintainer metadata before any public repository/submission.
