# IntrinsicScan 0.99.3

- Added GitHub `URL` and `BugReports` metadata for Bioconductor submission.
- Added an explicit AI-assisted development/provenance disclosure and source-code provenance comments.
- Replaced the vignette author placeholder with the package maintainer name.
- Updated the README for the public GitHub repository and submission stage.
- No changes to the partitioning, scoring, inference, ranking, or resampling algorithms.

# IntrinsicScan 0.99.2

- Added runnable examples for all exported functions.
- Scoped RNG seeding locally to preserve deterministic numerical results without altering the caller's global RNG state.
- Updated the R dependency for Bioconductor 3.24 development.
- Expanded the vignette with required Introduction, Installation, labeled code chunks, references, and session information.
- Added the `utils::head` namespace import identified by `R CMD check`.
- No changes to the partitioning, scoring, inference, ranking, or resampling algorithms.

# IntrinsicScan 0.99.1

- Fixed Bioconductor metadata access: global `SummarizedExperiment` metadata is accessed through `S4Vectors::metadata()`, where the accessor is defined.
- Simplified `NAMESPACE` because Bioconductor accessors are called with explicit namespace qualification.
- No changes to the numerical algorithms.

# IntrinsicScan 0.99.0

* Initial Bioconductor-oriented package implementation of the IntrinsicScan method.
* Matrix, SummarizedExperiment, and ExpressionSet inputs.
* SummarizedExperiment output with ranked feature statistics and optional state assay.
* Exact-J inference, optional stability, and fixed-partition permutation calibration preserved.
