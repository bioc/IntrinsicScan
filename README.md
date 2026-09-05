# IntrinsicScan

`IntrinsicScan` is a Bioconductor-oriented implementation of a label-blind,
feature-wise intrinsic binary partitioning framework.

For each continuous feature, `IntrinsicScan`:

1. fits an exact one-dimensional two-medoid partition under L1 loss using only
   the feature values;
2. reveals the binary phenotype only after the partition is fixed;
3. quantifies orientation-free phenotype concordance using Youden J;
4. provides an exact conditional hypergeometric test for the fixed-partition J;
5. optionally evaluates resampling stability with deterministic,
   feature-specific random-number streams; and
6. optionally performs fixed-partition genome-wide permutation calibration.

The package is designed to separate **intrinsic boundary discovery** from
**phenotype concordance**, **inferential evidence**, and **resampling
stability**. It is not a supervised cutoff optimizer, a general clustering
algorithm, or a bimodality test. A fitted binary partition does not by itself
establish the existence of two discrete biological states.

## Bioconductor-style input

The main function accepts a feature-by-sample matrix, `SummarizedExperiment`, or
`ExpressionSet`. It returns a `SummarizedExperiment` whose rows are ordered by
the IntrinsicScan ranking.

```r
fit <- intrinsicScan(se, phenotype = "group", positive = "case")
rowData(fit)
intrinsicResults(fit)
```

No normalization or missing-value imputation is performed. Input measurements
must already be preprocessed, continuous, complete, and finite.

## Installation

`IntrinsicScan` is currently being prepared for Bioconductor submission. The
source repository is available at:

<https://github.com/memedesimo98/IntrinsicScan>

After release through Bioconductor, installation will use:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("IntrinsicScan")
```

## Development status

Version `0.99.3` is a pre-submission Bioconductor development version. The
validated numerical engine is intentionally kept stable while package
infrastructure, documentation, and submission metadata are finalized.

## AI-assisted development disclosure

Development of **IntrinsicScan** was assisted by **OpenAI ChatGPT (GPT-5.6
Sol)**. The model contributed to iterative method critique and development, R
package engineering, testing and validation strategy, Bioconductor compliance,
documentation, and scientific drafting.

All methodological decisions, interpretation of results, validation of the
implementation, and responsibility for the accuracy, reliability, licensing,
and maintenance of the package remain with the human authors.

**Assisted-by:** OpenAI ChatGPT (GPT-5.6 Sol)

## License

IntrinsicScan is distributed under the MIT license. See `LICENSE` for details.
