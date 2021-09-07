
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CNpare

## Summary of the package

**CNpare** is designed to compare and contrast genome-wide tumour copy
number profiles. By using CNpare, one can:

-   Calculate similarities between copy-number profiles
-   Quantify the percentage of genome differences between copy-number
    profiles
-   Compare copy-number profiles based on copy-number signatures

## Installation

You can install the development version of CNpare from this github
repository using devtools package. Once installed, CNpare can be loaded:

``` r
library(CNpare)
```

## Workflow overview

### Input data

The most basic initial input to the CNpare package consists of a data
frame containing the copy number profile segmented from a sample.

This structure must contain the following columns:

-   chromosome – `chromosome`
-   start position – `start`
-   end position – `end`
-   copy number – `segVal`
-   sample identifier – `sample`

Alternatively, a `QDNAseqCopyNumbers` object can be also used as input.

The segmented profile used as input needs to be transformed to bin-level
copy number profiles. This can be done using `getCNbins()`. Genomic
positions of bins can be obtained using `getBinsStartsEnds()`.

This package includes `CellModels_ASCAT.RData`, a dataframe with
copy-number profiles of 1,417 human cancer cell lines from the Cancer
Cell Line Encyclopaedia (CCLE) and Genomics of Drug Sensitivity in
cancer (GDSC), which has been generated with ASCAT. This dataset must
also be transformed to bin-level copy numbers before to use for
comparisons.

### Compare two copy-number profiles

The function `getSimilarities()` can be used to get all similarity
metrics calculated by CNpare. This function needs at least two bin-level
copy number profiles to compare.

Four similarity metrics are computed for comparison by default, but the
user can select to compute only one by modifying the `method` parameter
as follows:

-   Pearson correlation – `method=pearson`
-   Manhattan distance – `method=manhattan`
-   Euclidean distance – `method=euclidean`
-   Cosine similarity – `method=cosine`

### Compute the % genome differences

Using the function `getDifference()`, the extent of genome differences
(%) between two profiles can be calculated

### Cluster profiles based on their exposure to copy-number signature

First, CNpare extracts the activity levels of the 7 ovarian copy-number
signatures using the computational approach developed by *Macintyre el
al. Nature Genetics, 2018*. Functions used for this purpose are:

-   `extractCopynumberFeatures()`
-   `generateSampleByComponentMatrix()`
-   `quantifySignatures()`

The, using the function `plotClusters()`, profiles are clustered
according to their exposure to each copy number signature.

## Manuscript analysis Rmarkdown

The two `Rmarkdown` documents can be compiled to reproduce the entire
analysis carried out in the manuscript:

-   `CNpare_AssessingTool` – this tutorial reproduces the analyses
    performed to check the ability of the tool to find the most similar
    profile
-   `CNpare_DemonstratingUtility` – this tutorial follows the typical
    workflow of the tool

To compile this document run knitr to html. Packages required to compile
the Rmarkdown document include: knitr, dplyr, ggplot2, reshape2,
data.table, factoextra, flexmix, NMF, stats, splitstackshape, kableExtra
and magrittr.
