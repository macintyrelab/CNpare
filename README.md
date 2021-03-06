
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CNpare

## Summary of the package

**CNpare** is designed to compare and contrast genome-wide tumour copy
number profiles. CNpare can be used for:

-   Calculating similarities between copy-number profiles
-   Quantifying the percentage of genome differences between copy-number
    profiles
-   Comparing copy-number profiles based on copy-number signatures

## Installation

You can install the development version of CNpare as follows:

``` r
git clone https://github.com/macintyrelab/CNpare.git
```

Once installed, CNpare must be loaded

``` r
library(CNpare)
```

## Workflow overview

### Input data

The most basic initial input to the CNpare package consists of a data
frame containing the segmented copy number profile from a tumor sample.

This structure must contain the following columns:

-   chromosome – `chromosome`
-   start position – `start`
-   end position – `end`
-   copy number – `segVal`
-   sample identifier – `sample`

Alternatively, a `QDNAseqCopyNumbers` object can be also used as input.

The segmented profile used as input needs to be transformed to bin-level
copy number profiles. This can be done using `getCNbins()`, which needs
the following arguments:

-   `posBins` – genomic positions of bins. These can be obtained using
    `getBinsStartsEnds()`
-   `data` – segmented copy number
-   `samples` – vector with the name of samples

### Package data

This package includes `cells_segcn.RData`, a dataframe with copy-number
profiles of 1,417 human cancer cell lines from the Cancer Cell Line
Encyclopaedia (CCLE) and Genomics of Drug Sensitivity in cancer (GDSC),
which has been generated with `ASCAT.sc`. The original file can be found
in the following github repository: `VanLoo-lab/ASCAT.sc`
(<https://github.com/VanLoo-lab/ASCAT.sc>). This dataset must also be
transformed to bin-level copy numbers before to be used for comparisons.

### Compare two copy-number profiles

The function `getSimilarities()` can be used to get all similarity
metrics calculated by CNpare. This function needs at least two bin-level
copy number profiles to compare.

Four similarity metrics are computed for comparison by default, but the
user can be interested in only computing one of them. It can be done by
modifying the `method` parameter as follows:

-   Pearson correlation – `method="pearson"`
-   Manhattan distance – `method="manhattan"`
-   Euclidean distance – `method="euclidean"`
-   Cosine similarity – `method="cosine"`

The user can choose to compute the empirical p-value of each comparison
by adding the argument `pvalue=TRUE`. To calculate the empirical
p-value, CNpare uses the distribution of similarities obtained after
comparing all samples included in the CNpare dataset. The function
`plot_simdensity` can be used to visually inspect the uniqueness of the
similarity value respect to the observed values.

### Compute the % genome differences

Using the function `getDifference()`, the extent of genome differences
(%) between two profiles can be calculated. This function provides the
user the option to calculate the % genome difference normalized by
ploidy status by adding the argument `method="normalized"`

### Cluster profiles based on their exposure to copy-number signature

First, CNpare extracts the activity levels of the 7 ovarian copy-number
signatures using the computational approach developed by *Macintyre el
al. Nature Genetics, 2018*. Functions used for this purpose are:

-   `extractCopynumberFeatures()`
-   `generateSampleByComponentMatrix()`
-   `quantifySignatures()`

Data needed for quantifying these signatures is included in CNpare,
except for `component_parameters.rds` which may be downloaded from here:
<https://bitbucket.org/britroc/cnsignatures/src/master/data/>. We are
preparing an experiment data package for including all data needed for
identifying copy-number signatures

The `getClusterSamples()` function can be used to get names of samples
included in the same (`include=TRUE`, default) or closest
(`include=FALSE`) cluster to the sample of interest according to the
exposure to each copy number signature. This function automatically
identifies the optimal number of clusters. Cosine similarity is used as
distance metric for clustering.

Finally, using the function `plotClusters()`, samples are plotted by
showing the two signatures with the highest variation across cluster
means. Clusters are also represented.

## Manuscript analysis Rmarkdown

The `CNpare_Workflow` rmarkdown in the `/vignettes` folder can be
compiled to reproduce a typical workflow of the tool.

The `Rmarkdown` documents in the `macintyrelab/CNpare_analyses`
repository (<https://github.com/macintyrelab/CNpare_analyses>) can be
compiled to reproduce the analyses performed in the manuscript.

To compile these documents run knitr to html. Packages required to
compile Rmarkdown documents include: knitr, dplyr, ggplot2, reshape2,
data.table, flexmix, NMF, stats, splitstackshape, kableExtra, qusage,
string and magrittr.
