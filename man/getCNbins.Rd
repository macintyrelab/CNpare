% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataprocessing_functions.R
\name{getCNbins}
\alias{getCNbins}
\title{Get copy-number per bins}
\usage{
getCNbins(posBins, data, samples)
}
\arguments{
\item{posBins}{list with genomic positions of bins. Each list contains
data from each chromosome. Obtained from getBinsStartsEnds}

\item{data}{segment table of copy numbers of all samples}

\item{samples}{vector with sample names}
}
\value{
bin table of copy numbers of all samples
}
\description{
This function transform segment tables to bin tables
}
\examples{
posBins <- lapply(1,function(chr)
    getBinsStartsEnds(window=1000000, chr, lengthChr[chr]))
cells_bin <- getCNbins(posBins=posBins,
    data=cells_segcn[cells_segcn$sample=="22RV1",], samples="22RV1")
}
