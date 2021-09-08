
#' Metadata of the CCLE dataset
#'
#' @description This data includes metadata of cell lines included in the
#' Cancer Cell Line Encyclopaedia (CCLE) project.
#'
#' @format his data frame has 1,046 rows (cell lines) and the following
#' 13 columns:
#' \describe{
#'   \item{CCLE.name}{CCLE id}
#'   \item{Cell.line.primary.name}{cell line id}
#'   \item{Cell.line.aliases}{}
#'   \item{Gender}{"M" or "F"}
#'   \item{Site.Primary}{primary site of tumour}
#'   \item{Histology}{histology type}
#'   \item{Hist.Subtype1}{histology subtype}
#'   \item{Notes}{}
#'   \item{Source}{}
#'   \item{Expression.arrays}{file id}
#'   \item{SNP.arrays}{file id}
#'   \item{Oncomap}{}
#'   \item{Hybrid.Capture.Sequencing}{}
#' }
"CCLE_metadata"


#' Information of dataset from where genomic data of cell lines are downloaded
#'
#' @description This data includes the study where cell lines are included.
#'
#' @format This data frame has 1,946 rows and 8 columns.
"cells_mapping"


#' Precomputed dataset of 1,417 human cancer cell line profiles.
#'
#' @description This dataset contains the absolute copy number profiles of
#' 1,417 human cancer cell lines. Profiles are derived using ASCAT.
#'
#' @format This data frame has 101,637 rows and 5 columns.
#' \describe{
#'   \item{chromosome}{}
#'   \item{start}{}
#'   \item{end}{}
#'   \item{segVal}{raw copy numbers}
#'   \item{sample}{name of the sample}
#' }
"cells_segcn"


#' Length of chromosomes
#'
#' @description This dataset contains the length of chromosomes (hg19 version)
#'
#' @format vector with the size of chromosomes.
"lengthChr"
