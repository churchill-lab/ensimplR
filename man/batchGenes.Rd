% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ensimpl_functions.R
\name{batchGenes}
\alias{batchGenes}
\title{Get batch gene information.}
\usage{
batchGenes(ids, species, release, details = FALSE, debug = FALSE)
}
\arguments{
\item{ids}{A list of Ensembl identifiers.}

\item{species}{Either 'Mm' or 'Hs'.}

\item{release}{ensimpl release, see releases()}

\item{details}{TRUE for extra information.}

\item{debug}{Display some debugging information.}
}
\value{
By default, with details=FALSE, a tibble will be returned.
  If details=TRUE, a list will be returned with each element being
  the matching "id".
}
\description{
Get batch gene information.
}
\examples{
genes <- batchGenes(c("ENSMUSG00000000001", "ENSMUSG00000000002"),
                    "Mm", "101")
}
