\name{ppi}
\alias{ppi}
\docType{data}
\title{
ppi
}
\description{
two-colomn matrix with each row shows interaction between a pair of proteins.
}
\usage{data("ppi")}
\format{
  The format is:
 chr [1:210410, 1:2] "KANSL1L" "SLC12A8" "RBM47" "RBM47" "RBM47" "RBM47" "RBM47" "RBM47" "TTC26" "TTC26" "TTC26" ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:2] "V1" "V2"
}
\references{
MoonFinder: a framework for the identification of moonlighting non-coding RNAs.
}
\examples{
data(ppi)
head(ppi)
}
\keyword{datasets}
