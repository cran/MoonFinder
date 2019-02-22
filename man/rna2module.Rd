\name{rna2module}
\alias{rna2module}
\docType{data}
\title{
RNA-module association matrix
}
\description{
A 0-1 binary matrix representing the association between RNAs and modules
}
\usage{data("rna2module")}
\format{
  The format is:
 num [1:776, 1:175] 0 0 0 0 0 0 0 0 0 0 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:776] "CRNDE" "NEAT1" "ZFAS1" "GAS5" ...
  ..$ : chr [1:175] "1" "2" "3" "4" ...
}
\references{
MoonFinder: a framework for the identification of moonlighting non-coding RNAs.
}
\examples{
data(rna2module)
head(rna2module)
}
\keyword{datasets}
