\name{simMatBP.Resnik}
\alias{simMatBP.Resnik}
\docType{data}
\title{
a semantic similarity (Resnik) matrix
}
\description{
An adjacent matrix of the semantic similarity (Resnik) among differnet protein modules.
}
\usage{data("simMatBP.Resnik")}
\format{
  The format is:
 num [1:175, 1:175] 0.586 0.498 0.542 0.505 0.296 0.455 0.401 0.442 0.451 0.46 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:175] "mapped_1" "mapped_2" "mapped_3" "mapped_4" ...
  ..$ : chr [1:175] "mapped_1" "mapped_2" "mapped_3" "mapped_4" ...
}
\references{
MoonFinder: a framework for the identification of moonlighting non-coding RNAs.
}
\examples{
data(simMatBP.Resnik)
head(simMatBP.Resnik)
}
\keyword{datasets}
