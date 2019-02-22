\name{simMatBP.Jiang}
\alias{simMatBP.Jiang}
\docType{data}
\title{
a semantic similarity (Jiang) matrix
}
\description{
An adjacent matrix of the semantic similarity (Jiang) among differnet protein modules.
}
\usage{data("simMatBP.Jiang")}
\format{
  The format is:
 num [1:175, 1:175] 1 0.857 0.935 0.868 0.547 0.78 0.713 0.775 0.779 0.807 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:175] "mapped_1" "mapped_2" "mapped_3" "mapped_4" ...
  ..$ : chr [1:175] "mapped_1" "mapped_2" "mapped_3" "mapped_4" ...
}
\references{
MoonFinder: a framework for the identification of moonlighting non-coding RNAs.
}
\examples{
data(simMatBP.Jiang)
dim(simMatBP.Jiang)
}
\keyword{datasets}
