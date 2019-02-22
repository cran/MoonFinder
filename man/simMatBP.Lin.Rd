\name{simMatBP.Lin}
\alias{simMatBP.Lin}
\docType{data}
\title{
a semantic similarity (Lin) matrix
}
\description{
An adjacent matrix of the semantic similarity (Lin) among differnet protein modules.
}
\usage{data("simMatBP.Lin")}
\format{
  The format is:
 num [1:175, 1:175] 1 0.877 0.945 0.886 0.561 0.808 0.733 0.793 0.805 0.825 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:175] "mapped_1" "mapped_2" "mapped_3" "mapped_4" ...
  ..$ : chr [1:175] "mapped_1" "mapped_2" "mapped_3" "mapped_4" ...
}
\references{
MoonFinder: a framework for the identification of moonlighting non-coding RNAs.
}
\examples{
data(simMatBP.Lin)
head(simMatBP.Lin)
}
\keyword{datasets}
