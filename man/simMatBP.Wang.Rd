\name{simMatBP.Wang}
\alias{simMatBP.Wang}
\docType{data}
\title{
a semantic similarity (Wang) matrix~
}
\description{
An adjacent matrix of the semantic similarity (Wang) among differnet protein modules.
}
\usage{data("simMatBP.Wang")}
\format{
  The format is:
 num [1:175, 1:175] 1 0.827 0.919 0.845 0.527 0.761 0.659 0.725 0.742 0.77 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:175] "mapped_1" "mapped_2" "mapped_3" "mapped_4" ...
  ..$ : chr [1:175] "mapped_1" "mapped_2" "mapped_3" "mapped_4" ...
}
\references{
MoonFinder: a framework for the identification of moonlighting non-coding RNAs.
}
\examples{
data(simMatBP.Wang)
head(simMatBP.Wang)
}
\keyword{datasets}
