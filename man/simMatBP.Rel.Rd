\name{simMatBP.Rel}
\alias{simMatBP.Rel}
\docType{data}
\title{
a semantic similarity (Rel) matrix
}
\description{
An adjacent matrix of the semantic similarity (Rel) 
among differnet protein modules.
}
\usage{data("simMatBP.Rel")}
\format{
  The format is:
 num [1:175, 1:175] 0.994 0.87 0.939 0.879 0.534 0.798 0.72 0.784 0.795 0.815 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:175] "mapped_1" "mapped_2" "mapped_3" "mapped_4" ...
  ..$ : chr [1:175] "mapped_1" "mapped_2" "mapped_3" "mapped_4" ...
}
\references{
MoonFinder: a framework for the identification of moonlighting non-coding RNAs.
}
\examples{
data(simMatBP.Rel)
head(simMatBP.Rel)
}
\keyword{datasets}
