\name{rna2protein}
\alias{rna2protein}
\docType{data}
\title{
rna2protein
}
\description{
two-column matrix with each row shows the binary association between RNA and protein.
}
\usage{data("rna2protein")}
\format{
  The format is:
 chr [1:12008, 1:2] "LINC00271" "MALAT1" "CRNDE" "PINK1-AS" "TUG1" "MALAT1" "CRNDE" "C11orf95" "ANCR" "LOC647323" ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:2] "Interactor1" "Interactor2"
}
\references{
MoonFinder: a framework for the identification of moonlighting non-coding RNAs.
}
\examples{
data(rna2protein)
head(rna2protein)
}
\keyword{datasets}
