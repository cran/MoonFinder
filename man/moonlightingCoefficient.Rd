\name{moonlightingCoefficient}
\alias{moonlightingCoefficient}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
moonlightingCoefficient
}
\description{
Calculate the moonlighting coefficient for RNAs.
}
\usage{
moonlightingCoefficient(rna2mod, simMat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rna2mod}{two-column matrix with each row shows the 
      binary association between RNA and protein.}
  \item{simMat}{adjacent matrix indicating the semantic similarity 
      among protein modules}
}
\value{
  \item{rnaStat}{a table including ncRNA symble, number of target modules, 
      and coresponding moonlighting coefficient.}
}
\references{
MoonFinder: a framework for the identification of moonlighting non-coding RNAs.
}
\author{
Lixin Cheng <easonlcheng@gmail.com>
}
\examples{
data(rna2module)
data(simMatBP.Rel)
mcOutput = moonlightingCoefficient(rna2module,simMat=simMatBP.Rel)
head(mcOutput)
}
