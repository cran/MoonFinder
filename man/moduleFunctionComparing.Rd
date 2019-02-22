\name{moduleFunctionComparing}
\alias{moduleFunctionComparing}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compare functions of target modules
}
\description{
Plot a figure to compare the functions enriched by different modules.
}
\usage{
moduleFunctionComparing(rna, rna2mod, modList)
}
\arguments{
  \item{rna}{
an RNA symbol.
}
  \item{rna2mod}{
N x M 0-1 binary matrix of N RNA and M module. 
}
  \item{modList}{
a list of protien modules.
}
}
\value{
a dot plot comparing the functions enriched by distinct modules.
}
\references{
MoonFinder: a framework for the identification of moonlighting non-coding RNAs.
}
\author{
Lixin Cheng <easonlcheng@gmail.com>
}
\examples{
library(clusterProfiler)
data(rna2module); data(moduleList)
# moduleFunctionComparing("DIRC3",rna2module,moduleList)
# moduleFunctionComparing("CRNDE",rna2module,moduleList)
}
