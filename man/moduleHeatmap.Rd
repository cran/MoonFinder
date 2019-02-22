\name{moduleHeatmap}
\alias{moduleHeatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
moduleHeatmap
}
\description{
Draw a heatmap to illustrate the protein modules mediated by an RNA.
}
\usage{
moduleHeatmap(rna, rna2mod, modList, net)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rna}{
RNA symbol.
}
  \item{rna2mod}{
N x M 0-1 binary matrix of N RNA and M module.
}
  \item{modList}{
a list of protien moudles.
}
  \item{net}{
a protein-protein interaction network compiled in igraph.
}
}
\value{
a heatmap figure indicating protein modules targeted by an RNA.
}
\references{
MoonFinder: a framework for the identification of moonlighting non-coding RNAs.
}
\author{
Lixin Cheng <easonlcheng@gmail.com>
}
\examples{
library(igraph)
data(ppi)
data(rna2module)
data(moduleList)
g <- graph.edgelist(ppi,directed=FALSE)
moduleHeatmap("DIRC3",rna2module,moduleList,g)
# moduleHeatmap("CRNDE",rna2module,moduleList,g)
}
