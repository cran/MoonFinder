\name{MoonFinder-package}
\alias{MoonFinder-package}
\alias{MoonFinder}
\docType{package}
\title{
MoonFinder
}
\description{
Identify moonlighting non-coding RNAs
}
\details{
MoonFinder is developed for the identification of mncRNAs at the genome level based on the functional annotation and interactome data of ncRNAs and proteins. 
MoonFinder provides an extensive framework for the detection of protein modules and establishment of RNA-module associations.
}
\author{
Maintainer: Lixin Cheng <lixincheng@cuhk.edu.hk>
}
\references{
MoonFinder: a framework for the identification of moonlighting non-coding RNAs.
}
\keyword{ package }
\examples{
## 1) Identify the moonlight lncRNAs
# install.packages("ggraph") 
data("ppi","simMatBP.Resnik","moduleList","rna2module")
rnaStatMat = moonlightingCoefficient(rna2module,simMatBP.Resnik)
head(rnaStatMat)

## 2) Set up the SDI threshold according to SDI distributions

sdi = as.numeric(rnaStatMat[,3])
dev.new(width=4, height=4)
hist(sdi,border=FALSE,xlab="MC",xlim=c(0,1),ylim=c(0,4),
    breaks=seq(0,1,0.02),prob=TRUE,main="Resnik")
lines(density(sdi),lwd=2,col="orange")
mlncRNA = rnaStatMat[sdi>0.4,1]

## 3) Heatmap 
# Illustration of the mediated modules for a given moonlighting ncRNA.
library(igraph)
g <- graph.edgelist(ppi,directed=FALSE)
moduleHeatmap("DIRC3",rna2module,moduleList,g)
# moduleHeatmap("CRNDE",rna2module,moduleList,g)

## 4) Dotplot
# Comparison of the functions enriched by different modules:
# library(clusterProfiler)
# moduleFunctionComparing("DIRC3",rna2module,moduleList)
# moduleFunctionComparing("CRNDE",rna2module,moduleList)
}
