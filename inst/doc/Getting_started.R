## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(MoonFinder)
data("ppi","simMatBP.Resnik","moduleList","rna2module")
rnaStatMat = moonlightingCoefficient(rna2module,simMatBP.Resnik)
head(rnaStatMat)

## ------------------------------------------------------------------------
sdi = as.numeric(rnaStatMat[,3])
dev.new(width=4, height=4)
hist(sdi,border=FALSE,xlab="MC",xlim=c(0,1),ylim=c(0,4),
    breaks=seq(0,1,0.02),prob=TRUE,main="Resnik")
lines(density(sdi),lwd=2,col="orange")
mlncRNA = rnaStatMat[sdi>0.4,1]

## ----fig.width = 6,fig.height = 6----------------------------------------
library(igraph)
g <- graph.edgelist(ppi,directed=FALSE)
moduleHeatmap("CRNDE",rna2module,moduleList,g)

## ----fig.width = 8-------------------------------------------------------
library(clusterProfiler)
dev.new(width=14, height=7)
moduleFunctionComparing("CRNDE",rna2module,moduleList)

