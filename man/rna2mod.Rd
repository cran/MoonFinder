\name{rna2mod}
\alias{rna2mod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
rna2mod
}
\description{
Construct a 0-1 binary matrix to represent the association between 
RNAs and protein modules.
}
\usage{
rna2mod(rna2prot, modlist, pCutoff = 0.01, bgProtNum)
}
\arguments{
  \item{rna2prot}{two-column matrix with each row shows the binary 
      association between RNA and protein.}
  \item{modlist}{a list of protein modules.}
  \item{pCutoff}{threshold of the P-value for enrichment analysis.}
  \item{bgProtNum}{the number of background proteins.}
}
\value{
  \item{rna2module}{N x M 0-1 binary matrix of N RNA and M module.}
  \item{moduleList}{a list of modules targeted by at least one RNA.}
}
\references{
MoonFinder: a framework for the identification of moonlighting non-coding RNAs.
}
\author{
Lixin Cheng <easonlcheng@gmail.com>
}
\examples{
data(rna2protein)
data(combinedModuleList)
data(uniGene)
uniGeneNum = length(uniGene)
result = rna2mod(rna2protein[1:100,],combinedModuleList,
    pCutoff = 0.01,uniGeneNum)
rna2module = result[[1]]
moduleList = result[[2]]
}
