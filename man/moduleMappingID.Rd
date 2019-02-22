\name{moduleMappingID}
\alias{moduleMappingID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
moduleMappingID
}
\description{
Map protein ID for a list of modules.
}
\usage{
moduleMappingID(inputList, ID1 = "UNIPROT", ID2 = "SYMBOL", thdMS = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputList}{a list of modules with proteins in ID1.}
  \item{ID1}{ID type of proteins of the input moduels.}
  \item{ID2}{ID type of protiens of the output moduels.}
  \item{thdMS}{mimimum module size.}
}
\value{
  \item{outputList}{a list of modules with proteins in ID2.}
}
\references{
MoonFinder: a framework for the identification of moonlighting non-coding RNAs.
}
\author{
Lixin Cheng <easonlcheng@gmail.com>
}
\examples{
m1 = c("Q15021", "Q9BPX3", "Q15003", "O95347", "Q9NTJ3")
m2 = c("Q92828", "Q13227", "O15379", "O75376", "O60907", "Q9BZK7")
m3= c("P61160", "P61158", "O15143", "O15144", "O15145", "P59998", "O15511")
uniprotModuleList = list(m1,m2,m3)
symbolModuleList = moduleMappingID(uniprotModuleList)
names(symbolModuleList) = paste0("M",1:length(symbolModuleList))
head(symbolModuleList)
}
