\name{moduleCombiding}
\alias{moduleCombiding}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Combine modules from different resources.
}
\description{
Combine modules with high overlapping score.
}
\usage{
moduleCombiding(inputList, OScutoff = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputList}{a list of protein modules.}
  \item{OScutoff}{the threshold of overlapping score.}
}
\details{
Combine modules from different resources. For example, one set is computational
predicted while the other sets are collected from databases.
}
\value{
  \item{outputList}{a list of combined modules.}
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
m4 = c("Q92828", "Q13227", "O15379", "O75376", "Q15021")
m5= c("P61160", "P61158", "O15143", "O15144", "O15145", "Q9NTJ3")
mList1 = list(m1,m2,m3,m4,m5)
mList2 = moduleCombiding(mList1,OScutoff = 0.5)
head(mList2)
}
