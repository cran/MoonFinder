\name{SDI}
\alias{SDI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Shannon diversity index
}
\description{
Calculate the Simpson's diversity index (SDI).
}
\usage{
SDI(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{weights of latent features.}
}
\value{
  \item{z}{moonlighting coefficient, or normalized Shannon's diversity index.}
}
\references{
MoonFinder: a framework for the identification of moonlighting non-coding RNAs.
}
\author{
Lixin Cheng <easonlcheng@gmail.com>
}
\examples{
a = sample(1:10,4)
b = a/sum(a)
SDI(b)
}
