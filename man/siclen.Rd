% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SICLEN.R
\name{siclen}
\alias{siclen}
\title{SICLEN: Single-cell Clustering based on effective noise reduction through ensenble similarity network}
\usage{
siclen(data, nens = 20, knum = 30, alpha = 0.7, npc = 10)
}
\arguments{
\item{data}{(M x N) dimensional matrix for gene expressions
M: The number of genes
N: The number of cells}

\item{nens}{The number of similarity estimations}

\item{knum}{The number of nearest neighbors for KNN networks}

\item{alpha}{Restart probability of the random walker}

\item{npc}{The number of principal components}
}
\description{
SICLEN: Single-cell Clustering based on effective noise reduction through ensenble similarity network
}
\examples{
library(SICLEN)
data(scdata)

logsc <- log2(1+cpm(scdata$counts))
cl_res <- SICLEN(data = logsc)
}
\references{
Accurate single-cell clustering through effective noise reduction over ensemble similarity network
}
\author{
Hyundoo Jeong <hdj@inu.ac.kr>
}
