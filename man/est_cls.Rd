% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SICLEN.R
\name{est_cls}
\alias{est_cls}
\title{SICLEN: Single-cell Clustering based on effective noise reduction through ensenble similarity network}
\usage{
est_cls(data)
}
\arguments{
\item{data}{(M x N) dimensional matrix for gene expressions
M: The number of genes
N: The number of cells}
}
\value{
Estimated number of clusters in the single-cell RNA sequencing data
}
\description{
This function estimates the number of clusters in the single-cell RNA sequencing data
}
\references{
Accurate single-cell clustering through effective noise reduction over ensemble similarity network
}
\author{
Hyundoo Jeong <hdj@inu.ac.kr>
}
