% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FDR.R
\name{FDR}
\alias{FDR}
\title{False discovery rate}
\usage{
FDR(QTL.true, QTL.detected, distance = 5)
}
\arguments{
\item{QTL.true}{\code{data.frame} of true QTL positions with marker id,
chromosome and position in cM.}

\item{QTL.detected}{\code{data.frame} of detected QTL positions with marker
id, chromosome and position in cM.}

\item{distance}{maximal distance between the true QTL position and the detected
signal. Default = 5.}
}
\value{
FDR
}
\description{
Function computing the false discovery rate (FDR) as the percentage of detected
QTL positions that are not in a neighbouring region (+-distance) from a true
QTL position.
}
\author{
Vincent Garin
}
