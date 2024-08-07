% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TPRdesignQTL.R
\name{TPRdesignQTL}
\alias{TPRdesignQTL}
\title{TPR for MPP design x QTL combination}
\usage{
TPRdesignQTL(
  QTL_true,
  QTL_detected,
  d_QTL,
  n_des = 9,
  n_mod = 4,
  n_QTL = 8,
  MPP_names
)
}
\arguments{
\item{QTL_true}{list of true QTL positions.}

\item{QTL_detected}{list results of detected QTLs.}

\item{d_QTL}{distance to the QTL.}

\item{n_des}{Number of MPP design. Default = 9.}

\item{n_mod}{Number of QTL detection models. Default = 4.}

\item{n_QTL}{Number of QTL. Default = 8.}

\item{MPP_names}{MPP design names.}
}
\value{
A matrix with the precentage of detection over the whole simulation
for each type of QTL and MPP design combination.
}
\description{
Computes and organises the TPR results for each combination of MPP design
and type of QTL. The TPR is averaged over the repetition and model.
The function returns the percentage of time a specific type
of QTL effect was detected by in a specific MPP design over a number of
repetitions and QTL detection models.
}
\author{
Vincent Garin
}
