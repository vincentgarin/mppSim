% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dQTLTab.R
\name{dQTLTab}
\alias{dQTLTab}
\title{Table with distance to the simulated QTL}
\usage{
dQTLTab(
  QTL_true = Q_sim[1:n_rep],
  QTL_detected = Q_res,
  des_names = des_names,
  n_par = n_par,
  n_des = n_des,
  n_mod = 4,
  n_QTL = 8,
  mod_names = c("cr", "par", "anc", "biall"),
  N = N
)
}
\arguments{
\item{QTL_true}{list of true QTL positions.}

\item{QTL_detected}{list results of detected QTLs.}

\item{des_names}{MPP design names. Use only: 'Diallel', 'Chessboard',
'Factorial', 'NAM' or 'Realized'. Default = NULL.}

\item{n_par}{Number of parents of the different designs}

\item{n_des}{Number of MPP design. Default = 9.}

\item{n_mod}{Number of QTL detection models. Default = 4.}

\item{n_QTL}{Number of simulated QTL. Default = 8.}

\item{mod_names}{Models names. Default = c("cr", "par", "anc", "biall").}

\item{N}{Total number of individuals in the population}
}
\value{
A matrix with the first column being the average RDQ over the
repetitions and the other being the incidence matrices of the factors
}
\description{
This function measure the distance between the simulated QTL and the largest
significant peak.
}
\author{
Vincent Garin
}
