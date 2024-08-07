% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TPRRegTableIndQTL.R
\name{TPRRegTableIndQTL}
\alias{TPRRegTableIndQTL}
\title{Table of TPR results for single QTLs}
\usage{
TPRRegTableIndQTL(
  QTL_true,
  QTL_detected,
  d_QTL,
  n_des = 9,
  n_mod = 4,
  n_QTL = 8,
  mod_names = c("cr", "par", "anc", "biall"),
  des_names = NULL,
  n_par,
  N
)
}
\arguments{
\item{QTL_true}{list of true QTL positions with frequency and number of effect
information added.}

\item{QTL_detected}{list results of detected QTLs.}

\item{d_QTL}{distance to the QTL.}

\item{n_des}{Number of MPP design. Default = 9.}

\item{n_mod}{Number of QTL detection models. Default = 4.}

\item{n_QTL}{Number of QTL. Default = 8.}

\item{mod_names}{Models names. Default = c("cr", "par", "anc", "biall").}

\item{des_names}{MPP design names. Use only: 'Diallel', 'Chessboard',
'Factorial', 'NAM' or 'Realized'. Default = NULL.}

\item{n_par}{Number of parents of the different designs.}

\item{N}{Total number of individuals in the population.}
}
\value{
A data.frame with first column indicating if the QTL was detected
or not (0, 1), the type of QTL effect, the size of the QTL effect,
the segregation rate (proportion of ind that receive a non zero effect from
the considered QTL), the number of QTL effects, the model, the MPP design,
the number of parents, the number of individuals per cross.
}
\description{
Organise the TPR results in a table that can be used for regression. The
independent variable (y) is the average TPR decomposed over each factor:
MPP design, number of parents, individual per cross (continuous variable),
QTL detection model, type of QTL, and the size of the QTL.
}
\author{
Vincent Garin
}
