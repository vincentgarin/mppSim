##########
# NindCr #
##########

#' Nb individual per cross
#'
#' Function to determine the average number of individual per cross given the
#' total size of the population, the type of MPP design and the number of parents
#'
#' @param N Total number of individual in the population
#'
#' @param Npar Number of parents
#'
#' @param MPP_des String indicating the type of MPP design among: 'Diallel',
#' 'Chessboard', 'Factorial', 'NAM', 'Realized'
#'
#' @return Number of individuals per cross
#'
#' @author Vincent Garin
#'

NindCr <- function(N, Npar, MPP_des){

  if(MPP_des == "Diallel"){

    ncr <- (Npar * (Npar - 1))/2

  } else if (MPP_des %in% c("Chessboard", "Factorial")){

    np_mod <- Npar %/% 2

    rest <- Npar - np_mod

    ncr <- np_mod * rest


  }  else if (MPP_des == "NAM") {

    ncr <- Npar - 1

  } else if (MPP_des == "Realized"){

    ncr <- 11

  }

  N/ncr

}
