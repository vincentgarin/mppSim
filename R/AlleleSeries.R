################
# AlleleSeries #
################

#' Form QTL incidence matrices
#'
#' The function for the QTL incidence matrices of position given different type
#' of effects: cross-specific, parental, ancestral and bi-allelic.
#'
#' @param na Number of alleles
#'
#' @param a parameter of the function
#'
#' @return vector with values of the allelic series
#'
#' @export


AlleleSeries <- function(na, a){

  x <- 1:na

  # determine if the number of a alleles is even

  na_even <- ((na %% 2) == 0)

  if(na_even){

    g.ser <- a^(1:(na/2))
    all.ser <- c(g.ser, -rev(g.ser))

  } else {

    g.ser <- a^(1:((na - 1)/2))
    all.ser <- c(g.ser, 0, -rev(g.ser))

  }

  return(all.ser)

}
