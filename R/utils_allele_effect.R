#' allele_effect
#'
#' @description Calculate allelic effect given formula from Mollandin et al. 2022
#'
#' @return marker additive effect
#'
#' @noRd


allele_effect <- function(p = 0.5, Sg, k, f = 0.5){

  q <- 1-p
  u <- sample(c(1, -1), size = 1)
  beta <- u * f *  sqrt((k * Sg)/(2*p*q))
  return(beta)

}
