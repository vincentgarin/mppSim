#' mk_eff_dist
#'
#' @description Calculate the marker effect distribution
#'
#' @return vector of value proportional to the marker additive effects
#'
#' @noRd


mk_eff_dist <- function(n = 20, theta = 0.9, plot = FALSE){
  y <- theta^{1:n}
  if(plot){

    print(plot(x = 1:n, y = y, main = "mk effect decay"))

  }
  return(y)
}
