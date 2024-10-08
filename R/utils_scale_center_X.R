#' scale_center_X
#'
#' @description Scale and center a genotype marker matrix
#'
#' @return Standardized marker matrix
#'
#' @noRd

scale_center_X <- function(X){

  p <- apply(X, 2, FUN = function(x) sum(x)/(2*length(x)))
  P <- rep(1, nrow(X)) %*% (2*t(p))
  Z <- X - P
  den <- sqrt(2 * p * (1 - p))
  Z <- t(t(Z) / den)
  # equivalent to Z %*% diag(1 / den)

  return(Z)

}
