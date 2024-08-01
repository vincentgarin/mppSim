#' cross_scheme_NAM
#'
#' @description Form crossing scheme for a NAM population with np parents.
#'
#' @param np Parent number in the design.
#'
#' @return parent per cross information
#'
#' @author Vincent Garin
#'
#' @export

cross_scheme_NAM <- function(np){

  n_cross <- np -1

  cr_index <- paste0("%0", nchar(n_cross), "d")
  crosses <- paste0("cr", sprintf(cr_index, 1:n_cross))

  P1 <- rep("P1", n_cross)
  P2 <- paste0("P", 2:np)

  par_per_cross <- cbind(crosses, P1, P2)

  return(par_per_cross)

}
