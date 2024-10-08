#' mc_const_inter
#'
#' @description Calcule le terme d'interaction entre deux croisements
#'
#' @return nombre : le terme d'interaction
#'
#' @noRd

mc_const_inter <- function(ind_i, ind_j){

  C_i_j <- 2 * ind_i$q_i * ind_j$q_i * ind_i$phi_i * ind_j$phi_i
  return(C_i_j)
}
