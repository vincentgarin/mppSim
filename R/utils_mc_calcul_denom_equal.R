#' mc_calcul_denom_equal
#'
#' @description Calcule le dénominateur de la formule généralisée
#' pour calculer les effets EGAUX de plusieurs croisements
#' (utilise mc_const_inter, et mc_param_cross)
#'
#' @return nombre : le dénominateur
#'
#' @noRd


mc_calcul_denom_equal <- function(nc, tab_param_cross){

  # partie 1 :
  somme_p1 <- tab_param_cross[[1]]$C_i
  for (i in 2:nc){
    # print(paste("somme = ", somme_p1))
    somme_p1 <- somme_p1 + tab_param_cross[[i]]$C_i
    # print(paste("somme = ", somme_p1))
  }

  # partie 2 :
  somme_p2 <- 0
  for (i in 1:nc){
    # print (paste("i =", i))
    for (j in (i+1):nc){
      # print(paste("j=", j))
      if ( (i!= j) & (j<=length(tab_param_cross))){
        # print("coucou")
        # print(paste("somme = ", somme_p2))
        inter <- mc_const_inter(tab_param_cross[[i]], tab_param_cross[[j]])
        somme_p2 <- somme_p2 + inter
        # print(paste("somme = ", somme_p2))
      }
    }
  }

  denom <- somme_p1 - somme_p2
  return(denom)

}
