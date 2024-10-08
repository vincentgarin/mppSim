#' mc_allele_effect_equal
#'
#' @description Calculate allelic effect for multiple crosses
#' with a effect equal for all crosses POUR 1 MARQUEUR
#'
#' @return list with the additive allele effect
#'  for each cross (identical)
#'
#' @noRd


mc_allele_effect_equal <- function(k, nc, sigma_G, denom,
                                   f = 0.5){

  u <- sample(c(1, -1), size = 1)

  # l'effet égal pour tous les cross :
  beta_i <- u * f * sqrt((k* sigma_G)/ denom)

  beta_effets<- vector(mode = "list", length = nc)
  for (i in 1:nc){
    beta_effets[[i]]<- beta_i
  }
  return(beta_effets)
  # les betas en 1 colonne avec un élément de chaque ligne
  # pour le nc correspondant
}
