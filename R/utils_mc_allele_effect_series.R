#' mc_allele_effect_series
#'
#' @description Calculate allelic effect for multiple crosses
#' with a effect following an allelic serie with the parameter s
#'
#' @return list with the additive allele effect
#'  for each cross
#'
#' @noRd


mc_allele_effect_series <- function(k, nc, sigma_G, denom,
                                    s = 0.9^2,
                                    f = 0.5){

  u <- sample(c(1, -1), size = 1)

  # l'effet 'initial'
  beta1_i <- u * f * sqrt((k* sigma_G)/ denom)

  # la liste des facteurs s selon le nombre de croisement :
  list_fact_s<- vector(mode = "list", length = nc)
  list_fact_s[[1]]<-1
  for (r in 2:nc){
    facteur_s <-s^(r-1)
    list_fact_s[[r]] <- facteur_s
  }

  # mélanger les effets pour qu'ils soient alloués
  # aléatoirement à 1 croisement
  facteur_s_bis <- sample(list_fact_s)

  beta_effets<- vector(mode = "list", length = nc)
  # les effets de la série allélique (dépendent de l'effet initial)
  for (i in 1:nc){
    beta_effets[[i]]<- beta1_i * facteur_s_bis[[i]]
  }

  return(beta_effets)
}

