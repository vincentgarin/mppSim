#' mc_calcul_pheno
#'
#' @description Calculate phenotypes with allelic effect given for
#' multiple crosses
#'
#' @return vector with phenotypes in one column
#' (rownames remains unchanged)
#'
#' @noRd
#'

# nb_mk = n_mk_f
# tab_res <- tab_res_f
# nc = 10
# geno = X_f
# geno_i<-geno[cr_ind %in% cr_id[i], ]

mc_calcul_pheno <- function(nc, nb_mk, tab_res){

  # Récuperer les betas pour chaque croisement
  # et chaque marqueur (en 1 colonne, avec
  # par ex si y'a 5 mk et 10 cross, d'abord les 5
  # betas du cross 1 puis les 5 beta du cross 2 etc...)

  tab_beta <- matrix(rep(NA, (nb_mk*nc)))
  compteur <-1
  for (j in 1:nc){

    for (i in 1:length(tab_res)){

      tab_beta[compteur] <- tab_res[[i]]$beta[[j]]
      compteur <- compteur +1
      }
  }

  pheno_i <-vector(mode = 'list', length = nc)
  compteur2 <-0
  for (i in 1:nc){
    tab_beta_i <- tab_beta[((compteur2)+1): (compteur2 +length(tab_beta)/(nc)),
                           drop = FALSE]

    X_cross <- tab_res[[1]]$tab_param_cross[[i]]$geno_i
    # le 1 c'est pour le 1er mk ( mais c'est
    # la même matrice pour les autres mk)
    # le i c'est pour le numero du cross
    # c'est la matrice
    # de genotypage mais avec uniquement les individus du cross
    # en question

    pheno_i[[i]] <- X_cross %*% matrix(tab_beta_i)
    compteur2 <- compteur2 + length(tab_beta)/(nc)
    # print(compteur2)
  }


  pheno <-pheno_i[[1]]
  for (i in 2:nc){
    pheno <-rbind(pheno, pheno_i[[i]])
  }

  colnames(pheno)<- c("pheno")


  return(list(pheno = pheno, beta = tab_beta))
}
