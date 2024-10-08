#' mc_param_cross
#'
#' @description Récupére les paramètres de chaque croisement
#'
#' @return list avec les paramètres de chaque croisement
#' 1 élément de la liste pour chaque croisement, et dans cet
#' élement, il y a la matrice de genotypage pour ce cross,
#' la matrice pour le marqueur selectionne, q_i, p_i, phi_i et
#' C_i ( la constante calculée suivant notre formule)
#'
#' @noRd


mc_param_cross <- function (geno, pos_i, cr_ind){

  cr_id <-unique(cr_ind)

  nc <- length(cr_id)
  # normalement la matrice geno qui est donnée,
  # c'est parce qu'on veut tous les croisements, donc
  # c'est pour ca que nc est à l'intérieur de la fct
  # et n'est pas un param qu'on demande en entrée

  tab_param_cross <- vector(mode = 'list', length = nc)

  for (i in 1:nc){

    # if (ncol(geno) == 1) {
    #   # pour conserver le format matrice quand il y a un seul QTL
    #   geno_i <- matrix(geno[cr_ind %in% cr_id[i]])
    # } else {
    #   geno_i <- geno[cr_ind %in% cr_id[i], ]
    # }

    # pour conserver le format matrice même quand il y a un seul QTL
    # y' a deux virgules avant drop : [ligne, colonne, drop]
    geno_i <- geno[cr_ind %in% cr_id[i], ,drop = FALSE ]

    # geno_i<-geno[cr_ind %in% cr_id[i], ]

    geno_i_sel <- geno_i[, pos_i, drop = FALSE]

    q_i <- as.numeric(apply(geno_i_sel, MARGIN = 2, FUN = function(x) sum(x)/(2*length(x))))

    p_i <- 1 - q_i

    phi_i <- (nrow(geno)/nc) / nrow(geno)

    # Constante Ci :

    C_i <- (0.5 * p_i * q_i * phi_i) + (q_i^2 * (phi_i * (1 - phi_i)))

    list_i <- list(geno_i = geno_i,
                   geno_i_sel = geno_i_sel,
                   q_i = q_i,
                   p_i = p_i,
                   phi_i = phi_i,
                   C_i = C_i)


    tab_param_cross[[i]]<-list_i


  }

  return(tab_param_cross)

}
