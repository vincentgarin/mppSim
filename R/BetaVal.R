###########
# BetaVal #
###########

# function to form the Beta vectors that correspond to the differents
# scenarios in term of type of QTL effect (cr, par, anc, biall) and
# distribution of the effect.


BetaVal <- function(QTL, QTL_inc){

  Beta_list <- vector(mode = "list", length = length(QTL_inc))

  for(i in 1:length(QTL_inc)){

    Qsc <- QTL[i, 5]
    Qmat <- QTL_inc[[i]]

    if(Qsc == 1){ # all crosses

      Beta <- AlleleSeries(na = dim(Qmat)[2], a = 0.9)
      Beta_list[[i]] <- Beta[sample(1:dim(Qmat)[2])]

    } else if (Qsc == 2){ # half of the crosses

      na <- dim(Qmat)[2]
      na_half <- round(na/2)

      Beta <- AlleleSeries(na = na_half, a = 0.8)
      Beta <- c(Beta, rep(0, (na - na_half)))
      Beta_list[[i]] <- Beta[sample(1:na)]

    } else if (Qsc == 3){ # all parents

      Beta <- AlleleSeries(na = dim(Qmat)[2], a = 0.7)
      Beta_list[[i]] <- Beta[sample(1:dim(Qmat)[2])]

    } else if (Qsc == 4){ # a single parents

      na <- dim(Qmat)[2]
      Beta_list[[i]] <- c(1, rep(0, (na - 1)))[sample(1:na)]

    } else if (Qsc == 5){ # all ancestors

      Beta <- AlleleSeries(na = dim(Qmat)[2], a = 0.5)
      Beta_list[[i]] <- Beta[sample(1:dim(Qmat)[2])]

    } else if (Qsc == 6){ # a single ancestor (largest)

      na <- dim(Qmat)[2]
      Beta_list[[i]] <- c(rep(0, (na - 1)), 1)

    } else if (Qsc == 7){ # bi-allelic minor allele

      Beta_list[[i]] <- c(-0.5, 0.5)

      ############### c(-0.5, 0.5)[sample[1:2]]

    } else if (Qsc == 8){ # bi-allelic minor allele

      Beta_list[[i]] <- c(-0.5, 0.5)

    }

  }

  return(Beta_list)

}
