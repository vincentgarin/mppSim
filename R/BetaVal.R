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

      # Beta <- AlleleSeries(na = dim(Qmat)[2], a = 0.9)
      Beta <- round(runif(n = dim(Qmat)[2], min = 1, max = 10), 2)
      Beta_sign <- sample(x = c(-1, 1), size = dim(Qmat)[2], replace = TRUE)

      Beta_list[[i]] <- Beta * Beta_sign

    } else if (Qsc == 2){ # half of the crosses

      na <- dim(Qmat)[2]
      na_half <- round(na/2)

      # Beta <- AlleleSeries(na = na_half, a = 0.8)
      Beta <- round(runif(n = na_half, min = 1, max = 10), 2)
      Beta_sign <- sample(x = c(-1, 1), size = na_half, replace = TRUE)
      Beta <- c(Beta * Beta_sign, rep(0, (na - na_half)))
      Beta_list[[i]] <- Beta[sample(1:na)]

    } else if (Qsc == 3){ # all parents

      # Beta <- AlleleSeries(na = dim(Qmat)[2], a = 0.7)
      Beta <- round(runif(n = dim(Qmat)[2], min = 1, max = 10), 2)
      Beta_sign <- sample(x = c(-1, 1), size = dim(Qmat)[2], replace = TRUE)

      Beta_list[[i]] <- Beta * Beta_sign

    } else if (Qsc == 4){ # a single parents

      na <- dim(Qmat)[2]

      Beta <- round(runif(n = 1, min = 1, max = 10), 2)
      Beta_sign <- sample(x = c(-1, 1), size = 1, replace = TRUE)
      Beta <- c(rep(0, (na - 1)), Beta * Beta_sign)

      Beta_list[[i]] <- Beta[sample(1:na)]

    } else if (Qsc == 5){ # all ancestors

      # Beta <- AlleleSeries(na = dim(Qmat)[2], a = 0.5)

      Beta <- round(runif(n = dim(Qmat)[2], min = 1, max = 10), 2)
      Beta_sign <- sample(x = c(-1, 1), size = dim(Qmat)[2], replace = TRUE)

      Beta_list[[i]] <- Beta * Beta_sign

    } else if (Qsc == 6){ # a single ancestor (largest)

      na <- dim(Qmat)[2]
      Beta <- round(runif(n = 1, min = 1, max = 10), 2)
      Beta_sign <- sample(x = c(-1, 1), size = 1, replace = TRUE)
      Beta <- c(rep(0, (na - 1)), Beta * Beta_sign)

      # Beta_list[[i]] <- Beta[sample(1:na)]
      Beta_list[[i]] <- Beta

    } else if (Qsc == 7){ # bi-allelic minor allele

      Beta <- round(runif(n = 1, min = 1, max = 10), 2)
      Beta_sign <- sample(x = c(-1, 1), size = 1, replace = TRUE)

      # Beta_list[[i]] <- c(-Beta, Beta)
      Beta_list[[i]] <- Beta * Beta_sign

      ############### c(-0.5, 0.5)[sample[1:2]]

    } else if (Qsc == 8){ # bi-allelic minor allele

      Beta <- round(runif(n = 1, min = 1, max = 10), 2)
      Beta_sign <- sample(x = c(-1, 1), size = 1, replace = TRUE)

      # Beta_list[[i]] <- c(-Beta, Beta)
      Beta_list[[i]] <- Beta * Beta_sign

    }

  }

  return(Beta_list)

}
