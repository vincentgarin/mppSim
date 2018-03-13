###################
# FDR alternative #
###################

#' False discovery rate 2
#'
#' Alternative way to measure the false discovery rate looking at the results on
#' the chromosome where no signal was simulated.
#'
#' @param QTL_true list of true QTL positions.
#'
#' @param QTL_detected list results of detected QTLs.
#'
#' @param n_des Number of MPP design. Default = 9.
#'
#' @param n_mod Number of QTL detection models. Default = 4.
#'
#' @param n_chr Number of chromosome.
#'
#' @return FDR in percent
#'
#' @export

# QTL_true = Q_sim[1:30]
# QTL_detected = Q_res
# n_des = 9
# n_mod = 4
# n_chr = 9

FDR2 <- function(QTL_true, QTL_detected, n_des = 9, n_mod = 4, n_chr){

  chr_ind <- 1:n_chr

  n_rep <- length(QTL_detected)

  N <- n_rep * n_des * n_mod

  TPR_res <- rep(0, N)

  ind.res <- 1

  for(i in 1:n_rep){

    emp_chr <- chr_ind[!(1:n_chr %in% Q_sim[[i]]$chr)] # determine the empty chr

    QTL_det_i <- QTL_detected[[i]]

    for(j in 1:n_des){

      QTL_det_ij <- QTL_det_i[[j]]

      for(k in 1:n_mod){

        if(is.data.frame(QTL_det_ij[[k]])){

          TPR_res[ind.res] <- emp_chr %in% QTL_det_ij[[k]]$chr * 1
          ind.res <- ind.res + 1

        }



      }

    }

  }

  FDR <- sum(TPR_res)/N * 100

  return(FDR)

}
