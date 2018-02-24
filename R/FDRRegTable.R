###############
# FDRRegTable #
###############

#' FDR regression table results
#'
#'
#' @param QTL_true list of true QTL positions.
#'
#' @param QTL_detected list results of detected QTLs.
#'
#' @param d_QTL distance to the QTL.
#'
#' @param n_des Number of MPP design. Default = 9.
#'
#' @param n_mod Number of QTL detection models. Default = 4.
#'
#' @param mod_names Models names. Default = c("cr", "par", "anc", "biall").
#'
#' @param des_names MPP design names. Use only: 'Diallel', 'Chessboard',
#' 'Factorial', 'NAM' or 'Realized'. Default = NULL.
#'
#' @param n_par Number of parents of the different designs.
#'
#' @param N Total number of individuals in the population. Default = 800.
#'
#' @return ...
#'
#' @author Vincent Garin
#'
#' @export

# QTL_true = Q_sim[1:30]
# QTL_detected = Q_res
# MPP_names = MPP_names
# d_QTL = 10
# n_des = 9
# n_mod = 4
# des_names <- c("Diallel", "Diallel", "Chessboard", "Chessboard",
#                "Factorial", "Factorial", "NAM", "NAM", "Realized")
# n_par <- c(9, 5, 9, 5, 9, 5, 9, 5, 9)
# N <- 800

FDRRegTable <- function(QTL_true, QTL_detected, d_QTL, n_des = 9, n_mod = 4,
                        mod_names = c("cr", "par", "anc", "biall"),
                        des_names, n_par, N=800){

  # list to store the results per QTLs

  n_rep <- length(QTL_detected)

  MPP_names <- paste(des_names, n_par, sep = "_")

  des_var <- rep(des_names, each = n_mod)
  Npar_var <- rep(n_par, each = n_mod)
  NindCr_var <- mapply(FUN = function(x, y, N) NindCr(N = N, Npar = x, MPP_des = y),
                       x = Npar_var, y = des_var, MoreArgs = list(N = N))
  mod_var <- rep(mod_names, n_des)

  ResTable <- c()

  for(r in 1:n_rep){

    QTL_true_r <- QTL_true[[r]] # fix the true QTL

    QTL_det_r <- QTL_detected[[r]]

    # create a list to store the results for each MPP design

    FDR_MPP <- matrix(0, nrow = n_mod, ncol = n_des)
    TPR_MPP <- matrix(0, nrow = n_mod, ncol = n_des)

    for(i in 1:n_des){

      QTL_det_i <- QTL_det_r[[i]]

      # iterate over the different models

      FDR_i <- rep(0, n_mod)
      TPR_i <- rep(0, n_mod)

      for(j in 1:n_mod){

        if(!is.data.frame(QTL_det_i[[j]])){

          if(QTL_det_i[[j]] == "no_QTL"){ FDR_i[j] <- NA
          TPR_i[j] <- 0

          } else if (QTL_det_i[[j]] == "error") { FDR_i[j] <- NA
          TPR_i[j] <- NA

          }

        } else {

          FDR_i[j] <- FDR(QTL.true = QTL_true_r, QTL.detected = QTL_det_i[[j]],
                          distance = d_QTL)

          TPR_i[j] <- TPR(QTL.true = QTL_true_r, QTL.detected = QTL_det_i[[j]],
                          distance = d_QTL)

        }

      }

      FDR_MPP[, i] <- FDR_i * 100
      TPR_MPP[, i] <- TPR_i * 100

      colnames(FDR_MPP) <- MPP_names
      rownames(FDR_MPP) <- mod_names

      colnames(TPR_MPP) <- MPP_names
      rownames(TPR_MPP) <- mod_names

    }

    data_r <- data.frame(c(FDR_MPP), c(TPR_MPP), des_var, Npar_var, NindCr_var,
                         mod_var, stringsAsFactors = FALSE)

    ResTable <- rbind(ResTable, data_r)

  }


  colnames(ResTable) <- c("FDR", "TPR",  "MPP_des", "N_par", "N_ind_cr", "Model")

  TPR_FDR <- ResTable$TPR/ResTable$FDR
  max_sc <- max(TPR_FDR[!is.infinite(TPR_FDR)])
  TPR_FDR[is.infinite(TPR_FDR)] <- max_sc

  ResTable <- data.frame(ResTable, TPR_FDR)

  ResTable$Model <- factor(ResTable$Model, levels = c("cr", "par", "anc", "biall"))
  ResTable$MPP_des <- factor(ResTable$MPP_des, levels = c("Diallel", "Chessboard",
                                                          "Factorial", "NAM", "Realized"))

  return(ResTable)


}

