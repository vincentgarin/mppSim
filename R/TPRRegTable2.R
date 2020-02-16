###############
# TPRRegTable2 #
###############

#' TPR results for single type of QTL effect simulation
#'
#' Organise the TPR results in a table that can be used for regression. The
#' independent variable (y) is the average TPR decomposed over each factor:
#' MPP design, number of parents, individual per cross (continuous variable),
#' QTL detection model, type of QTL, and the size of the QTL.
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
#' @param n_QTL Number of QTL. Default = 2.
#'
#' @param QTL_ind Numeric indicator for the type of QTL between 1 to 8.
#'
#' @param mod_names Models names. Default = c("cr", "par", "anc", "biall").
#'
#' @param des_names MPP design names. Use only: 'Diallel', 'Chessboard',
#' 'Factorial', 'NAM' or 'Realized'. Default = NULL.
#'
#' @param n_par Number of parents of the different designs
#'
#' @param N Total number of individuals in the population
#'
#' @return A matrix with the first column being the average TPR over the
#' repetitions and the other being the incidence matrices of the factors
#'
#' @author Vincent Garin
#'
#' @export


TPRRegTable2 <- function(QTL_true, QTL_detected, d_QTL,
                        n_des = 9, n_mod = 4, n_QTL = 2, QTL_ind,
                        mod_names = c("cr", "par", "anc", "biall"),
                        des_names = NULL, n_par, N){

  n_rep <- length(QTL_detected)

  MPP_names <- paste(des_names, n_par, sep = "_")

  TPR_res <- vector(mode = "list", length = n_rep)

  for(r in 1:n_rep){

    QTL_true_r <- QTL_true[[r]] # fix the true QTL

    # determine the sub-type of QTL: Q1 small, Q1 large, Q2 small, etc.

    QtypeId <- 1:(2*n_QTL)
    QTL_comb <- paste0(rep(QTL_ind, each = 2), rep(c(2, 6), n_QTL))
    names(QtypeId) <- QTL_comb

    # Qtype <- QtypeId[paste0(QTL_true_r$Qeff, QTL_true_r$Qsize)]

    QTL_det_r <- QTL_detected[[r]] # select the detected QTLs

    # create a list to store the results for each MPP design

    TPR_MPP <- vector(mode = "list", length = n_des)

    for(i in 1:n_des){

      QTL_det_i <- QTL_det_r[[i]]

      res <- matrix(0, n_mod, 2*n_QTL)

      # need to select a single combination of QTL effect (type + size)

      QTL_comb_i <- paste0(QTL_true_r$Qeff, QTL_true_r$Qsize)


      for(j in 1:length(QTL_comb)){

        QTL_true_rj <- QTL_true_r[QTL_comb_i %in% QTL_comb[j], ]

        # iterate over the different models

        TPR_j <- rep(0, n_mod)

        for(k in 1:n_mod){

          if(!is.data.frame(QTL_det_i[[k]])){

            if(QTL_det_i[[k]] == "no_QTL"){ TPR_j[k] <- 0

            } else if (QTL_det_i[[k]] == "error") { TPR_j[k] <- NA

            }

          } else {

            TPR_j[k] <- TPR(QTL.true = QTL_true_rj, QTL.detected = QTL_det_i[[k]],
                            distance = d_QTL)

          }


        }

        index <- QtypeId[QTL_comb[j]]

        res[, index] <- TPR_j

      }

      rownames(res) <- mod_names
      colnames(res) <- paste0('Q', QTL_comb)

      TPR_MPP[[i]] <- res

    }

    names(TPR_MPP) <- MPP_names

    TPR_res[[r]] <- TPR_MPP

  }

  # sum over the Replication

  MPP_res <- vector(mode = "list", length = n_des)

  for(y in 1:n_des){

    # combine in a single list all results of MPP_i

    MPP_res_i <- lapply(X = TPR_res, FUN = function(x) x[[y]])

    res_0 <- matrix(0, n_mod, 2*n_QTL)

    for(z in 1:n_rep){

      res_0 <- res_0 + MPP_res_i[[z]]

    }

    MPP_res[[y]] <- res_0/n_rep * 100

  }

  # Unfold the results

  # For the first design

  ResTable <- c()

  n_el <- (n_mod * n_QTL * 2)
  ModVar <- rep(mod_names, n_QTL * 2) #
  Qlabel <- c("Q1_cr", "Q2_cr", "Q3_par", "Q4_par", "Q5_anc", "Q6_anc",
              "Q7_biall", "Q8_biall")

  Qlabel <- Qlabel[QTL_ind]

  QtypeVar <- rep(Qlabel, each = 2*n_mod)
  QsizeVar <- rep(rep(c("small", "big"), each = n_mod), n_QTL)

  for(x in 1:n_des){

    y_x <- c(MPP_res[[x]])

    des_x <- rep(des_names[x], n_el)
    Npar_x <- rep(n_par[x], n_el)

    res_x <- data.frame(y_x, des_x, Npar_x, ModVar, QtypeVar, QsizeVar,
                        stringsAsFactors = FALSE)

    colnames(res_x) <- c("TPR", "MPP_des", "N_par", "Model", "QTL_type",
                         "QTL_size")

    ResTable <- rbind(ResTable, res_x)

  }

  NindCrVar <- mapply(FUN = function(x, y, N) NindCr(N = N, Npar = x, MPP_des = y),
                      x = ResTable$N_par, y = ResTable$MPP_des,
                      MoreArgs = list(N = N))

  ResTable <- data.frame(ResTable, NindCrVar,stringsAsFactors = FALSE)
  colnames(ResTable)[7] <- "N_ind_cr"

  # Factorize to obtain the results in the correct order

  ResTable$Model <- factor(ResTable$Model, levels = c("cr", "par", "anc", "biall"))
  ResTable$MPP_des <- factor(ResTable$MPP_des, levels = c("Diallel", "Chessboard",
                                                          "Factorial", "NAM", "Realized"))
  ResTable$QTL_type <- factor(ResTable$QTL_type, levels = Qlabel)

  ResTable$QTL_size <- factor(ResTable$QTL_size, levels = c("small", "big"))

  return(ResTable)

}
