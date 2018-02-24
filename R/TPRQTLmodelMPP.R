##################
# TPRQTLmodelMPP #
##################

#' TPR for QTL x model x MPP combination
#'
#' Computes and organises the TPR results for each combination of QTL, model,
#' and MPP design. The TPR is averaged over the repetitions. The function
#' returns the average and the standard deviation of the TPR for a specific
#' combination of QTL model and MPP design.
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
#' @param n_QTL Number of QTL. Default = 8.
#'
#' @param mod_names Models names. Default = c("cr", "par", "anc", "biall").
#'
#' @param MPP_names MPP design names.
#'
#' @return A matrix with the precentage of detection over the repetition of the
#' simulation procedure.
#'
#' @author Vincent Garin
#'
#' @export

# QTL_true = Q_sel[1:20]
# QTL_detected = Q_res
# d_QTL = 10
# n_des = 9
# n_mod = 4
# n_QTL = 8
# mod_names = c("cr", "par", "anc", "biall")
# MPP_names <- c("NAM_9", "NAM_5", "Dia_9", "Dia_5", "Chess_9", "Chess_5",
#                "Fact_9", "Fact_5", "Real_9")


TPRQTLmodelMPP <- function(QTL_true, QTL_detected, d_QTL,
                        n_des = 9, n_mod = 4, n_QTL = 8,
                        mod_names = c("cr", "par", "anc", "biall"),
                        MPP_names = NULL){

  n_rep <- length(QTL_detected)

  if(is.null(MPP_names)){ MPP_names <- paste0("MPP_", 1:n_des)}

  TPR_res <- vector(mode = "list", length = n_rep)

  for(r in 1:n_rep){

    QTL_true_r <- QTL_true[[r]] # fix the true QTL

    QTL_det_r <- QTL_detected[[r]]

    # create a list to store the results for each MPP design

    TPR_MPP <- vector(mode = "list", length = n_des)

    for(i in 1:n_des){

      QTL_det_i <- QTL_det_r[[i]]

      res <- matrix(0, n_mod, n_QTL)

      for(j in 1:n_QTL){

        QTL_true_rj <- QTL_true_r[QTL_true_r$Qeff == j, ]

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

        res[, j] <- TPR_j

      }

      rownames(res) <- mod_names
      colnames(res) <- paste0("Q", 1:n_QTL)

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

    res_0 <- matrix(0, n_mod, n_QTL)

    for(z in 1:n_rep){

      res_0 <- res_0 + MPP_res_i[[z]]

    }

    MPP_res[[y]] <- res_0

  }

  # Organise per QTL type

  QTL_res <- vector(mode = "list", length = n_QTL)

  for(l in 1:n_QTL){

  res_Q_i <- matrix(0, n_mod, n_des)

  for(m in 1:n_des){

    res_Q_i[, m] <- MPP_res[[m]][, l]

  }

  rownames(res_Q_i) <- mod_names
  colnames(res_Q_i) <- MPP_names

  QTL_res[[l]] <- (res_Q_i/n_rep) * 100

  }

  names(QTL_res) <- paste0("Q", 1:n_QTL)

  return(QTL_res)

}

