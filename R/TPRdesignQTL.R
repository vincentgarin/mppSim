################
# TPRdesignQTL #
################

#' TPR for MPP design x QTL combination
#'
#' Computes and organises the TPR results for each combination of MPP design
#' and type of QTL. The TPR is averaged over the repetition and model.
#' The function returns the percentage of time a specific type
#' of QTL effect was detected by in a specific MPP design over a number of
#' repetitions and QTL detection models.
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
#' @param MPP_names MPP design names.
#'
#' @return A matrix with the precentage of detection over the whole simulation
#' for each type of QTL and MPP design combination.
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
# MPP_names = MPP_names

TPRdesignQTL <- function(QTL_true, QTL_detected, d_QTL, n_des = 9, n_mod = 4,
                         n_QTL = 8, MPP_names){

  n_rep <- length(QTL_detected)

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

      res <- colSums(res, na.rm = TRUE)
      names(res) <- paste0("Q", 1:n_QTL)

      TPR_MPP[[i]] <- res

    }

    names(TPR_MPP) <- MPP_names

    TPR_res[[r]] <- TPR_MPP

  }

  # sum over the Replication and MPP design

  glb_res <- matrix(0, n_des, n_QTL)

  TPR_res_mat <- lapply(X = TPR_res,
                        FUN = function(x, n_des) matrix(unlist(x),
                                                        nrow = n_des,
                                                        byrow = TRUE),
                        n_des = n_des)

  for(y in 1:n_rep){

    glb_res <- glb_res + TPR_res_mat[[y]]

  }

  rownames(glb_res) <- MPP_names
  colnames(glb_res) <- paste0("Q", 1:n_QTL)

  N_tot <- n_rep * n_mod

  av_res <- (glb_res/N_tot) * 100


  return(av_res)

}

