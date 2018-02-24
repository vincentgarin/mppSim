###############
# TPRmodelQTL #
###############

#' TPR for model x QTL combination
#'
#' Computes and organises the TPR results for each combination of model and type
#' of QTL. The TPR is averaged over the repetition and MPP design. The function
#' returns the percentage of time a specific type
#' of QTL effect was detected by a specific model over a number of repetitions
#' and MPP design sampling.
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
#' @return A matrix with the precentage of detection over the whole simulation
#' for each type of QTL and model combination.
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

TPRmodelQTL <- function(QTL_true, QTL_detected, d_QTL,
                              n_des = 9, n_mod = 4, n_QTL = 8,
                              mod_names = c("cr", "par", "anc", "biall")){

  n_rep <- length(QTL_detected)

  MPP_names <- paste0("MPP_", 1:n_des)

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

  # sum over the Replication and MPP design

  glb_res <- matrix(0, n_mod, n_QTL)

  for(y in 1:n_rep){

    TPR_res_y <- TPR_res[[y]]

    glb_res_y <- matrix(0, n_mod, n_QTL)

    for(i in 1:length(TPR_res_y)){

      glb_res_y <- glb_res_y + TPR_res_y[[i]]

    }

    glb_res <- glb_res + glb_res_y

  }

  N_tot <- n_rep * n_des

  av_res <- (glb_res/N_tot) * 100

  # sd_res <- apply(X = glb_res, MARGIN = c(1, 2),
  #                 FUN = function(x, y) sd(c(rep(1, x), rep(0, (y-x) ))),
  #                 y = N_tot)
  #
  # sd_res2 <- apply(X = av_res, MARGIN = c(1, 2),
  #                  FUN = function(x) sqrt( x * (1-x) ) )

  return(av_res)

}

