###############
# TPRTable #
###############

#' TPR table
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
#' @param n_QTL Number of QTL. Default = 2.
#'
#' @param mod_names Models names. Default = c("cr", "par", "anc", "biall").
#'
#' @param des_names MPP design names.
#'
#' @return  Average TPR per MPP design and QTL detection model.
#'
#' @author Vincent Garin
#'
#' @export


TPRTable <- function(QTL_true, QTL_detected, d_QTL,
                         n_des = 9, n_mod = 4, n_QTL = 2,
                         mod_names = c("cr", "par", "anc", "biall"),
                         des_names = NULL, n_par, N){

  n_rep <- length(QTL_detected)

  MPP_names <- des_names

  TPR_res <- vector(mode = "list", length = n_rep)

  for(r in 1:n_rep){

    QTL_true_r <- QTL_true[[r]] # fix the true QTL

    QTL_det_r <- QTL_detected[[r]] # select the detected QTLs

    # create a list to store the results for each MPP design

    TPR_MPP <- vector(mode = "list", length = n_des)

    for(i in 1:n_des){

      QTL_det_i <- QTL_det_r[[i]]

      TPR_i <- rep(0, n_mod)

      for(k in 1:n_mod){

        if(!is.data.frame(QTL_det_i[[k]])){

          if(QTL_det_i[[k]] == "no_QTL"){ TPR_i[k] <- 0

          } else if (QTL_det_i[[k]] == "error") { TPR_i[k] <- NA

          }

        } else {

          TPR_i[k] <- TPR(QTL.true = QTL_true_r, QTL.detected = QTL_det_i[[k]],
                          distance = d_QTL)

        }


      }

      names(TPR_i) <- mod_names

      TPR_MPP[[i]] <- TPR_i

    }

    names(TPR_MPP) <- MPP_names

    TPR_res[[r]] <- TPR_MPP

  }

  # sum over the Replication

  MPP_res <- vector(mode = "list", length = n_des)

  for(y in 1:n_des){

    # combine in a single list all results of MPP_i

    MPP_res_i <- lapply(X = TPR_res, FUN = function(x) x[[y]])

    res_0 <- matrix(0, 1, n_mod)

    for(z in 1:n_rep){

      res_0 <- res_0 + MPP_res_i[[z]]

    }

    MPP_res[[y]] <- round(res_0/n_rep * 100, 1)

  }

  # Unfold the results

  res <- matrix(unlist(MPP_res), 4, n_des)
  rownames(res) <- mod_names
  colnames(res) <- des_names

  return(res)

}
