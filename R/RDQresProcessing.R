####################
# RDQresProcessing #
####################

#' Processing RDQ results
#'
#' Function to organise the resolution distance to the QTL (RDQ) results obtained
#' from a list of "true" QTL positions and a list of detected QTLs in several
#' MPP designs with different models. The function compute the average RDQ in
#' each situation (Rep * MPP design * model). Then it return a table with the
#' average RDQ over the x replicate for each combination of MPP design and QTL
#' detection model. It return a similar table with the standard deviation.
#' Finally, it also return the global average per model and per MPP design.
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
#' @param MPP_names MPP design names.
#'
#' @return list with vector of RDQ average per MPP design, and vector of RDQ
#' average per model, matrix of average RDQ per MPP and model, sandard deviation
#' of the RDQ per MPP and model, vector of RDQ average per MPP design, and
#' vector of RDQ average per model.
#'
#' @author Vincent Garin
#'
#' @export

# QTL_true = Q_sel[1:20]
# QTL_detected = Q_res
# MPP_names = MPP_names
# d_QTL = 10
# n_des = 9
# n_mod = 4

RDQresProcessing <- function(QTL_true, QTL_detected, d_QTL, n_des = 9, n_mod = 4,
                             mod_names = c("cr", "par", "anc", "biall"),
                             MPP_names){

  # list to store the results per QTLs

  n_rep <- length(QTL_detected)

  RDQ_res <- vector(mode = "list", length = n_rep)

  for(r in 1:n_rep){

    QTL_true_r <- QTL_true[[r]] # fix the true QTL

    QTL_det_r <- QTL_detected[[r]]

    # create a list to store the results for each MPP design

    RDQ_MPP <- vector(mode = "list", length = n_des)

    for(i in 1:n_des){

      QTL_det_i <- QTL_det_r[[i]]

      # iterate over the different models

      RDQ_i <- rep(0, n_mod)

      for(j in 1:n_mod){

        if(!is.data.frame(QTL_det_i[[j]])){

          if(QTL_det_i[[j]] == "no_QTL"){ RDQ_i[j] <- NA

          } else if (QTL_det_i[[j]] == "error") { RDQ_i[j] <- NA

          }

        } else {

          dist_j <- RDQ(QTL.true = QTL_true_r, QTL.detected = QTL_det_i[[j]],
                        distance = d_QTL)

          RDQ_i[j] <- mean(dist_j, na.rm = TRUE)

        }

        names(RDQ_i) <- mod_names

      }

      RDQ_MPP[[i]] <- RDQ_i

    }

    names(RDQ_MPP) <- MPP_names

    RDQ_res[[r]] <- RDQ_MPP

  }

  # Now average and standard deviation over x reps

  RDQ_res_mat <- vector(mode = "list", length = n_des)

  for(i in 1:n_des){

    # select all the results from MPPi

    MPP_i_res <- lapply(X = RDQ_res, FUN = function(x) x[[i]])

    MPP_i_res <- matrix(unlist(MPP_i_res), ncol = n_mod, byrow = TRUE)

    colnames(MPP_i_res) <- mod_names
    rownames(MPP_i_res) <- paste0("Rep_", 1:n_rep)


    RDQ_res_mat[[i]] <- MPP_i_res

  }

  # produce the average per model

  RDQ_res_av <- lapply(X = RDQ_res_mat, FUN = function(x) colMeans(x, na.rm = TRUE))

  RDQ_res_av <- matrix(unlist(RDQ_res_av), ncol = n_mod, byrow = TRUE)

  colnames(RDQ_res_av) <- mod_names
  rownames(RDQ_res_av) <- MPP_names

  # produce the standard deviation per model

  RDQ_res_sd <- lapply(X = RDQ_res_mat,
                       FUN = function(x) apply(X = x, MARGIN = 2,
                                               FUN = function(x) sd(x = x, na.rm = TRUE)))

  RDQ_res_sd <- matrix(unlist(RDQ_res_sd), ncol = n_mod, byrow = TRUE)

  colnames(RDQ_res_sd) <- mod_names
  rownames(RDQ_res_sd) <- MPP_names

  # average power per model over all design

  RDQ_mod_av <- colMeans(RDQ_res_av)

  # average RDQ per design over models.

  RDQ_des_av <- rowMeans(RDQ_res_av)

  res <- list(av_des = RDQ_des_av, av_mod = RDQ_mod_av, av_des_mod = RDQ_res_av,
              sd_des_mod = RDQ_res_sd)

}

