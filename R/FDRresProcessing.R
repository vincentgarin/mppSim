####################
# FDRresProcessing #
####################

#' Processing FDR results
#'
#' Function to organise the FDR results obtained from a list of "true" QTL
#' positions and a list of detected QTLs in several MPP designs with different
#' models. The function compute the FDR in each situation (Rep * MPP design *
#' model). Then it return a table with the average FDR over the x replicate for
#' each combination of MPP design and QTL detection model. It return a similar
#' table with the standard deviation. Finally, it also return the global average
#' per model and per MPP design.
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
#' @return list with vector of FDR average per MPP design, and vector of FDR
#' average per model, matrix of average FDR per MPP and model, sandard deviation
#' of the FDR per MPP and model, vector of FDR average per MPP design, and
#' vector of FDR average per model.
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

FDRresProcessing <- function(QTL_true, QTL_detected, d_QTL, n_des = 9, n_mod = 4,
                             mod_names = c("cr", "par", "anc", "biall"),
                             MPP_names){

  # list to store the results per QTLs

  n_rep <- length(QTL_detected)

  FDR_res <- vector(mode = "list", length = n_rep)

  for(r in 1:n_rep){

    QTL_true_r <- QTL_true[[r]] # fix the true QTL

    QTL_det_r <- QTL_detected[[r]]

    # create a list to store the results for each MPP design

    FDR_MPP <- vector(mode = "list", length = n_des)

    for(i in 1:n_des){

      QTL_det_i <- QTL_det_r[[i]]

      # iterate over the different models

      FDR_i <- rep(0, n_mod)

      for(j in 1:n_mod){

        if(!is.data.frame(QTL_det_i[[j]])){

          if(QTL_det_i[[j]] == "no_QTL"){ FDR_i[j] <- NA

          } else if (QTL_det_i[[j]] == "error") { FDR_i[j] <- NA

          }

        } else {

          FDR_i[j] <- FDR(QTL.true = QTL_true_r, QTL.detected = QTL_det_i[[j]],
                          distance = d_QTL)

        }

        names(FDR_i) <- mod_names

      }

      FDR_MPP[[i]] <- 100 * FDR_i

    }

    names(FDR_MPP) <- MPP_names

    FDR_res[[r]] <- FDR_MPP

  }

  # Now average and standard deviation over x reps

  FDR_res_mat <- vector(mode = "list", length = n_des)

  for(i in 1:n_des){

    # select all the results from MPPi

    MPP_i_res <- lapply(X = FDR_res, FUN = function(x) x[[i]])

    MPP_i_res <- matrix(unlist(MPP_i_res), ncol = n_mod, byrow = TRUE)

    colnames(MPP_i_res) <- mod_names
    rownames(MPP_i_res) <- paste0("Rep_", 1:n_rep)


    FDR_res_mat[[i]] <- MPP_i_res

  }

  # produce the average per model

  FDR_res_av <- lapply(X = FDR_res_mat, FUN = function(x) colMeans(x, na.rm = TRUE))

  FDR_res_av <- matrix(unlist(FDR_res_av), ncol = n_mod, byrow = TRUE)

  colnames(FDR_res_av) <- mod_names
  rownames(FDR_res_av) <- MPP_names

  # produce the standard deviation per model

  FDR_res_sd <- lapply(X = FDR_res_mat,
                       FUN = function(x) apply(X = x, MARGIN = 2,
                                               FUN = function(x) sd(x = x, na.rm = TRUE)))

  FDR_res_sd <- matrix(unlist(FDR_res_sd), ncol = n_mod, byrow = TRUE)

  colnames(FDR_res_sd) <- mod_names
  rownames(FDR_res_sd) <- MPP_names

  # average power per model over all design

  FDR_mod_av <- colMeans(FDR_res_av)

  # standard deviation over the design

  sd_des <- lapply(X = FDR_res_mat, FUN = function(x) sd(c(x), na.rm = TRUE))
  sd_des <- unlist(sd_des)
  names(sd_des) <- MPP_names

  # average FDR per design over models.

  FDR_des_av <- rowMeans(FDR_res_av)

  # standard deviation TPR per model over the design

  sd_mod <- rep(0, n_mod)

  for(i in 1:n_mod){

    sd_mod[i] <- sd(unlist(lapply(X = FDR_res_mat,
                                  FUN = function(x, i) x[, i], i=i)), na.rm = TRUE)

  }

  res <- list(av_des = FDR_des_av, sd_des = sd_des, av_mod = FDR_mod_av,
              sd_mod = sd_mod,av_des_mod = FDR_res_av, sd_des_mod = FDR_res_sd)

}

