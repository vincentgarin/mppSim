####################
# TPRresProcessing #
####################

#' Processing TPR results
#'
#' Function to organise the TPR results obtained from a list of "true" QTL
#' positions and a list of detected QTLs in several MPP designs with different
#' models. The function compute the TPR in each situation (Rep * MPP design *
#' model). Then it return a table with the average TPR over the x replicate for
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
#' @return list with vector of TPR average per MPP design, and vector of TPR
#' average per model, matrix of average TPR per MPP and model, sandard deviation
#' of the TPR per MPP and model, vector of TPR average per MPP design, and
#' vector of TPR average per model.
#'
#' @author Vincent Garin
#'
#' @export

# QTL_true = Q_sel[1:20]
# QTL_detected = Q_res
# d_QTL = 10

TPRresProcessing <- function(QTL_true, QTL_detected, d_QTL, n_des = 9, n_mod = 4,
                             mod_names = c("cr", "par", "anc", "biall"),
                             MPP_names){

  # list to store the results per QTLs

  n_rep <- length(QTL_detected)

  TPR_res <- vector(mode = "list", length = n_rep)

  for(r in 1:n_rep){

    QTL_true_r <- QTL_true[[r]] # fix the true QTL

    QTL_det_r <- QTL_detected[[r]]

    # create a list to store the results for each MPP design

    TPR_MPP <- vector(mode = "list", length = n_des)

    for(i in 1:n_des){

      QTL_det_i <- QTL_det_r[[i]]

      # iterate over the different models

      TPR_i <- rep(0, n_mod)

      for(j in 1:n_mod){

        if(!is.data.frame(QTL_det_i[[j]])){

          if(QTL_det_i[[j]] == "no_QTL"){ TPR_i[j] <- 0

          } else if (QTL_det_i[[j]] == "error") { TPR_i[j] <- NA

          }

        } else {

          TPR_i[j] <- TPR(QTL.true = QTL_true_r, QTL.detected = QTL_det_i[[j]],
                          distance = d_QTL)

        }

        names(TPR_i) <- mod_names

      }

      TPR_MPP[[i]] <- 100 * TPR_i

    }

    names(TPR_MPP) <- MPP_names

    TPR_res[[r]] <- TPR_MPP

  }

  # Now average and standard deviation over x reps

  TPR_res_mat <- vector(mode = "list", length = n_des)

  for(i in 1:n_des){

    # select all the results from MPPi

    MPP_i_res <- lapply(X = TPR_res, FUN = function(x) x[[i]])

    MPP_i_res <- matrix(unlist(MPP_i_res), ncol = n_mod, byrow = TRUE)

    colnames(MPP_i_res) <- mod_names
    rownames(MPP_i_res) <- paste0("Rep_", 1:n_rep)


    TPR_res_mat[[i]] <- MPP_i_res

  }

  # produce the average per model

  TPR_res_av <- lapply(X = TPR_res_mat, FUN = function(x) colMeans(x, na.rm = TRUE))

  TPR_res_av <- matrix(unlist(TPR_res_av), ncol = n_mod, byrow = TRUE)

  colnames(TPR_res_av) <- mod_names
  rownames(TPR_res_av) <- MPP_names

  # produce the standard deviation per model

  TPR_res_sd <- lapply(X = TPR_res_mat,
                       FUN = function(x) apply(X = x, MARGIN = 2,
                                               FUN = function(x) sd(x = x, na.rm = TRUE)))

  TPR_res_sd <- matrix(unlist(TPR_res_sd), ncol = n_mod, byrow = TRUE)

  colnames(TPR_res_sd) <- mod_names
  rownames(TPR_res_sd) <- MPP_names

  # average power per model over all designs

  TPR_mod_av <- colMeans(TPR_res_av)

  # standard deviation over the design

  sd_des <- lapply(X = TPR_res_mat, FUN = function(x) sd(c(x), na.rm = TRUE))
  sd_des <- unlist(sd_des)
  names(sd_des) <- MPP_names

  # average TPR per design over models.

  TPR_des_av <- rowMeans(TPR_res_av)

  # standard deviation TPR per model over the design

  sd_mod <- rep(0, n_mod)

  for(i in 1:n_mod){

    sd_mod[i] <- sd(unlist(lapply(X = TPR_res_mat,
                                  FUN = function(x, i) x[, i], i=i)), na.rm = TRUE)

  }

  res <- list(av_des = TPR_des_av, sd_des = sd_des, av_mod = TPR_mod_av,
              sd_mod = sd_mod, av_des_mod = TPR_res_av, sd_des_mod = TPR_res_sd)

}

