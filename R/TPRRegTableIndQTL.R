#####################
# TPRRegTableIndQTL #
#####################

#' Table of TPR results for single QTLs
#'
#' Organise the TPR results in a table that can be used for regression. The
#' independent variable (y) is the average TPR decomposed over each factor:
#' MPP design, number of parents, individual per cross (continuous variable),
#' QTL detection model, type of QTL, and the size of the QTL.
#'
#' @param QTL_true list of true QTL positions with frequency and number of effect
#' information added.
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
#' @param des_names MPP design names. Use only: 'Diallel', 'Chessboard',
#' 'Factorial', 'NAM' or 'Realized'. Default = NULL.
#'
#' @param n_par Number of parents of the different designs.
#'
#' @param N Total number of individuals in the population.
#'
#' @return A data.frame with first column indicating if the QTL was detected
#' or not (0, 1), the type of QTL effect, the size of the QTL effect,
#' the segregation rate (proportion of ind that receive a non zero effect from
#' the considered QTL), the number of QTL effects, the model, the MPP design,
#' the number of parents, the number of individuals per cross.
#'
#' @author Vincent Garin
#'
#' @export

# load("./results/QTL_detection/An1_N_800/Q_res_1_50.RData")
#
# load("./data/simulation/Pop_36_x_300_100_reps/Q_sim_ext.RData")
#
# Q_res <- Q_res_1_50
# rm(Q_res_1_50)
#
# QTL_true = Q_sim_ext[1:50]
# QTL_detected = Q_res
# d_QTL = 10
# n_des = 9
# n_mod = 4
# n_QTL = 8
# mod_names = c("cr", "par", "anc", "biall")
# des_names <- c("Diallel", "Diallel", "Chessboard", "Chessboard",
#                "Factorial", "Factorial", "NAM", "NAM", "Realized")
# n_par <- c(9, 5, 9, 5, 9, 5, 9, 5, 9)
# N <- 800

TPRRegTableIndQTL <- function(QTL_true, QTL_detected, d_QTL,
                              n_des = 9, n_mod = 4, n_QTL = 8,
                              mod_names = c("cr", "par", "anc", "biall"),
                              des_names = NULL, n_par, N){

  n_rep <- length(QTL_detected)

  MPP_names <- paste(des_names, n_par, sep = "_")

  TPR_res <- vector(mode = "list", length = n_rep)

  for(r in 1:n_rep){

    QTL_true_r <- QTL_true[[r]] # fix the true QTL

    # determine the sub-type of QTL: Q1 small, Q1 large, Q2 small, etc.

    QtypeId <- 1:(2*n_QTL)
    names(QtypeId) <- paste0(rep(1:8, each = 2), rep(c(2, 6), 8))

    # Qtype <- QtypeId[paste0(QTL_true_r$Qeff, QTL_true_r$Qsize)]

    QTL_det_r <- QTL_detected[[r]] # select the detected QTLs

    # create a list to store the results for each MPP design

    TPR_MPP <- vector(mode = "list", length = n_des)

    for(i in 1:n_des){

      QTL_det_i <- QTL_det_r[[i]]

      res <- matrix(0, n_mod, 2*n_QTL)

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

        index <- QtypeId[paste0(QTL_true_rj$Qeff, QTL_true_rj$Qsize)]

        res[, index] <- TPR_j

      }

      rownames(res) <- mod_names
      colnames(res) <- paste(rep(paste0("Q", 1:n_QTL), each = 2),
                             rep(c(2, 6), n_QTL), sep = "_")

      TPR_MPP[[i]] <- res

    }

    names(TPR_MPP) <- MPP_names

    TPR_res[[r]] <- TPR_MPP

  }

  # unfold the results

  ResTable <- c()

  for(i in 1:n_rep){

    sim_conf <- paste(paste0("Q", QTL_true[[i]]$Qeff), QTL_true[[i]]$Qsize,
                      sep = "_")

    col_ind <- paste(rep(paste0("Q", 1:n_QTL), each = 2),
                     rep(c(2, 6), n_QTL), sep = "_")

    sel_col <- col_ind %in% sim_conf

    TPR_res_i <- TPR_res[[i]]
    TPR_res_i <- lapply(X = TPR_res_i, FUN = function(x, ref) x[, ref],
                        ref = sel_col)

    # extract the frequency and the number of alleles for each QTL

    ind_QTL <- match(colnames(TPR_res_i[[1]]), sim_conf)

    fract_Q <- rep(rep(QTL_true[[i]]$fract_n_0[ind_QTL], each = n_mod),
                   times = n_des)
    n_eff <- rep(rep(QTL_true[[i]]$n_alleles[ind_QTL], each = n_mod),
                 times = n_des)
    SizeQ <- rep(rep(QTL_true[[i]]$Qsize[ind_QTL], each = n_mod),
                 times = n_des)
    Qeff <- rep(rep(QTL_true[[i]]$Qeff[ind_QTL], each = n_mod),
                times = n_des)
    Model_var <- rep(rep(mod_names, times = n_QTL), times = n_des)

    Des_var <- rep(des_names, each = n_mod * n_QTL)

    Npar_var <-  rep(n_par, each = n_mod * n_QTL)

    TPR_scores <- unlist(TPR_res_i)

    res_x <- data.frame(TPR_scores, Qeff, SizeQ, fract_Q, n_eff, Model_var,
                        Des_var, Npar_var)

    colnames(res_x) <- c("TPR", "QTL_type", "QTL_size", "seg_rate", "n_eff",
                         "Model", "MPP_des", "N_par")

    ResTable <- rbind(ResTable, res_x)

  }

  NindCrVar <- mapply(FUN = function(x, y, N) NindCr(N = N, Npar = x, MPP_des = y),
                      x = ResTable$N_par, y = ResTable$MPP_des,
                      MoreArgs = list(N = N))

  ResTable <- data.frame(ResTable, NindCrVar, stringsAsFactors = FALSE)
  colnames(ResTable)[9] <- "N_ind_cr"

  # Transform and factorize the results

  Q_ind <- c("Q1_cr", "Q2_cr", "Q3_par", "Q4_par", "Q5_anc", "Q6_anc",
             "Q7_biall", "Q8_biall")

  names(Q_ind) <- as.character(1:n_QTL)
  ResTable$QTL_type <- Q_ind[ResTable$QTL_type]
  ResTable$QTL_type <- factor(ResTable$QTL_type, levels = Q_ind)

  Qs_ind <- c("small", "large")
  names(Qs_ind) <- as.character(c(2, 6))
  ResTable$QTL_size <- Qs_ind[as.character(ResTable$QTL_size)]
  ResTable$QTL_size <- factor(ResTable$QTL_size, levels = Qs_ind)

  ResTable$Model <- factor(ResTable$Model, levels = c("cr", "par", "anc", "biall"))
  ResTable$MPP_des <- factor(ResTable$MPP_des, levels = c("Diallel", "Chessboard",
                                                          "Factorial", "NAM", "Realized"))

  par_ind <- c("5_par", "9_par")
  names(par_ind) <- as.character(c(5, 9))
  ResTable$N_par <- par_ind[as.character(ResTable$N_par)]
  ResTable$N_par <- factor(ResTable$N_par, levels = par_ind)


  return(ResTable)

}
