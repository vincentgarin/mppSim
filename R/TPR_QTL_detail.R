##################
# TPR_QTL_detail #
##################

#' Detail results TPR QTLs
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
#' @param n_par Number of parents of the different designs
#'
#' @param N Total number of individuals in the population
#'
#' @return A data.frame with for each TPR QTL the replication, the design,
#' the number of parents of the design, the number of indidual per cross,
#' the model, the -log10(p.value) and the distance to a simulated QTL,
#' (999 if the QTL was on a chromosome with no simulated QTL). For the
#' run (combination of rep, design and model) where no TPR QTL was present,
#' the -log10pval and the distance to QTL are set to NA.
#'
#' @author Vincent Garin
#'
#' @export

TPR_QTL_detail <- function(QTL_true, QTL_detected, d_QTL, n_des = 9, n_mod = 4,
                           mod_names = c("cr", "par", "anc", "biall"),
                           des_names, n_par, N){

  # list to store the results per QTLs

  n_rep <- length(QTL_detected)

  # Table to store the results

  QTL_TPR_tab <- c()

  col_names <- c("Rep", "MPP_des", "N_par", "Model", "logp", "dQTL")

  for(r in 1:n_rep){

    QTL_true_r <- QTL_true[[r]] # fix the true QTL

    QTL_det_r <- QTL_detected[[r]]

    for(i in 1:n_des){

      QTL_det_i <- QTL_det_r[[i]]

      # iterate over the different models

      # TPR_i <- vector(mode = 'list', length =  n_mod)

      for(j in 1:n_mod){

        if(is.data.frame(QTL_det_i[[j]])){

          QTL_TPR_j <- QTL_TPR(QTL.true = QTL_true_r, QTL.detected = QTL_det_i[[j]],
                               distance = d_QTL)

          if(is.character(QTL_TPR_j)){

            QTL_TPR_tab_j <- data.frame(r, des_names[i], n_par[i], mod_names[j],
                                        NA, NA, stringsAsFactors = FALSE)

            colnames(QTL_TPR_tab_j) <- col_names

            QTL_TPR_tab <- rbind(QTL_TPR_tab, QTL_TPR_tab_j)

          } else {

            QTL_TPR_tab_j <- data.frame(r, des_names[i], n_par[i], mod_names[j],
                                        QTL_TPR_j$log10pval, QTL_TPR_j$d_T_QTL,
                                        stringsAsFactors = FALSE)

            colnames(QTL_TPR_tab_j) <- col_names

            QTL_TPR_tab <- rbind(QTL_TPR_tab, QTL_TPR_tab_j)

          }

        }

      }

    }

  }

  NindCrVar <- mapply(FUN = function(x, y, N) NindCr(N = N, Npar = x, MPP_des = y),
                      x = QTL_TPR_tab$N_par, y = QTL_TPR_tab$MPP_des,
                      MoreArgs = list(N = N))

  ResTable <- data.frame(QTL_TPR_tab, NindCrVar,stringsAsFactors = FALSE)
  colnames(ResTable)[7] <- "N_ind_cr"

  # Factorize to obtain the results in the correct order

  ResTable$Model <- factor(ResTable$Model, levels = c("cr", "par", "anc", "biall"))
  ResTable$MPP_des <- factor(ResTable$MPP_des, levels = c("Diallel", "Chessboard",
                                                          "Factorial", "NAM", "Realized"))

  return(ResTable)

}
