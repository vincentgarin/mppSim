###########
# dQTLTab #
###########

#' Table with distance to the simulated QTL
#'
#' This function measure the distance between the simulated QTL and the largest
#' significant peak.
#'
#' @param QTL_true list of true QTL positions.
#'
#' @param QTL_detected list results of detected QTLs.
#'
#' @param n_des Number of MPP design. Default = 9.
#'
#' @param n_mod Number of QTL detection models. Default = 4.
#'
#' @param n_QTL Number of simulated QTL. Default = 8.
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
#' @return A matrix with the first column being the average RDQ over the
#' repetitions and the other being the incidence matrices of the factors
#'
#' @author Vincent Garin
#'
#' @export


dQTLTab <- function(QTL_true = Q_sim[1:n_rep], QTL_detected = Q_res,
                    des_names = des_names, n_par = n_par, n_des = n_des,
                    n_mod = 4, n_QTL = 8,
                    mod_names = c('cr', 'par', 'anc', 'biall'), N = N){

  n_rep <- length(QTL_detected)

  ResTab <- c()

  for(r in 1:n_rep){

    QTL_true_r <- QTL_true[[r]]

    QTL_det_r <- QTL_detected[[r]]

    # For each simulated QTL iterate over the design and model from QTL_detected

    for(q in 1:n_QTL){

      QTL_true_rq <- QTL_true_r[q, ]

      # start the iteration over designs and models

      for(i in 1:n_des){

        QTL_det_i <- QTL_det_r[[i]]

        for(j in 1:n_mod){

          if(is.data.frame(QTL_det_i[[j]])){

            # Test if something has been detected on the sim chromosome

            if(QTL_true_rq$chr %in% QTL_det_i[[j]]$chr){

              # select the detected QTL on the chromosome

              Q_sel <- QTL_det_i[[j]][QTL_det_i[[j]]$chr %in% QTL_true_rq$chr, ]
              Q_sel <- Q_sel[which.max(Q_sel$log10pval), ]

              dist_ij <- round(abs(QTL_true_rq$pos.cM - Q_sel$pos.cM), 3)
              Q_type <- paste(QTL_true_rq$Qtype, QTL_true_rq$Qeff, sep = "_")

              res_ij <- c(dist_ij, Q_type, QTL_true_rq$Qsize, des_names[i],
                          n_par[i], mod_names[j])

              ResTab <- rbind(ResTab, res_ij)

            }
          }

        }

      }

     }

  }

  # process the results in a data.frame

  colnames(ResTab) <- c('distance', 'QTL_type', 'QTL_size', 'MPP_des', 'N_par',
                        'Model')
  rownames(ResTab) <- as.character(1:dim(ResTab)[1])

  ResTab <- data.frame(ResTab, stringsAsFactors = FALSE)

  ResTab$distance <- as.numeric(ResTab$distance)

  Q_size <- c('small', 'big')
  names(Q_size) <- c('2', '6')

  Q_size_vect <- Q_size[ResTab$QTL_size]

  ResTab$QTL_size <- Q_size_vect

  ResTab$N_par <- as.numeric(ResTab$N_par)

  NindCrVar <- mapply(FUN = function(x, y, N) NindCr(N = N, Npar = x, MPP_des = y),
                      x = ResTab$N_par, y = ResTab$MPP_des,
                      MoreArgs = list(N = N))

  ResTab <- data.frame(ResTab, NindCrVar, stringsAsFactors = FALSE)
  colnames(ResTab)[7] <- "N_ind_cr"

  # convert into factor for regression

  ResTab$Model <- factor(ResTab$Model, levels = c("cr", "par", "anc", "biall"))
  ResTab$MPP_des <- factor(ResTab$MPP_des, levels = c("Diallel", "Chessboard",
                                                          "Factorial", "NAM", "Realized"))

  ResTab$QTL_type <- factor(ResTab$QTL_type, levels = unique(ResTab$QTL_type))

  ResTab$QTL_size <- factor(ResTab$QTL_size, levels = c("small", "big"))

  return(ResTab)

}
