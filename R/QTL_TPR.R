###########
# QTL_TPR #
###########

#' Details QTL TPR
#'
#' Function that return the list of TPR QTL with the -log10(p-value) and the
#' distance to the simulated QTL.
#'
#' @param QTL.true \code{data.frame} of true QTL positions with marker id,
#' chromosome and position in cM.
#'
#' @param QTL.detected \code{data.frame} of detected QTL positions with marker
#'  id, chromosome and position in cM.
#'
#' @param distance maximal distance between the true QTL position and the detected
#' signal. Default = 5.
#'
#' @return FDR
#'
#' @author Vincent Garin
#'
#' @export

# QTL.true = QTL_true_r
# QTL.detected = QTL_det_i[[j]]
# distance = d_QTL

QTL_TPR <- function(QTL.true, QTL.detected, distance = 5){

  N.QTL <- dim(QTL.detected)[1]

  test <- rep(FALSE, N.QTL)

  for(i in 1:N.QTL){

    Q_i <- QTL.detected[i, ]

    # restrict the positions of the current QTL chromosome

    if(Q_i[, 2] %in% QTL.true[, 2]){

      Q.t.pos_i <- QTL.true[QTL.true[, 2] == Q_i[, 2], 3]
      test[i] <- sum(abs(Q.t.pos_i - Q_i[, 4]) <= distance) == 1

    }

  }

  # subsets the undetected QTLs

  if(sum(test) == 0){QTL_TPR <- 'no_TP'} else {

    QTL_TPR <- QTL.detected[test, ]

    d_T_QTL <- rep(0, dim(QTL_TPR)[1])

    for(j in 1:dim(QTL_TPR)[1]){

      chr_j <- QTL_TPR[j, 2]

        d_T_QTL[j] <- abs(QTL.true[QTL.true$chr == chr_j, 3] - QTL_TPR[j, 4])

    }

    QTL_TPR <- data.frame(QTL_TPR, d_T_QTL, stringsAsFactors = FALSE)

  }

  return(QTL_TPR)

}
