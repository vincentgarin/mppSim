###########
# QTL_FDR #
###########

#' Details QTL FDR
#'
#' Function that return the list of FDR QTL with the distance to the nearest
#' simlated QTL. If no QTL was simulated on the chromosome return 999.
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

QTL_FDR <- function(QTL.true, QTL.detected, distance = 5){

  N.QTL <- dim(QTL.detected)[1]

  test <- rep(TRUE, N.QTL)

  for(i in 1:N.QTL){

    Q_i <- QTL.detected[i, ]

    # restrict the positions of the current QTL chromosome

    if((Q_i[, 2] %in% QTL.true[, 2])){

      Q.t.pos_i <- QTL.true[QTL.true[, 2] == Q_i[, 2], 3]
      test[i] <- sum(abs(Q.t.pos_i - Q_i[, 4]) <= distance) == 0

    }

  }

  # subsets the undetected QTLs

  if(sum(test) == 0){QTL_FDR <- 'no_FDR'} else {

    QTL_FDR <- QTL.detected[test, ]

    d_T_QTL <- rep(0, dim(QTL_FDR)[1])

    for(j in 1:dim(QTL_FDR)[1]){

      chr_j <- QTL_FDR[j, 2]

      if(!(chr_j %in% QTL.true$chr)){

        d_T_QTL[j] <- 999

      } else {

        d_T_QTL[j] <- abs(QTL.true[QTL.true$chr == chr_j, 3] - QTL_FDR[j, 4])

      }

    }

    QTL_FDR <- data.frame(QTL_FDR, d_T_QTL, stringsAsFactors = FALSE)

  }

  return(QTL_FDR)

}
