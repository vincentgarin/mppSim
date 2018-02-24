########################
# False discovery rate #
########################

#' False discovery rate
#'
#' Function computing the false discovery rate (FDR) as the percentage of detected
#' QTL positions that are not in a neighbouring region (+-distance) from a true
#' QTL position.
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

FDR <- function(QTL.true, QTL.detected, distance = 5){

  N.QTL <- dim(QTL.detected)[1]

  if(distance == 0){

    n.undetect <- sum(!(QTL.detected[, 1] %in% QTL.true[, 1]))

  } else {

    test <- rep(1, N.QTL)

    for(i in 1:N.QTL){

      Q_i <- QTL.detected[i, ]

      # restrict the positions of the current QTL chromosome

      if((Q_i[, 2] %in% QTL.true[, 2])){

        Q.t.chr_i <- QTL.true[QTL.true[, 2] == Q_i[, 2], ]
        test[i] <- (sum(abs(Q.t.chr_i[, 3] - Q_i[, 4]) <= distance) != 1) * 1

      }

    }

    n.undetect <- sum(test)

  }

  FDR <- n.undetect/N.QTL

  return(FDR)

}
