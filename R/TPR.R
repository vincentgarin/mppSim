######################
# True positive rate #
######################

#' True positive rate
#'
#' Function computing the true positive rate (TPR) as the percentage of true
#' QTL positions that have been detected at a certain distance
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
#' @return TPR
#'
#' @author Vincent Garin
#'
#' @export

TPR <- function(QTL.true, QTL.detected, distance = 5) {

  N.QTL <- dim(QTL.true)[1]

  if(distance == 0){

    n.detect <- sum(QTL.true[, 1] %in% QTL.detected[, 1])

  } else {

    test <- rep(0, N.QTL)

    for(i in 1:N.QTL){

      Q_i <- QTL.true[i, ]

      # restrict the positions of the current QTL chromosome

      if((Q_i[, 2] %in% QTL.detected[, 2])){

        Q.det.chr_i <- QTL.detected[QTL.detected[, 2] == Q_i[, 2], ]
        test[i] <- any(abs(Q.det.chr_i[, 4] - Q_i[, 3]) <= distance) * 1
        # test[i] <- (sum(abs(Q.det.chr_i[, 4] - Q_i[, 3]) <= distance) == 1) * 1

      }

    }

    n.detect <- sum(test)

  }

  TPR <- n.detect/N.QTL

  return(TPR)

}
