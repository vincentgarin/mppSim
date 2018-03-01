##########
# TPorFP #
##########

#' TP or FP QTL
#'
#' Given a list of TRUE QTL and a list of detected QTL, the function return
#' for each QTL if it was a true positive or a false positive.
#'
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
#' @return A logical vector: TRUE for TP QTLs and FALSE for FP QTLs
#'
#' @author Vincent Garin
#'
#' @export


TPorFP <- function(QTL.true, QTL.detected, distance = 5) {

  N.QTL <- dim(QTL.detected)[1]

  if(distance == 0){

    test <- QTL.detected[, 1] %in% QTL.true[, 1]

  } else {

    test <- rep(FALSE, N.QTL)

    for(i in 1:N.QTL){

      Q_i <- QTL.detected[i, ]

      # restrict the positions of the current QTL chromosome

      if((Q_i[, 2] %in% QTL.true[, 2])){

        Q.t.chr_i <- QTL.true[QTL.true[, 2] == Q_i[, 2], ]
        test[i] <- any(abs(Q.t.chr_i[, 3] - Q_i[, 4]) <= distance)

      }

    }

  }

  return(test)

}
