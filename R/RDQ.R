##############################
# Resolution distance to QTL #
##############################

#' Resolution distance to QTL
#'
#' Function measuring the distance between a "true" QTL position and a detected
#' QTL. The function looks if there is a signal at +- distance of the QTL. If a
#' signal was detected in this region, it measures the distance between the QTL
#' position and the detected position.
#'
#' @param QTL.true \code{data.frame} of true QTL positions with marker id,
#' chromosome and position in cM.
#'
#' @param QTL.detected \code{data.frame} of detected QTL positions with marker
#'  id, chromosome and position in cM.
#'
#' @param distance maximal distance between the true QTL position and the detected
#' signal. Default = 10.
#'
#' @return Vector with distance to the position. NA if nothing was detected in
#' the surrounding of the QTL.
#'
#' @author Vincent Garin
#'
#' @export

# QTL.true <- Q_sel[[1]]
# QTL.detected <- Q_res[[1]][[3]]$QTL_par
# distance <- 10

RDQ <- function(QTL.true, QTL.detected, distance = 10) {

  N.QTL <- dim(QTL.true)[1]

  d_QTL <- rep(NA, N.QTL)

    for(i in 1:N.QTL){

      Q_i <- QTL.true[i, ]

      # restrict the positions of the current QTL chromosome

      if((Q_i[, 2] %in% QTL.detected[, 2])){

        Q.det.chr_i <- QTL.detected[QTL.detected[, 2] == Q_i[, 2], ]

        # Select the position with the largest -log10(p-val)

        Q.det.chr_i <- Q.det.chr_i[Q.det.chr_i$log10pval == max(Q.det.chr_i$log10pval), ]

        # check if the distance of the position is in the tolerance distance

        d_QTL_i <- abs(Q.det.chr_i[, 4] - Q_i[, 3])

        if(d_QTL_i < distance) {d_QTL[i] <- d_QTL_i}

      }

    }

  return(d_QTL)

}
