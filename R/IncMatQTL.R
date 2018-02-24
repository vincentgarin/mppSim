#############
# IncMatQTL #
#############

#' Form QTL incidence matrices
#'
#' The function for the QTL incidence matrices of position given different type
#' of effects: cross-specific, parental, ancestral and bi-allelic.
#'
#' @param QTL list of QTL position obtained with function \code{\link{QposSelect}}
#'
#' @param mppData IBD \code{mppData} object
#'
#' @param mppData_bi IBS \code{mppData} object
#'
#' @param par_clu parent clustering object.
#'
#' @return list of QTL incidence matrices
#'
#' @export

IncMatQTL <- function(QTL, mppData, mppData_bi, par_clu){

  Q.pos <- QTL[, 1]
  Q.eff <- QTL[, 4]

  cross.mat <- IncMat_cross(cross.ind = mppData$cross.ind)
  parent.mat <- IncMat_parent(mppData)

  Q.pos <- vapply(X = Q.pos,
                  FUN = function(x, mppData) which(mppData$map[, 1] == x),
                  FUN.VALUE = numeric(1), mppData = mppData)


  Q.list <- mapply(FUN = IncMatQTLMQE, x = Q.pos, Q.eff = Q.eff,
                   MoreArgs = list(mppData = mppData, mppData_bi = mppData_bi,
                                   par.clu = par_clu, cross.mat = cross.mat,
                                   par.mat = parent.mat, order.MAF = TRUE),
                   SIMPLIFY = FALSE)

  return(Q.list)

}
