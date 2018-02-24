#############
# QCMarkers #
#############

#' Quality control on markers
#'
#' Perform a quick QC on marker to detect the one that do not segregate in at
#' least on of the cross with 25%. The monomorphic markers are therefore also
#' tagged.
#'
#' @param geno genotype marker matrix
#'
#' @param cross_ind Cross indicator
#'
#' @return Numeric vector with problematic markers.
#'
#' @author Vincent Garin
#'
#' @export

QCMarkers <- function(geno, cross_ind, parallel = FALSE, cluster = NULL){

  off.MAF <- QC_MAF(mk.mat = geno, cross.ind = cross_ind, parallel = parallel,
                    cluster = cluster)

  # segregation within crosses

  MAF.cr <- off.MAF[[2]]

  lim <- rep(.25, length(unique(cross_ind)))
  mk.sel <- rep(FALSE, dim(geno)[2])

  for(i in 1:dim(geno)[2]){

    test <- MAF.cr[, i] > lim
    test[is.na(test)] <- FALSE

    if(sum(test) > 0){mk.sel[i] <- FALSE } else {mk.sel[i] <- TRUE }

  }

  return(c(which(mk.sel)))

}
