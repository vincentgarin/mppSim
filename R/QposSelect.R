##############
# QposSelect #
##############

#' Select QTL positions
#'
#' The function select QTL positions that will be used for simulation. First
#' the positions are peaked on different chromosome. Then, the next positions
#' are chosen with a minimum distance in between.
#'
#' @param N Number of positions. Default = 8
#'
#' @param map Three columns map with marker id, chromosome and position in cM.
#'
#' @param dist Minimum distance between the positions on a same chromosome.
#' Default = 30.
#'
#' @return Three columns \code{data.frame} with the selected positions
#'
#' @author Vincent Garin
#'
#' @examples
#'
#' # Later
#'
#' @export

QposSelect <- function(N = 8, map, dist = 30){

  chr.id <- unique(map[, 2])
  n.chr <- length(chr.id)

  int.div <- N %/% n.chr
  rest <- N %% n.chr

  chr.samp <- c(rep(chr.id, int.div), sample(x = chr.id, size = rest,
                                             replace = FALSE))
  chr.samp <- sort(chr.samp)
  chr.id <- unique(chr.samp)
  chr.samp.nb <- table(chr.samp)

  mk.samp <- rep("", N)
  mk.id <- 1

  for(i in 1:length(chr.id)){

    if(chr.samp.nb[i] == 1){

      map_i <- map[map[, 2] == chr.id[i], ]
      mk.samp[mk.id] <- map_i[sample(x = 1:dim(map_i)[1], size = 1), 1]

      mk.id <- mk.id + 1

    } else if (chr.samp.nb[i] > 1) {

      j <- 0

      cand.list <- map[map[, 2] == chr.id[i], c(1, 3)]

      while(j != chr.samp.nb[i]){

        # select a position

        cand.j <- cand.list[sample(x = 1:dim(cand.list)[1], size = 1), 1]

        mk.samp[mk.id] <- cand.j

        # Remove too close positions

        rem.pos <- abs(cand.list[, 2] - cand.list[cand.list[, 1] == cand.j, 2]) < dist
        cand.list <- cand.list[!rem.pos, ]

        j <- j+1
        mk.id <- mk.id + 1

      }

    }

  }

  res <- map[map[, 1] %in% mk.samp, ]
  colnames(res) <- c("mk.id", "chr", "pos.cM")

  return(res)

}
