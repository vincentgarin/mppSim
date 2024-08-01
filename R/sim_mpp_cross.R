#################
# sim_mpp_cross #
#################

#' Simulate MPP crosses population
#'
#' The function simulates a multi-parent population composed of N F2 crosses.
#' It uses the genotypes data provided in the argument geno_par. The cross
#' scheme is specified in the argument cross_scheme. The functions call
#' function from the simcross package (Broman, 2019)
#'
#' @param geno_par Character genotype marker matrix of the parents used for the simulation.
#' The allele must be coded with one letter per allele. For example, AA, CC, GG,
#' etc. The rownames must specify the parent identifiers used in the crossing
#' scheme (cross_scheme).
#'
#' @param map Data.frame. Correspondig three columns genetic map
#' (marker, chr, pos in cM).
#'
#' @param cross_scheme Three column marker matrix specifying: the cross
#' indicator, the first parent (male), and the second parent (female).
#' The parents identifiers must be the same as the rownames of geno_par.
#'
#' @param n_ind_cr numeric vector indicating the number of simulated genotype
#' per cross. by default 100 genotypes per cross.
#'
#' @param geno_score_format \code{Character} string specifying the format of the
#' marker scores. "numeric_IBD": 0, 1, 2 (parent 1, heterozygous, parent 2 within cross meaning)
#' "parents_IBS": translate the (1, 2, 3) IBD information into
#' IBS parent marker scores. Default = "numeric_IBD".
#'
#' @return Genotype marker matrix of all simulated crosses
#'
#' @author Vincent Garin
#'
#' @references
#'
#' Karl W. Broman (2019) R/simcross: an R package for simulating and plotting general
#' experimental crosses. \url{https://github.com/kbroman/simcross}
#'
#' @examples
#'
#' data("geno_par")
#' rownames(geno_par) <- paste0('P', 1:nrow(geno_par))
#'
#' data("EUNAM_map")
#' rownames(map) <- map[, 1]
#'
#' # NAM crossing scheme with 9 parents
#' cross_scheme <- cross_scheme_NAM(9)
#'
#' geno <- sim_mpp_cross(geno_par = geno_par, map = map,
#'                      cross_scheme = cross_scheme,
#'                      n_ind_cr = rep(100, nrow(cross_scheme)),
#'                      geno_score_format = "numeric_IBD")
#'
#'
#' @import simcross
#' @import LDcorSV
#' @import ggplot2
#'
#' @export
#'

sim_mpp_cross <- function(geno_par, map, cross_scheme, n_ind_cr = NULL,
                          geno_score_format = "parents_IBS"){

  # iteration over the number of crosses

  pop_res <- c()

  n_cross <- nrow(cross_scheme)
  crosses_id <- cross_scheme[, 1]

  if(is.null(n_ind_cr)){n_ind_cr <- rep(100, n_cross)}

  map_pos <- split(x = map[, 3], f = factor(map[, 2]))

  n_chr <- length(map_pos)
  chr_lgh <- unlist(lapply(map_pos, max))

  for(i in 1:n_cross){

    # select the two parents scores
    cr_par <- cross_scheme[i, 2:3]
    geno_p1p2 <- geno_par[cr_par, ]

    # determine the heterozygous

    het_score <- apply(X = geno_p1p2, MARGIN = 2, FUN = het_det)
    geno_pipj <- rbind(geno_p1p2[1, ], het_score, geno_p1p2[2, ])

    # simulate the cross genotypes

    f2 <- sim_ril_pedigree(ngen = 1, selfing = TRUE, parents = 1:2)
    ext_ind <- data.frame(id = 5:(n_ind_cr[i]+4), mom=3, dad=3, sex=0,gen=2)
    f2_ext <- rbind(f2, ext_ind)

    xodat_ext <- sim_from_pedigree(f2_ext, L=chr_lgh, xchr = rep(FALSE, n_chr))
    geno_sim <- convert2geno(xodat_ext, map = map_pos)

    # fill with the actual parent genotypes scores

    geno_res_i <- c()

    for(j in 1:n_chr){

      geno_sim_chr_j <- geno_sim[[j]][5:(n_ind_cr[i]+4), ]

      if(geno_score_format == "parents_IBS"){

        for(k in 1:dim(geno_sim_chr_j)[2]){

          geno_sim_chr_j[, k] <- geno_pipj[as.numeric(geno_sim_chr_j[, k]), k]

        }

      } else if (geno_score_format == "numeric_IBD"){

        geno_sim_chr_j <- geno_sim_chr_j - 1

      }

      geno_res_i <- cbind(geno_res_i, geno_sim_chr_j)

    }

    # add labels

    geno_n_id <- paste0("%0", nchar(n_ind_cr[i]),"d")
    geno_id <- paste0("g", sprintf(geno_n_id, 1:n_ind_cr[i]))
    nm_geno_i <- paste0(crosses_id[i], '_', geno_id)

    rownames(geno_res_i) <- nm_geno_i

    pop_res <- rbind(pop_res, geno_res_i)

  }

  colnames(pop_res) <- map[, 1]
  return(pop_res)

}
