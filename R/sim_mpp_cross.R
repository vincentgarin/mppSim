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
#' @return List with genotype marker matrices of all simulated crosses in
#' different format:
#'
# geno_num_IBD: numeric IBD: 0, 1, 2: copy of the cross second parent allele.
# geno_chr_IBD: character IBD: "A", "H" "B": within cross indication of parent origin
# geno_chr_IBS: character IBS: "AA", "AC", "CC": marker score
# geno_num_IBS: numeric IBS: 0, 1, 2: number of copies of the minor allele
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
#'                      n_ind_cr = rep(100, nrow(cross_scheme)))
#'
#'
#' @import simcross
#' @import LDcorSV
#' @import ggplot2
#'
#' @export
#'

sim_mpp_cross <- function(geno_par, map, cross_scheme, n_ind_cr = NULL){

  # space to store the results
  geno_num_IBD <- geno_chr_IBD <- geno_chr_IBS <- c()

  n_cross <- nrow(cross_scheme)
  crosses_id <- cross_scheme[, 1]

  if(is.null(n_ind_cr)){n_ind_cr <- rep(100, n_cross)}

  map_pos <- split(x = map[, 3], f = factor(map[, 2]))

  n_chr <- length(map_pos)
  chr_lgh <- unlist(lapply(map_pos, max))

  # iteration over the number of crosses
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

    geno_num_IBD_i <- geno_chr_IBD_i <- geno_chr_IBS_i <- c()

    for(j in 1:n_chr){

      geno_sim_chr_ij <- geno_sim[[j]][5:(n_ind_cr[i]+4), ]

      # numeric IBD
      geno_num_IBD_ij <- geno_sim_chr_ij - 1

      # character IBD
      geno_chr_IBD_ij <- apply(X = geno_num_IBD_ij, MARGIN = 2, as.character)
      geno_chr_IBD_ij[geno_chr_IBD_ij == "0"] <- "A"
      geno_chr_IBD_ij[geno_chr_IBD_ij == "1"] <- "H"
      geno_chr_IBD_ij[geno_chr_IBD_ij == "2"] <- "B"


      # character IBS
      geno_chr_IBS_ij <- geno_sim_chr_ij
      for(k in 1:dim(geno_sim_chr_ij)[2]){

        geno_chr_IBS_ij[, k] <- geno_pipj[as.numeric(geno_chr_IBS_ij[, k]), k]

      }

      # combine the different matrix
      geno_num_IBD_i <- cbind(geno_num_IBD_i, geno_num_IBD_ij)
      geno_chr_IBD_i <- cbind(geno_chr_IBD_i, geno_chr_IBD_ij)
      geno_chr_IBS_i <- cbind(geno_chr_IBS_i, geno_chr_IBS_ij)

    }

    # add labels

    geno_n_id <- paste0("%0", nchar(n_ind_cr[i]),"d")
    geno_id <- paste0("g", sprintf(geno_n_id, 1:n_ind_cr[i]))
    nm_geno_i <- paste0(crosses_id[i], '_', geno_id)

    rownames(geno_num_IBD_i) <- rownames(geno_chr_IBD_i) <-
      rownames(geno_chr_IBS_i) <- nm_geno_i

    # colnames(geno_num_IBD_i) <- colnames(geno_chr_IBD_i) <- map[, 1]

    # add it to the global matrices
    geno_num_IBD <- rbind(geno_num_IBD, geno_num_IBD_i)
    geno_chr_IBD <- rbind(geno_chr_IBD, geno_chr_IBD_i)
    geno_chr_IBS <- rbind(geno_chr_IBS, geno_chr_IBS_i)



  }

  # calculate geno num IBS (minor allele frequency)
  geno_num_IBS <- geno_012(mk.mat = geno_chr_IBS)[[1]]

  colnames(geno_num_IBD) <- colnames(geno_chr_IBD) <- map[, 1]
  colnames(geno_num_IBS) <-colnames(geno_chr_IBS) <- map[, 1]

  return(list(geno_num_IBD = geno_num_IBD, geno_num_IBS = geno_num_IBS,
              geno_chr_IBD = geno_chr_IBD, geno_chr_IBS = geno_chr_IBS))

}
