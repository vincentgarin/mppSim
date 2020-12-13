###############
# formMppData #
###############

#' Form mppData object
#'
#' Function that subset a global population composed from all the cross of a
#' Dialleld structure and form mppData object for realisation of particular
#' MPP design.
#'
#' @param geno_sim genotype matrix of the global population.
#'
#' @param nPar Number of parents in the original diallel.
#'
#' @param map Correspondig three column genetic map (mk.id, chr, pos cM).
#'
#' @param MPP_sub subset information obtained as output of the function
#' \code{\link{SampleMPPDesign}}.
#'
#' @param pheno Simulated phenotype values.
#'
#' @param type Type of population.
#'
#' @param F.gen Number of F generations.
#'
#' @param BC.gen Number of BC generations.
#'
#' @param window size of the window used by clusthaplo to cluster the parents.
#'
#' @return mppData object
#'
#' @author Vincent Garin
#'
#' @export
#'


formMppData <- function(geno_sim, nPar, map, MPP_sub, pheno,
                         type = NULL, F.gen = NULL, BC.gen = NULL, window = 20){

  # equalize map and marker matrix

  map <- map[map[, 1] %in% colnames(geno_sim), ]
  geno_sim <- geno_sim[, colnames(geno_sim) %in% map[, 1]]
  geno_sim <- geno_sim[, map[, 1]]

  # prepare mppData object

  geno.off <- geno_sim[(nPar+1):dim(geno_sim)[1], ]
  geno.par <- geno_sim[1:nPar, ]

  geno.off <- geno.off[MPP_sub$sel_ind, ]
  sel_par <- paste0("P", MPP_sub$sel_par)
  geno.par <- geno.par[sel_par, ]

  pheno <- pheno[MPP_sub$sel_ind]
  pheno <- as.matrix(pheno)
  rownames(pheno) <- rownames(geno.off)

  cross_ind <- rep(MPP_sub$cr_ind_sel, times = MPP_sub$ind_cr)

  # ppc <- form_par_per_cross(length(MPP_sub$sel_par))
  ppc <- form_par_per_cross(nPar)

  ppc <- ppc[ppc[, 1] %in% MPP_sub$cr_ind_sel, ]

  mppData <- create.mppData(geno.off = geno.off, geno.par = geno.par,
                            map = map, pheno = pheno, cross.ind = cross_ind,
                            par.per.cross = ppc)

  mppData <- QC.mppData(mppData = mppData, verbose = FALSE)
  mppData <- IBS.mppData(mppData = mppData)
  mppData <- IBD.mppData(mppData = mppData, type = type, F.gen = F.gen,
                         BC.gen = BC.gen)

  mppData <- parent_cluster.mppData(mppData = mppData, window = window,
                                    plot = FALSE, method = 'clusthaplo')


  # par_clu_i <- par_clu[, mppData$parents]
  # par_clu_i <- par_clu_i[mppData$map[, 1], ]

  # mppData <- par_clu_chg.mppData(mppData = mppData, par.clu = par_clu_i)


  return(mppData)

}
