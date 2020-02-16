###############
# formMppData #
###############

#' Form mppData objects
#'
#' Function that subset a global population composed from all the cross of a
#' Dialleld structure between 9 parents to for mppData objects. The function
#' start from a particular realization (global population) composed of 5000
#' markers.
#'
#' @param geno_ABH genotype matrix of the global population in ABH format.
#'
#' @param geno_012 genotype matrix of the global population in 012 format.
#'
#' @param map Correspondig three column genetic map (mk.id, chr, pos cM).
#'
#' @param MPP_sub subset information obtained as output of the function
#' \code{\link{SampleMPPDesign}}.
#'
#' @param mk_prob Problematic markers in the subset of the global population
#'
#' @param pheno Simulated phenotype values.
#'
#' @param type Type of population. Default = "F"
#'
#' @param nb_gen Number of generation. Default = 2.
#'
#' @return list containing two type of mppData object
#'
#' @author Vincent Garin
#'
#' @export
#'


formMppData <- function(geno_ABH, geno_012, map, MPP_sub, mk_prob, pheno,
                        type = "F", nb_gen = 2){

  geno_ABH_red <- geno_ABH[MPP_sub$sel_ind, ]
  geno012_red <- geno_012[MPP_sub$sel_ind, ]

  sel_par <- paste0("P", MPP_sub$sel_par)
  geno_par_red <- geno_par[sel_par, ]

  if(length(mk_prob) > 0){

    geno_ABH_red <- geno_ABH_red[, -mk_prob]
    geno012_red <- geno012_red[, -mk_prob]
    map <- map[- mk_prob, ]
    geno_par_red <- geno_par_red[, -mk_prob]

  }

  pheno <- pheno[MPP_sub$sel_ind]

  # Various cross and genotype indicators

  crosses <- paste0("c", 1:36)

  # par.per.cross

  P1 <- paste0("P" , c(rep(1, 8), rep(2, 7), rep(3, 6), rep(4, 5), rep(5, 4),
                       rep(6, 3), rep(7, 2), rep(8, 1)))

  P2 <- paste0("P", c(2:9, 3:9, 4:9, 5:9, 6:9, 7:9, 8:9, 9))

  par_per_cross <- cbind(crosses, P1, P2)

  par_per_cross <- par_per_cross[crosses %in% MPP_sub$cr_ind_sel, ]

  cross_ind <- rep(MPP_sub$cr_ind_sel, times = MPP_sub$ind_cr)

  trait <- data.frame(rownames(geno_ABH_red), pheno, stringsAsFactors = FALSE)


  mppData <- mppData_form(geno.off = geno_ABH_red, geno.par = geno_par_red,
                          type = type, F.gen = nb_gen, map = map, trait = trait,
                          cross.ind = cross_ind, par.per.cross = par_per_cross)

  mppData_bi <- mppData_form(geno.off = geno012_red, geno.par = geno_par_red,
                             type = type, F.gen = nb_gen, map = map,
                             trait = trait, IBS = TRUE, IBS.format = "012",
                             cross.ind = cross_ind,
                             par.per.cross = par_per_cross)

  return(list(mppData = mppData, mppData_bi = mppData_bi))


}
